#-------------------------------------------------------------------------------
# Name:        orthogonality.py
# Purpose:     Determine the orthogonality of building footprints
#
# Author:      Adebayo .Y. Ishola
#
# Created:     February 2020
#
# Note:        Buildings have to be closed and not with rings within. (Simple)
#              [Rectangular buildings]
#-------------------------------------------------------------------------------

'''
 Data:        https://www.geofabrik.de/data/shapefiles.html (OSM)

 Note:        Ensure buildings vector dataset is projected to UTM
'''



#packages

import os, arcpy, math, matplotlib, numpy as np, warnings, time

warnings.simplefilter(action='ignore', category=FutureWarning)                                                          # ignore FutureWarning from Numpy due to version installed with ArcGIS 10.7

# validate extension of feature data :: [Adapted for class exercise -- Dr. Nick]
def controlExtension(inName, ext):                                                                                      # checking for feature class input file type to ensure it end with right format
    if inName.rfind('.') > 0:
        return inName[:inName.find('.')] + ext                                                                          # return file name as well as extension
    else:
        return inName + ext                                                                                             # if not, still return file name and extension

# validate file existence in path :: [Adapted for class exercise -- Dr. Nick]
def checkExistence(pathList):                                                                                           # function to check for existence of input feature file
    check = True
    for data in pathList:                                                                                               # looping through the path list to validate if input FeatureClass exists
        if not arcpy.Exists(data):
            check = False
            print '! dataset ' + data + ' is missing'                                                                   # print this statement if data set can't be found in the input folder
            break
    return check

# extracts the field names :: [Adapted for class exercise -- Dr. Nick]
def getFieldNames(table):                                                                                               # function reads attribute table of FeatureClass and gets fieldNames
    fnames = []
    fields = arcpy.ListFields(table)                                                                                    # arcpy function that handles listing of fields in attribute table
    if fields:
        for field in fields:                                                                                            # loop through fields and append to empty array of fnames declared
            fnames.append(field.name)
    return fnames

# create subdirectory to save outputs :: [Adapted for class exercise -- Dr. Nick]
def createSubdir(workspace, subdirList):                                                                                # function takes workspace and list of subdirectories to be created
    for subdir in subdirList:                                                                                           # loop through list of subdirectories to be created and join to workspace
        if not os.path.isdir(workspace + '/' + subdir):                                                                 # append path with / for path completion
            os.mkdir(os.path.join(workspace, subdir))

# join project directory, subdirs and file names into complete paths :: [Adapted for class exercise -- Dr. Nick]
def completePath(workspace, subdir, nameList):                                                                          # function to create a complete path with workspace, subdirectory and name of files
    for ix in range(len(nameList)):
        nameList[ix] = workspace + '/' + subdir + '/' + nameList[ix]                                                    # create a list of the output files with a complete path
    return nameList

# creates output files
def outputFiles(outputList):                                                                                            # function that handles the creation of geoprocessing output files
    file = []
    for output in outputList:                                                                                           # loop through output list and append .shp extension to each file
        file.append(str(output) + '.shp')
    return file

# Generalize/simplify polygon -- remove unnecessary polygon vertices -- Douglas-Peucker Algorithm
def simplifyBuilding(inFC, outFC):                                                                                      # function to simplify building polygons by removing vertices at a tolerance of 1meter
    arcpy.SimplifyPolygon_cartography(inFC, outFC, algorithm="POINT_REMOVE", tolerance=1, minimum_area=0)

# convert polygon to line
def polygonToLine(inFC, outFC):                                                                                         # function to convert polygon/building to continous polyline using default paramaters
    arcpy.PolygonToLine_management(inFC, outFC)

# split polygon polylines at vertices
def polylineToSegments(inFC, outFC):                                                                                    # function to convert building polyline to line segments 
    arcpy.SplitLine_management(inFC, outFC)

# get building coordinates in arrays
def _building_arr_(geometry):                                                                                           # function to get building vertex coordinates and store in an array
    """Return building coordinates."""
    def _split_(part):
        yield [(p.X, p.Y) for p in part if p]
        
    coords =[]                                                                                                          # closed array of coordinates (first vertex squeezed at end of array).
    for part in geometry:                                                                                               # loop through the geometry parts and push coordinates in an array using numpy
        output = []
        w = np.where(np.in1d(part, None, invert=False))[0]
        bits = np.split(part, w)
        for bit in bits:
            geom_sub = _split_(bit)
            output.append(np.array(*geom_sub).squeeze())
        coords.append(np.asarray(output).squeeze())
    return np.asarray(coords)                                                                                           # return numpy array of coordinates

# get angles of polygons or polylines (polylines here should not be split but rather continuous)
def _building_angles_(a, inside=True, in_degrees=True):                                                                 # function to compute building vertex angles

    def _xy_(a):
        ba = a - np.concatenate((a[-1, None], a[:-1]), axis=0)
        bc = a - np.concatenate((a[1:], a[0, None]), axis=0)
        return np.cross(ba, bc), ba, bc                                                                                 # cross product for vector computation

    if np.allclose(a[0], a[-1]):                                                                                        # closed loop, remove duplicates
        a = a[:-1]                                                                                                      # vertex coordinates
    cr, ba, bc = _xy_(a)
    dt = np.einsum('ij,ij->i', ba, bc)
    angle = np.arctan2(cr, dt)                                                                                          # angle here in radian
    
    if inside:                                                                                                          # with "inside" True is returned for computation of interior angles.
        angles = np.where(angle < 0, angle + (np.pi * 2.), angle)
    else:
        angles = np.where(angle > 0, (np.pi * 2.) - angle, angle)
    if in_degrees:                                                                                                      # with "in_degrees" True for degrees while False for radians angle computation
        angles = np.degrees(angles)                                                                                     # radians angle to degrees computation
    return angles

# get polygon angles
def geomAngles(inFC):                                                                                                   # function to compute actual angles of geometry
    with arcpy.da.SearchCursor(inFC, "SHAPE@") as cur:
        angles = []                                                                                                     # store all angles in array
        for row in cur:
            geom = row[0]                                                                                               # field index 0 contains shape geometry
            arr = _building_arr_(geom).squeeze()                                                                        # adds an extra list
            ang = _building_angles_(arr, inside=True, in_degrees=True).tolist()                                         # to list to convert to standard python list
            
            angles.append(ang)
        return angles

# get building angles at each vertex
def buildingVertexAngle(table, angles):                                                                                 # function to assign building interior angles to a field

    arcpy.AddField_management(table,"Angle","DOUBLE")                                                                   # create a new field called "Angle" and set to data type "DOUBLE"

    pointer = 0                                                                                                         # pointer start at "0" to update each row iteratively

    with arcpy.da.UpdateCursor(table, 'Angle') as uCur:                                                                 # use update cursor to add angle value to each row

        # store all angles in a list
        orient = []                                                                                                     # array to store angles
        for angle in angles:                                                                                            # looping through for each angle
            for i in angle:
                orient.append(i)
            
        for row in uCur:                                                                                                # looping through each row
            row[0] = orient[pointer]

            pointer += 1
            uCur.updateRow(row)                                                                                         # update each row

        return orient

# get building angles error (!= 90 degrees) at each vertex
def buildingVertexError(table, angles):                                                                                 # function to compute difference in vertex angle not orthogonal
    arcpy.AddField_management(table,"Angle_Err","DOUBLE")

    pointer = 0

    with arcpy.da.UpdateCursor(table, 'Angle_Err') as uCur:                                                             # writing of output vertex angle error via update cursor

        # store all angles in a list
        angleError = []                                                                                                 # list to store angle error
        for angle in angles:
            error = 90 - angle
            angleError.append( error )
            
        for row in uCur:
            row[0] = angleError[pointer]

            pointer += 1
            uCur.updateRow(row)                                                                                         # update by row the difference to right angle (90 degrees)

        return angleError

# split polyline length computation algorithm
def _length_(geom):                                                                                                     # function to compute length
    def _calc_(b):
        diff = b[:-1] - b[1:]                                                                                           # compute difference in length using coordinates of line segments
        return np.sqrt(np.einsum('ij,ij->i', diff, diff))
    return _calc_(geom)

# get length of split polyline
def geomLength(inFC):                                                                                                   # function to compute length of geometry
    with arcpy.da.SearchCursor(inFC, "SHAPE@") as cur:                                                                  # searching for each segment to compute length
        _lengths = []
        _arrs = []
        for row in cur:
            geom = row[0]
            arr = _building_arr_(geom).squeeze()                                                                        # squeeze to remove external appending brackets
            length = _length_(arr)

            _lengths.append(length)                                                                                     # append length to array of lengths
            _arrs.append(arr)

        return _lengths

# get building edge/segment length
def buildingLengths(table, lengths):                                                                                    # function to compute building length and store in table

    arcpy.AddField_management(table,"Length","DOUBLE")                                                                  # add a field called "Length" of type "Double"

    pointer = 0                                                                                                         # pointer to go through fields iteratively

    with arcpy.da.UpdateCursor(table, 'Length') as uCur:                                                                # updating Length field

        # store all lengths in a list
        buildLengths = []                                                                                               # store building lengths in a list                                   
        for length in lengths:
            for i in length:
                buildLengths.append(i)
            
        for row in uCur:                                                                                                # update length field row by row
            row[0] = buildLengths[pointer]

            pointer += 1
            uCur.updateRow(row)                                                                                         # update cursor to update length field by row
        return buildLengths

# building segments coordinates
def geomCoords(inFC):                                                                                                   # function to store coordinates of each segment for later reuse
    fieldList = ["SHAPE@", "RIGHT_FID", "Angle", "Length"]
    with arcpy.da.SearchCursor(inFC, fieldList) as cur:                                                                 # use a cursor to search each segment

        Coords = []
        for row in cur:
            geom = row[0]
            arr = _building_arr_(geom).squeeze()
            Coords.append( arr )
        return Coords                                                                                                   # returns an array of the coordinates of each line segment

# building edge longest side
def geoMaxLength(inFC):                                                                                                 # function to compute maximum length of each segment to be used as reference side
    coords = geomCoords(inFC)                                                                                           # using function defined above to get coordinates of building line segments
    groups = []                                                                                                         # store each maximum building segment and reuse to get it index location 
    buildings = dict()                                                                                                  # dictionary to store building info as key-value pair
    items = 0

    for i, coord in enumerate(coords):
        groups.append(coord)

        if i % 4 == 0:
            groups = []
            groups.append(coords[i])

        if len(groups) == 4:
            lengths = []
            for j in groups:
                length = _length_(j)
                lengths.append(length)
            lengths = np.array(lengths)

            maxlen = np.max(lengths)                                                                                    # get segment with maximum length 
            loc = np.argmax(lengths)                                                                                    # get index location of that segment with maximum length

            buildings[items] = dict()
            buildings[items]['max_coord'] = groups[loc]
            buildings[items]['max_len'] = maxlen
            buildings[items]['max_index'] = loc
        
            items = items + 1

    return buildings                                                                                                    # return the buildings

# building edge shortest side
def geoMinLength(inFC):                                                                                                 # function to compute minimum length of each segment to be used as reference side
    coords = geomCoords(inFC)                                                                                           # using function defined previously to get coordinates of building line segments    
    groups = []                                                                                                         # store each minimum building segment and reuse to get it index location 
    buildings = dict()                                                                                                  # dictionary to store building info as key-value pair
    items = 0

    for i, coord in enumerate(coords):
        groups.append(coord)

        if i % 4 == 0:
            groups = []
            groups.append(coords[i])

        if len(groups) == 4:
            lengths = []
            for i in groups:
                length = _length_(i)
                lengths.append(length)
            lengths = np.array(lengths)

            minlen = np.min(lengths)                                                                                    # get segment with minimum length 
            loc = np.argmin(lengths)                                                                                    # get index location of that segment with minimum length

            buildings[items] = dict()
            buildings[items]['min_coord'] = groups[loc]
            buildings[items]['min_len'] = minlen
            buildings[items]['min_index'] = loc
        
            items = items + 1

    return buildings

def calculateNewCoord(a, b, alpha, min_length):                                                                         # function to compute new vertex of buildings
    aa = abs(min_length * (np.cos(alpha + np.pi/2.)) - a)
    bb = abs(min_length * (np.sin(alpha + np.pi/2.)) - b)

    return a, b, aa, bb

# bring it all together: maximum length, minimum length, their coordinates and index of building segments
def bringAllTogether(inFC):                                                                                            # function to bring new vertex positions of buildings together and in ther order

    geo_sides = 4                                                                                                       # fixed building geometry constraint of 4 vertices, that is, rectangle
    max_sides = geoMaxLength(inFC)
    min_sides = geoMinLength(inFC)

    all_coords = []

    for a in max_sides:
        max_coord, max_index, max_length = max_sides.get(a)['max_coord'], max_sides.get(a)['max_index'], max_sides.get(a)['max_len']
        min_coord, min_index, min_length = min_sides.get(a)['min_coord'], min_sides.get(a)['min_index'], min_sides.get(a)['min_len']

        x1, y1 = max_coord[0]
        x2, y2 = max_coord[1]
        xdiff = x2 - x1
        ydiff = y2 - y1
        alpha = np.arctan2(ydiff, xdiff)

        shape_order = None                                                                                              # variable to store order of building segments

        new_coords = [[] for i in range(geo_sides)]
        for i in range(geo_sides):                                                                                      # check to see on what index is the new vertex coordinates computed located
            temp_coord = []
            if i == 0:
                a, b, aa, bb =   calculateNewCoord(x1, y1, alpha, min_length)
                new_coords[i] = np.array([[a, b], [aa, bb]]).tolist()

            if i == max_index:
                new_coords[i] = max_coord.tolist()

            if i == 2:
                a, b, aa, bb =   calculateNewCoord(x2, y2, alpha, min_length)
                new_coords[i] = np.array([[a, b], [aa, bb]]).tolist()

        if max_index == 1:                                                                                              # check to see if maximum length falls on segments coordinates with index == 1
            new_coords[3] = np.array([new_coords[0][1], new_coords[2][1]]).tolist()
            shape_order = 1
            
        elif max_index == 3:                                                                                            # check to see if maximum length falls on segments coordinates with index == 1
            new_coords[1] = np.array([new_coords[0][1], new_coords[2][1]]).tolist()
            shape_order = 3


        ordered_coords = dict()
        ordered_coords['coords'] = new_coords
        ordered_coords['order'] = shape_order
        all_coords.append(ordered_coords)

    return all_coords                                                                                                   # all coordinates of the buildings

# create an empty shapefile to store resulting orthogonal buidlings
def createEmptyShapefile(output_path, fc_name, spatReference):
    arcpy.CreateFeatureclass_management(output_path, fc_name, "POLYGON", "",
                                        "DISABLED", "DISABLED",
                                        spatReference)

# coupling of fixed building segments back into a polygon (rectangle)
def fixedSegments(coords, outFC):                                                                                       # function to form new buildings with true orthogonal segments
    arcpy.AddField_management(outFC,"RIGHT_ID","LONG")                                                                  # add a field "RIGHT_ID" with data type "LONG" to store building polygon index
    
    with arcpy.da.InsertCursor(outFC, ["SHAPE@", "RIGHT_ID"]) as cursor:
        for id,coord in enumerate( range( len(coords) ) ):

            first = coords[id]['coords'][0][1]
            
            if coords[id]['order'] == 1:
                second = coords[id]['coords'][0][0]
                third = coords[id]['coords'][1][1]
            
            elif coords[id]['order'] == 3:
                second = coords[id]['coords'][0][0]
                third = coords[id]['coords'][2][0]
            
            fourth = coords[id]['coords'][2][1]
            

            try:                                                                                                        # try catch to check for building geometries not considered in the algorith (i.e > 4)                        
                array = arcpy.Array([ arcpy.Point(first[0], first[1]),                                                  # arrays of building new coordinates to form a polygon of 4 sides
                                    arcpy.Point(second[0], second[1]),
                                    arcpy.Point(third[0], third[1]), 
                                    arcpy.Point(fourth[0], fourth[1]) ])
                
                polygon = arcpy.Polygon(array)                                                                          # polygon formation

                cursor.insertRow([polygon,id])                                                                          # insert new building "SHAPE" geometry row by row
            except UnboundLocalError as e:
                pass
            



#******************** FUNCTION CALL SECTION *******************************
# output
arcpy.env.overwriteOutput = True

# get file current directory -- change forward slash to backslash
currDirectory = ( os.path.dirname(os.path.realpath(__file__)) ).replace(os.sep, '/')

# TODO: create output directory programmatically and see how to handle input directory also
subdirList = ['temp', 'output']

# set current directory as workspace
workspace = currDirectory

# create subdirectory
createSubdir(workspace, subdirList)

# input
# TODO: allow users to enter input feature in "input" folder
inFCName = 'test_buildings.shp'

# validate file extension
inFCName = controlExtension(inFCName, '.shp')

# get complete path of input feature
inFC = completePath(workspace, 'input', [inFCName])[0]

# print "FieldNames of inputFC: ", getFieldNames(inFC)

#=========== processing files output list
outputFCNames = ['simplifyBuildings', 'polylines', 'splitPolyline', 'fixedBuildings']

# output function call
outputFile = outputFiles(outputFCNames)


simplifyBuildingsFC = completePath(workspace, 'temp', [ outputFile[0] ] )[0]
polyLinesFC = completePath(workspace, 'temp', [ outputFile[1] ] )[0]
splitLinesFC = completePath(workspace, 'temp', [ outputFile[2] ] )[0]
fixedBuildingsFC = completePath(workspace, 'output', [ outputFile[3] ] )[0]

fixedBuildingsFC = "fixedBuildings.shp"



#************** FUNCTION CALLS

start_time = time.time()                                                                                    # function to compute time of python file execution

print "\t"
arcpy.AddMessage("******** CLOSE ALL FILES EXCEPT INPUT IN ARCMAP BEFORE RUNNING SCRIPT ********")

print "\t"
check = checkExistence([inFC])
print "Is input FeatureClass available?: ", check, "\n", "InFC: \t \t \t \t  ", inFC
print "\t"
arcpy.AddMessage(">> Input FeatureClass check DONE!")

print "\n"

simplifyBuilding(inFC, simplifyBuildingsFC)
arcpy.AddMessage(">> Building polygon simplification -- Douglas-Peucker Algorithm -- DONE!")

print "\n"

polygonToLine(simplifyBuildingsFC, polyLinesFC)
arcpy.AddMessage(">> Polygon to polyline DONE!")

print "\n"

polylineToSegments(polyLinesFC, splitLinesFC)
arcpy.AddMessage(">> Split polyline to segments DONE!")

print "\n"

angles = geomAngles(polyLinesFC)
arcpy.AddMessage(">> Get building geometry angles DONE!")

print "\n"

measured_angle = buildingVertexAngle(splitLinesFC, angles)
arcpy.AddMessage(">> Assign building angles DONE!")

print "\n"

buildingVertexError(splitLinesFC, measured_angle)
arcpy.AddMessage(">> Vertex building error angles DONE!")

print "\n"

lengths = geomLength(splitLinesFC)
arcpy.AddMessage(">> Assign building lengths DONE!")

print "\n"

buildingLengths(splitLinesFC, lengths)
arcpy.AddMessage(">> Get building lengths DONE!")

print "\n"

geomCoords(splitLinesFC)
arcpy.AddMessage(">> Get building coordinates DONE!")

print "\n"

geoMaxLength(splitLinesFC)
arcpy.AddMessage(">> Get maximum lengths DONE!")

print "\n"

geoMinLength(splitLinesFC)
arcpy.AddMessage(">> Get minimum lengths DONE!")

print "\n"

spatRef = arcpy.Describe(inFC).spatialReference
arcpy.AddMessage(">> Spatial Reference acquisition DONE!")

print "\n"

createEmptyShapefile(workspace + "/output", fixedBuildingsFC, spatRef)
arcpy.AddMessage(">> Empty shapefile creation DONE!")

print "\n"

fixedCoord = bringAllTogether(splitLinesFC)
arcpy.AddMessage(">> Bring it all together DONE!")

print "\n"

fixedSegments(fixedCoord, workspace + "/output/" + fixedBuildingsFC)
arcpy.AddMessage(">> Creating fixed segments DONE!")
print "Output is here: \t \t \t ", workspace + "/output/" + fixedBuildingsFC

print "\n"


X = str(time.time() - start_time)
print "-"*40
print "Time of Execution: ", X, "secs"
print "-"*40

print "\t"