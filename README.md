# ORTHOGONALITY

## Enforcing Orthogonality of Segments in Building Footprints
- - - - 

## TASK DESCRIPTION: ##
Create new building footprints dataset that enforces rectangular corners, wherever the original ones fall in between 90°-ε and 90°+ε. ε is an inaccuracy measure.

## PROBLEM STATEMENT: ##

Often, after buildings have been extracted from LiDAR images, they produce building outlines that are not correctly orthogonal at its corners whereas they are in reality. It is then important to regularize these outlines for further use for mapmaking, further analysis and eventually, decision making.


The algorithm developed has been modularized using functional programming approaches within the context of valid (building) polygons as input.
A valid polygon in this case, has unambiguous exteriors, no touching segments except at vertices, non-zero length, clockwise direction of outer rings, non-zero area and with no overlap.

A valid polygon is said to be simple.

![polygon](https://github.com/bayoishola20/buildings-orthogonality/blob/master/assets/polygons.png "Building geometry")

![polygon](https://github.com/bayoishola20/buildings-orthogonality/blob/master/assets/polygons_.png "Building geometry")

**CONSTRAINT DEFINITION**: Four-cornered buildings with no inner rings and no self-intersections. Only rectangle-shaped buildings considered.

## WORKFLOW

![WORKFLOW](https://github.com/bayoishola20/buildings-orthogonality/blob/master/assets/workflow.png "Computation workflow")

## REQUIREMENTS: ##

ArcGIS 10.7.1 ArcPy with advanced license has been used in developing this solution. Script runs in ArcMap itself or python in ArcGIS system folder which looks like this: `C:/Python27/ArcGIS10.7/python.exe`. All packages used were those provided by the ArcPy API and so no additional installation is needed.

## RESULTS ##

| Number of buildings |	Time of execution |
|-------------|-----------------|
|4	|7.98099994659 secs|
|50	|22.1856946297 secs|
|100	|36.2016902341 secs|
|-------------|-----------------|


<u>**Test data**</u>: Data (test_building.shp) is a mini extraction of buildings gotten from OSM here: https://www.geofabrik.de/data/shapefiles.html. Data is projected to UTM before use.

## REFERENCES ##

1.	ArcPy ESRI Developers community: https://community.esri.com/community/developers/gis-developers/python
2.	Douglas, David and Peucker, Thomas, "Algorithms for the reduction of the number of points required to represent a digitized line or its caricature," The Canadian Cartographer 10(2), 112–122 (1973).


**PS: This was an MSc course project @TUDresden.**

