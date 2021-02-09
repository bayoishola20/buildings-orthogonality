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

[polygon](https://github.com/bayoishola20/buildings-orthogonality/blob/master/assets/polygons.png "Building geometry")


**PS: This was an MSc course project @TUDresden.**

