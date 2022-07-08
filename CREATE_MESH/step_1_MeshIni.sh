#!/bin/bash

## make geo file
python Contour2geo.py -r 10000.0 -i shapefile_antarctica_simplify2/contour_antarctica_simplify3.shp -o contour.geo
## mesh
gmsh -2 -format msh2 contour.geo

## convert to Elmer
ElmerGrid 14 2 contour.msh -autoclean

# command to generate vtu from msh file
ElmerGrid 14 5 Contour.msh -autoclean

## convert to vtu
ElmerGrid 2 5 contour
