Netcdf restart
==============

Convert .result => .nc
----------------------
A tools has been built to assist in this task
*Tools location*: TOOLS/RSLT2NC

Modification to the sif
-----------------------
- *1*: remove call the .result restart in the simulation section of the sif
- *2*: add a solver for restart with UgridReader (make sure restart contain only one time steps)
- *3*: if it failed to read some variable, you need to export them somewhere.

Checks did
----------
- Compare output between .nc and .result restart => differences on the partition edges: Need to be clarified

