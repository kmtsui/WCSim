###############################################################
## L.Anthony 01/12/2021                                      ##
## Readme for WCTE version of WCSim with multiple geometries ##
###############################################################

You can control the version of geometry in WCSim.mac. 
There are 3 new versions of WCTE geometry that can be controlled here:
18c = 18 columns 5 rows. Same as the original design but with dimensions from CAD.
16c = 16 columns 5 rows with dimensions taken from CAD
16cShort = 16 columns 4 rows with dimensions from CAD

The code which defines the geometries is src/WCSimDetectorConfigs.C
Here the size and number of columns and rows are set.
