# Approach to use riverconn with CMM2 data

This repository contains instructions for prepocessing of CCM2 data for river network connectivity calculations using the R package riverconn.

## functions 
The folder "functions" contains the main important functions for the preprocessing

**dam_include.R** - Function to create new segments at the location of a barrier and removes the old segment prior barrier
  - only for segment attributes that remain the same when including the dam (categorical attributes?..) or can be calulated (length)

**dam_snap_ccm.R** -  Function to relocate dam positions to overlap with the network
  -in progress
  - try to change network aswell and include dam as point of LINESTRING (propper snapping!)
  
 **edgelist.R**  Function to create the edge list for an igraph 
   - can be made shorter 

