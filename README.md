# Approach to use riverconn with CMM2 data

This repository contains instructions for prepocessing of CCM2 data for river network connectivity calculations using the R package riverconn.

## Functions 
The folder "functions" contains the main important functions for the preprocessing

**dam_include.R** - Function to create new segments at the location of a barrier and removes the old segment prior barrier
  - only for segment attributes that remain the same when including the dam (categorical attributes?..) or can be calulated (length)
  - for all the other attributes, where the ccm2 segment informations can not be changed the next confluence can be used.
  - Problem: What is with multiple dams per Segment?
  - 

**dam_snap_ccm.R** -  Function to relocate dam positions to overlap with the network
  -in progress (right now next node of Linestring)
  - try to change network aswell and include dam as point of LINESTRING (propper snapping!)
  
 **edgelist.R**  Function to create the edge list for an igraph object 
   - can be made shorter, does what I shall do

## Additional information
 
  **CCM2_Windows2b.JPG** ccm2 data windows
  
  **riverconn_vickyR.R** this is my code where I try new things 
  
  **Formulas_indices.pdf** this i a list of all the indices which shall be included in RivConnect 
  
  
## Whats next?
- Clean up the code
- summarize Rivernodes in one Database (compare with Riversegments!)
- create dashboard/app Rshiny
- 
