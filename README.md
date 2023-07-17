# Approach to use riverconn with CMM2 data

This repository contains instructions for preprocessing CCM2 data for river network connectivity calculations using the R package riverconn. You can find the riverconn package on [GitHub](https://github.com/cran/riverconn).

## Functions
The `functions` folder contains the main important functions for the preprocessing of CCM2 data.

**dam_include.R**
- Function to create new segments at the location of a barrier and remove the old segment prior to the barrier.
- This function is used for segment attributes that remain the same when including the dam, or for attributes that can be calculated (e.g., length).
- For other attributes where the CCM2 segment information cannot be changed, the next confluence can be used.
- Note: This function may need further development to handle multiple dams per segment.

**dam_snap_ccm.R**
- Function to relocate dam positions to overlap with the network.
- Currently a work in progress, where dams are snapped to the next node of a Linestring.
- The aim is to change the network as well and include the dam as a point of the LINESTRING (proper snapping).

**edgelist.R**
- Function to create the edge list for an igraph object.
- This function can be shortened, as it performs the required tasks.

## Additional Information
- **CCM2_Windows2b.JPG**: CCM2 data windows screenshot.
- **riverconn_vickyR.R**: This file contains my code where I experiment with new ideas.
- **Formulas_indices.pdf**: A list of all the indices that should be included in RivConnect.

## What's Next?
- Clean up the code.
- Summarize Rivernodes in one database and compare with Riversegments.
- Create a dashboard or R Shiny app.

