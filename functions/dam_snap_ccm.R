#' Function to relocate dam positions to overlap with the network
#'
#' @param dams POINT sf with barrier locations  
#' 
#' @param rivershape with segments  of class sf containing multiple ‘LINESTRING’ geometries
#'
#' @return POINT sf with new geometry overlapping rivershape when in max_dist
#' 
#' @export
#'
#' @examples
#'
dam_snap_ccm<-  function(dams, rivershape, max_dist = 100) {
  dams <- dams
  lines <- rivershape %>% st_as_sf()
  nf <- lines[st_nearest_feature(dams, lines), ]
  ##change here to identifier of segment in dam dataset (probably WSO_ID). test it whether it works!
  #nf <- lines[lines$WSO_ID== ndam$WSO_ID,] 
  
  #plot(st_geometry(nf))
  #plot(st_geometry(dams),add=TRUE, col="blue")
  
  if ( as.numeric(st_distance(dams, nf)) < max_dist) {
    nf.points <- nf  %>%  st_as_sf %>%   st_cast("POINT") 
    #nf.lines <- nf  %>%  st_as_sf %>%   st_cast("LINESTRING") 

    nf.points.min <- nf.points[ which.min(st_distance(dams, nf.points)),]
    st_geometry(dams) <- st_geometry(nf.points.min)
    
    #ignore this:
    #%>%     st_buffer(1)

    #nf_network <- sfnetwork(nf.points, nf.lines, edges_as_lines = TRUE)
    #nf_blend  = st_network_blend(nf_network, dams, tolerance = 90)

    #this also don't work!
    # nf.lines <- nf  %>%  st_as_sf %>%   st_cast("LINESTRING")
    # nf.segments <- nf.lines %>% st_segmentize( max_dist) %>%   st_cast("POINT")
    # nf.segments_min <- nf.segments[ which.min(st_distance(nf.segments,dams)),]
    # st_geometry(dams) <- st_geometry(nf.segments_min)
    #plot(st_geometry(dams),add=TRUE, col="red")
    
    
    dams <- st_as_sf(dams)
  }
  return(dams)
}
