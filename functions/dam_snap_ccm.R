#' Function to create new segments 
#'
#' @param dams POINT shapephile   
#' 
#' @param rivershape a polygon shapefile 
#' (must be of class sf)
#'
#' @return dams with new geometry overlapping rivershape when in max_dist
#' 
#' @export
#'
#' @examples
#'
dam_snap_ccm<-  function(dams, rivershape, max_dist = 10) {
  dams <- dams
  lines <- rivershape %>% st_as_sf()
  nf <- lines[st_nearest_feature(dams, lines), ]
  if ( as.numeric(st_distance(dams, nf)) < max_dist) {
    nf.points <- nf  %>%  st_as_sf %>%   st_cast("POINT")
    nf.points.min <- nf.points[ which.min(st_distance(nf.points,dams)),]
    st_geometry(dams) <- st_geometry(nf.points.min)
    dams <- st_as_sf(dams)
  }
  return(dams)
}
