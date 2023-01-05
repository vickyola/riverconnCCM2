#' Function to create new segments at the location of a barrier and removes the old segment prior barrier
#'
#' @param linksvector a sf dataframe "type"       "id_barrier" "pass_u"     "pass_d"     "geometry"   "id_links"  
#' 
#' @param river_net_simplified a polygon shapefile LINESTRING with EdgeID, FROMNODE, TONODE attributes
#' (must be of class sf)
#'
#' @return river_net_simplified with new segments
#' 
#' @export
#'
#' @examples
#'
dam_include <- function(networklinks,rivernetwork){ 
  links <- networklinks
  segments <- rivernetwork 
  dams  <-links[links$type=="dam",]
  
  #new Segments
  for (i in 1:nrow(dams)){
    ndam <- dams[i,]  
    seg <- segments[st_nearest_feature(ndam, segments), ] #%>%  st_combine()
    
    #plot(st_geometry(seg))
    
    newsegs <- lwgeom::st_split(seg,ndam) %>%
      st_collection_extract(.,"LINESTRING") %>% 
      mutate(LENGTH= as.integer( st_length(.))) 
    
    #  plot(st_geometry(newsegs[1,] ), add= TRUE, col = "blue")
    # plot(st_geometry(newsegs[2,] ), add= TRUE, col = "green")
    
    ##TODO more attributes
    bindnewsegs <-rbind( newsegs[1,] %>%
                           mutate(TONODE = as.integer(ndam$id_barrier ), EdgeID =as.character(paste0(EdgeID,"_1"))),
                         newsegs[2,] %>%
                           mutate(FROMNODE = as.integer(ndam$id_barrier ), EdgeID =as.character(paste0(EdgeID,"_2"))))
    
    segments <- rbind(bindnewsegs, segments)  %>%  filter(EdgeID != seg$EdgeID )
  }
  return(segments)
}
