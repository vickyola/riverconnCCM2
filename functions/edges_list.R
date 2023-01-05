#' Function to create the edge list for an igraph 
#'
#' @param networklinks sf data.frame with POINT geometrys containing the confluence and barriers  
#' 
#' @param rivernetwork a polygon shapefile LINESTRING with EdgeID, FROMNODE, TONODE attributes (hast do have the new segments created by dams inclueded!)
#' (must be of class sf)
#'
#' @return river_net_simplified with new segments
#' 
#' @export
#'
#' @examples
#'
edges_list<- function(networklinks,rivernetwork){ 
  links <- networklinks
  segments <- rivernetwork
  #Connections Links
  dirlist <- list()#edges
  for (i in unique(links$id_barrier)){
    dirlist[[i]] <- data.frame(
      "from" =  segments$EdgeID[segments$TONODE == i ][1] ,
      "from2" =  segments$EdgeID[segments$TONODE == i][2] ,
      "from3" =  segments$EdgeID[segments$TONODE == i][3] , ##cases with more than 3 from segments??
      "to" =   segments$EdgeID[segments$FROMNODE ==i][1],
      "type" = links$type[links$id_barrier == i ],
      "id_links" = links$id_links[links$id_barrier == i ], 
      "id_barrier" = i ,
      "pass_u" = links$pass_u[links$id_barrier == i ],
      "pass_d" = links$pass_d[links$id_barrier == i ]
    )
  }
  #clean it 
  dirdf <- do.call(rbind, dirlist) %>%  st_drop_geometry( )
  dirdfclean <- rbind( dirdf %>% dplyr::select(-from2, -from3, from , to),
                       dirdf %>% dplyr::select(-from, -from3 ,from2, to) %>%
                         rename(from = from2),
                       dirdf %>% dplyr::select(-from2, -from , from3, to ) %>% 
                         rename(from = from3))
  dfreturn <- dirdfclean[!(is.na(dirdfclean$to)) & !(is.na(dirdfclean$from)), ] ##include in pipe
  
  return(dfreturn)
}
