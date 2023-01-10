library("rgdal")
library("tidyverse")
library("sf")
library("raster")
library("ggspatial")
library("viridis")
library("igraph")
library("riverconn")
library("elevatr")
library("gridExtra")
library("ggnetwork")
#remotes::install_github("briatte/ggnetwork") ##install from github ! fixed some bugs 6 for plotting!
library("lwgeom")
library("gridExtra")
#library("corrmorant")  ##version problem try to make virt env with R version under which riverconn was writen
library("RANN")
library("ggpubr")
library("cowplot")
library("data.table")
library("rgeos")
#install.packages("remotes")
#remotes::install_github("michaeldorman/nngeo") 
#for parallel
library("foreach")
library("doParallel")
library("sfnetworks")
remotes::install_github("luukvdmeer/sfnetworks")


#check packages list.functions.in.file() 

#INPUT
##################################################################################
#Input formats
#shape_river   #sf LINESTRING
#shape_basin  #sf POLYGON
#shape_dams #sf MULTIPOLYGON

#download CCM2 gdb from https://ccm.jrc.ec.europa.eu/php/index.php?action=view&id=24
##################################################################################
#INPUT (change to your Location of input)

#window <- "2004"
#path <- paste0("C:/Users/Nutzer/Documents/phd/dammedfish/River_data/CCM/CCM21_LAEA_window",window,"/ccm21/") 
#C:\Users\Vicky\Documents\phd\CCM2
path <- paste0("C:/Users/Vicky/Documents/phd/CCM2/") 
dbname <- "WGS84_W2003.gdb"
#dbname <- paste0("LAEA_W",window,".gdb") 

dir.exists(paste0(path,dbname)) #check path
#ogrListLayers(paste0(path,dbname))

seaout <- st_read(paste0(path,dbname), layer = "SEAOUTLETS" )
segments  <- st_read(paste0(path,dbname), layer = "RIVERSEGMENTS" )
node <- st_read(paste0(path,dbname), layer = "RIVERNODES" )
#cat <- st_read(paste0(path,dbname), layer = "CATCHMENTS" )

#Barriers
#amberpath <- "C:/Users/Nutzer/Documents/phd/dammedfish/River_data/AMBER/atlas.csv"
amberpath <-  "C:/Users/Vicky/Documents/phd/Amber/atlas.csv"
#file.exists(amberpath) #check path

amber <- fread(amberpath)
dam_in <- st_as_sf(amber, coords = c("Longitude_WGS84","Latitude_WGS84"))
st_crs(dam_in) <- "+proj=longlat +datum=WGS84"

##################################################################################
#Chose River Basin (Seaoutlet)

##for some segments no River name is available 
#WSO_ID unique identifier 
seaout_df <- seaout %>% as.data.frame()

##in some cases there where no segments/basins to the WSO_ID of Rivers (Rivers with no Seaoutlet?)
##Example "Oude Rijn" in 2003
unique(seaout_df$NAME)
basinname <- "Seudre"
basin_id <- unique(na.omit(seaout_df[seaout_df$NAME == basinname,]$WSO_ID)) #WSO_ID River Basin ID

#alternative with unique Basin identifier WSO_ID
#for all WSO_IDs
#basin_id <- seaout_df$WSO_ID[642]
#seaout_df[seaout_df$WSO_ID == basin_id,]$NAME #name of river basin if available

#Basin
basn <- seaout[seaout$WSO_ID== basin_id, ] #SEAOUTLETS
#Segments
riv <- segments[segments$WSO_ID== basin_id,] #RIVERSEGMENTS
#Nodes
nodes <- node[node$WSO_ID== basin_id, ] #RIVERNODES

#Match crs
#Basin
shape_basin <- st_transform(basn,st_crs(dam_in))
#Segments
shape_river <- st_transform(riv,st_crs(dam_in))
#Confluences
river_joins <- st_transform(nodes,st_crs(dam_in))#for elevation at node compare with segments!
#river_joins <- st_as_sf(river_joins, coords = c("X_LAEA","Y_LAEA"))
#st_crs(river_joins) <- "+proj=longlat +datum=WGS84"

#Select dams in basin and in Buffer area 
##TODO there is a faster way!
dam <- st_intersection(dam_in, st_buffer(shape_basin,1)) #crop dams with basin and a buffer of 1

# ##plot input
plot(st_geometry(shape_basin))
plot(st_geometry(shape_river),add =TRUE)
plot(st_geometry(river_joins),add =TRUE, col="red")
plot(st_geometry(dam),add =TRUE, col="green")

#2.1 Shapefile Preprocessing  ##kann man kürzen
##################################################################################
#rename attributes
shape_river$UP_CELLS = shape_river$CONT_PIXELS# UP_CELLS The number of accumulated cells is a proxy of the upstream catchment area. 
#CONT_PIXELS Area upstream the From Node drained by the river segment in 100x100 grid cells
shape_river$alt = shape_river$ALT_GRADIENT 
st_geometry(shape_river) <- "geometry"

shape_river_simple <- shape_river %>%
  st_as_sf %>%
  st_union()

shape_river_line <- shape_river %>% st_cast(.,"LINESTRING" ) %>% 
  mutate(id = WSO1_ID)  # 1:nrow(.)

#2.3 Dams preprocessing  ##TODO kürzen
##################################################################################
#dams
dams_to_points <- dam %>% mutate(id = GUID) %>% dplyr::select(id) 
 
x <-dams_to_points

#################################################################################

dam_snap_ccm<-  function(dams, rivershape, max_dist = 10) {
  # 
  #
  max_dist=90
   i= 15
   dams <- x[i,]
   rivershape <- shape_river_line
  
  dams <- dams
  lines <- rivershape %>% st_as_sf()
  nf <- lines[st_nearest_feature(dams, lines), ]
  #str(nf)
  plot(st_geometry(nf))
  
  startnf = st_startpoint (nf)
  endnf = st_endpoint (nf)
  
  plot(st_geometry(endnf),add=TRUE, col="green")
  plot(st_geometry(startnf),add=TRUE, col="pink")
  
  plot(st_geometry(dams),add=TRUE, col="blue")

  if ( as.numeric(st_distance(dams, nf)) < max_dist) {
    nf.points <- nf  %>%  st_as_sf %>%   st_cast("POINT")
    nf.lines <- nf  %>%  st_as_sf %>%   st_cast("LINESTRING")

     nf.lines <- nf  %>%  st_as_sf %>%   st_cast("LINESTRING")
     nf.segments <- nf.lines %>% st_segmentize( max_dist) %>%   st_cast("POINT")
     nf.segments_min <- nf.segments[ which.min(st_distance(nf.segments,dams)),]
     plot(st_geometry(nf.segments_min),add=TRUE, col="grey",pch =9)
     st_geometry(dams) <- st_geometry(nf.segments_min)
     
     
     buff <- st_buffer(dams,20)  %>%  st_as_sf

     nf.splits <- st_split(nf.lines, buff) %>%  st_cast()%>% 
       mutate(LENGTH= as.integer( st_length(.))) %>%  subset(. , LENGTH > 50 )
  str(nf.splits)
     plot(st_geometry(nf.splits),add=TRUE, col="yellow")
     plot(st_geometry(st_startpoint(nf.splits[1,])),add=TRUE, col="green")
     plot(st_geometry(st_endpoint(nf.splits[1,])),add=TRUE, col="pink")
     

     # st_geometry(dams) <- st_geometry(nf.segments_min)
    #plot(st_geometry(nf.segments_min),add=TRUE, col="red")
    
    
    
    
    dams <- st_as_sf(dams)
  }
  return(dams)
}

##################################################################################
#Speed up by parrallel processing
#change to your R lib location and change to location of function (snap function outsorced!)
.libPaths() #find lib path
#setup parallel backend to use many processors
cores=detectCores()
#print(cores)
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# Loop
#clusterEvalQ(cl, .libPaths("C:/Users/Nutzer/AppData/Local/R/win-library/4.2" )) #Change  here your R library path
clusterEvalQ(cl, .libPaths("C:/Users/Vicky/AppData/Local/R/win-library/4.2" )) #Change  here your R library path

funpat <- "functions/" #path of function to run parallel 

#dir.exists(funpat)
#file.exists(paste0(funpat,"dam_snap_ccm.R"))

output <- foreach(i=1:nrow(x), .combine = rbind, 
                  .packages= c('dplyr', 'sf')) %dopar% {
  source(paste0(funpat,"dam_snap_ccm.R")) # That is the main point. Source your Function File here.
  temp <- dam_snap_ccm(x[i,],shape_river_line,max_dist= 90) # use your custom function after sourcing 
  temp
}

stopCluster(cl)

dams_snapped <- output
#str(dams_snapped)
##################################################################################
# #no parallel
# snapped <- list()
# for (i in 1:nrow(x)) {
#   snapped[[i]] <- dam_snap_ccm(x[i,],shape_river_line) #optional change here max dist.
#   cat(paste0(i, " "))
# }
# dams_snapped  <- do.call(rbind, snapped)
##################################################################################

plot(st_geometry(shape_river_line))
plot(st_geometry(x),add =TRUE, col= "red")
plot(st_geometry(dams_snapped),add =TRUE, col= "green")


# Retain dams that were snapped
dams_snapped_reduced <-
  dams_snapped[st_contains(shape_river_simple %>% st_sf(), dams_snapped, prepared = FALSE)[[1]],]

str(dams_snapped_reduced)

plot(st_geometry(dams_snapped_reduced),add =TRUE, col= "blue")

#check:
st_distance(dams_snapped_reduced,shape_river_simple ) %>% sum
#should be 0

##dobblesnap
dams_snapped_reduced_joined <- dams_snapped_reduced %>%
  mutate(cluster =
           st_equals(dams_snapped_reduced, dams_snapped_reduced) %>%
           sapply(., FUN = function(x) paste0(x, collapse = "_"))) %>%
  group_by(cluster) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(id_dam = as.character(1:nrow(.)), pass_u = 0.1, pass_d = 0.8) %>%
  as.data.frame %>%  st_as_sf()

##################################################################################
##Slicing with dams 

# Create junction point shapefile
network_links <- rbind(
  dams_snapped_reduced_joined %>%
    mutate(type = "dam", id_barrier = id_dam) %>%
    dplyr::select(type, id_barrier, pass_u, pass_d),
  river_joins %>%
    dplyr::select(ID) %>%
    mutate(pass_u = NA, pass_d = NA, type = "joint") %>%
    rename(id_barrier = ID, geometry = SHAPE)) %>% ##enter here name of geometry column
  mutate(id_links = 1:nrow(.))
str(network_links)

#2.5 Adding additional attributes and rename geometry column if necessary 
##################################################################################
# #Elevation
##alt = ALT_GRADIENT
river_net_simplified <- shape_river_line %>% 
  rename(alt = ALT_GRADIENT, geometry = SHAPE)  %>%
  mutate(EdgeID = as.character(1:nrow(.)))  #as.character

#plot(st_geometry(river_net_simplified))
##TODO more attributes?

#3. igraph creation
##################################################################################
# Create vertices vector

#' Function to create new segments 
#'
#' @param linksvector a sf dataframe "type"       "id_barrier" "pass_u"     "pass_d"     "geometry"   "id_links"  
#' 
#' @param river_net_simplified a polygon shapefile LINESTRING with EdgeID, FROMNODE, TONODe attributes
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

river_net_sliced <- dam_include(network_links, river_net_simplified)
edgedf <-edges_list(network_links,river_net_sliced)
#str(river_net_sliced)

vertices <- river_net_sliced %>%
  st_drop_geometry %>%
  rename(name = EdgeID) %>%
  mutate(across(c(everything()), as.character))  %>% #, -geom
  mutate(length = as.numeric(SHAPE_Length), alt = as.numeric(alt)) %>% dplyr::select(name, everything())%>%
  mutate(Conname = as.character(1:nrow(.)))  #as.character
head(vertices)

##################################################################################
#HIER try better include
dam_include <- function(networklinks,rivernetwork){ 
  #test 
  #links <- network_links
  #segments<-  river_net_simplified
  
  links <- networklinks
  segments <- rivernetwork 
  
  dams  <-links[links$type=="dam",]
  
  #new Segments
  for (i in 1:nrow(dams)){
    #i <-1
    ndam <- dams[i,]  
    seg <- segments[st_nearest_feature(ndam, segments), ] #%>%  st_combine()
    
    #plot(st_geometry(seg))
    #plot(st_geometry(ndam), add = TRUE, col ="red")
    
    newsegs <- lwgeom::st_split(seg,ndam) %>%
      st_collection_extract(.,"LINESTRING") %>% 
      mutate(LENGTH= as.integer( st_length(.))) 
    
    #plot(st_geometry(newsegs[1,] ), add= TRUE, col = "blue")
    #plot(st_geometry(newsegs[2,] ), add= TRUE, col = "green")
    
    ##TODO more attributes
    bindnewsegs <-rbind( newsegs[1,] %>%
                           mutate(TONODE = as.integer(ndam$id_barrier ), EdgeID =as.character(paste0(EdgeID,"_1"))),
                         newsegs[2,] %>%
                           mutate(FROMNODE = as.integer(ndam$id_barrier ), EdgeID =as.character(paste0(EdgeID,"_2"))))
    
    segments <- rbind(bindnewsegs, segments)  %>%  filter(EdgeID != seg$EdgeID )
  }
  return(segments)
}

##################################################################################
river_graph <- graph_from_data_frame(
  d = edgedf,
  directed = TRUE,
  v = vertices)

# check river_graph edges
igraph::get.edge.attribute(river_graph) %>% names
## [1] "type"       "id_links"   "id_barrier" "pass_u"     "pass_d"

# check river_graph vertices
igraph::get.vertex.attribute(river_graph) %>% names
## [1] "name"     "length"   "ARCID"    "UP_CELLS" "alt"
g2 = river_graph
l = c("MAINDRAIN_ID","LENGTH","CUM_LEN","PIXELS_100M","CATCHMENT_AREA","CONT_PIXELS","DRAIN_KM2","BURNED","CONFIDENCE","WINDOW" , "SHAPE_Length"  )

for (i in l){
  g2<-delete_vertex_attr(g2, i)
}
igraph::get.vertex.attribute(g2) %>% names
river_graph<-g2

##3.1 Editing the igraph object
##################################################################################
# update length attribute
#V(river_graph)$length <- V(river_graph)$length / 10000
hist(V(river_graph)$length)##TODO checken welche length

##save graph 
savepath <-"C:/Users/Nutzer/Documents/phd/dammedfish/River_data/CCM_Vicky/"
write_graph(river_graph, paste0(savepath,basinname,"_igraph.txt"))

#Habitat suitability

# Function for organism that prefers high elevation
suit_fun_high <- function(x){dnorm(x, mean = 1500, sd = 500)*410*sqrt(3*pi)}

# Function for organism that prefers low elevation
suit_fun_low <- function(x){exp(- 0.001*x)}

# Calculate HSI for the igraph nodes
V(river_graph)$HSI_low <- suit_fun_low(V(river_graph)$alt)
V(river_graph)$HSI_high <- suit_fun_high(V(river_graph)$alt)

# Calculate weighted usable length for igraph nodes
V(river_graph)$WUL_low <- V(river_graph)$HSI_low * V(river_graph)$length
V(river_graph)$WUL_high <- V(river_graph)$HSI_high * V(river_graph)$length

# Plot the two response functions
1:10:2500 %>%
  data.frame("Elevation" = .,
             "Low" = suit_fun_low(.), 
             "High" = suit_fun_high(.)) %>%
  pivot_longer(c("Low", "High"), 
               names_to = "Type", values_to = "HSI") %>%
  ggplot() + 
  geom_line(aes(x = Elevation, y = HSI, color = Type))+
  theme_bw()


#3.2 Plotting the igraph object
##################################################################################
#here
#names(river_net_sliced)
# Extract reaches centroids
river_net_simplified_centroids <- river_net_sliced %>% #  mutate(.,geometry = x) %>% 
  st_as_sf() %>% #st_transform(st_crs("+proj=longlat +datum=WGS84")) %>% 
  st_centroid() 

# get the centroids coordinates
coordinates <- river_net_simplified_centroids %>% # %>% 
  dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2]) %>%
  dplyr::select(lat, lon) %>%
  st_set_geometry( NULL)

# give postion as x , y to nodes (this is what ggnetwork does)
V(river_graph)$lat <- coordinates$lat
V(river_graph)$lon <- coordinates$lon

# fortify the igraph object
#sessionInfo()

edges_for_plot <- edgedf %>%
  inner_join(as_data_frame(river_graph, what = "vertices") %>% dplyr::select(name, lon, lat), by = c('from' = 'name')) %>%
  rename(x = lon, y = lat) %>%
  inner_join(as_data_frame(river_graph, what = "vertices") %>% dplyr::select(name, lon, lat), by = c('to' = 'name')) %>%
  rename(xend = lon, yend = lat)

##############
#läuft
ggplot() + geom_sf(data =river_net_sliced)+
  geom_curve(aes(x = x, y = y, xend = xend, yend = yend),
             data = edges_for_plot, curvature = 0.1,
             alpha = 0.5, position = "identity")
###########
#Hier bastel

#alle  egenschaften
edges_for_plot <- edgedf %>%
  inner_join(as_data_frame(river_graph, what = "vertices") , by = c('from' = 'name')) %>%
  rename(x = lon, y = lat) %>%
  inner_join(as_data_frame(river_graph, what = "vertices") , by = c('to' = 'name')) %>%
  rename(xend = lon, yend = lat) ##keep the attributes


ggplot(edges_for_plot, aes(x = x, y = y, xend = xend, yend = yend)) +
coord_fixed() +
 # geom_nodes(data = vertices) +
  geom_edges(alpha = 0.5) +
  scale_color_viridis()+
  theme_blank()+
  ggtitle("High-altitude organism") +
  labs(caption = "Network directionality not shown") 

#to test whether it is running 
ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_fixed() +
  geom_nodes(aes(color = HSI_high)) +
  geom_edges(alpha = 0.5) +
  scale_color_viridis()+
  theme_blank()+
  ggtitle("High-altitude organism") +
  labs(caption = "Network directionality not shown")

###########
#läuft net weils irgendwie die coordinaten zersägt
gg0 <- ggnetwork(river_graph, layout = coordinates %>% as.matrix(), scale = FALSE)
grid.arrange(
  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = HSI_high)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("High-altitude organism") +
    labs(caption = "Network directionality not shown"),

  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = WUL_high)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("High-altitude organism") +
    labs(caption = "Network directionality not shown"),

  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = HSI_low)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("Low-altitude organism") +
    labs(caption = "Network directionality not shown"),

  ggplot(gg0, aes(x = x, y = y, xend = xend, yend = yend)) +
    coord_fixed() +
    geom_nodes(aes(color = WUL_low)) +
    geom_edges(alpha = 0.5) +
    scale_color_viridis()+
    theme_blank()+
    ggtitle("Low-altitude organism") +
    labs(caption = "Network directionality not shown"),
  ncol=2, nrow=2)
##################################################################################

#HIER
#FIND OUT WHY NOT CONNECTED
gsize(river_graph)#number of edges
gorder(river_graph) #number of vertices #Segments
Isolated <- which(degree(river_graph)==0)  #degree of v number of its adjacent edges.
lonely_graph <- subgraph(river_graph, Isolated)
gsize(lonely_graph)#number of edges

gorder(lonely_graph) #number of vertices #Segments

##remove lonely graph if lonely

# river_graph <- delete_vertices(river_graph, V(lonely_graph)$name)
# igraph::get.edge.attribute(lonely_graph) %>% names
# igraph::get.vertex.attribute(lonely_graph) %>% names
##################################################################################
#4. Indices calculations

###HIIER für CONEFOR
confldf <-edgedf  ##edgedf ?
head(confldf) 
fro_to<- confldf[,1:2]
fro_to$pascom <- confldf$pass_u * confldf$pass_d 
head(fro_to)

names(vertices)
head(vertices)
#where are the names!!????
verts <- vertices %>%  dplyr::select(Conname, LENGTH) %>%  as.data.frame()
#TODO translate conname to edge ID and change in edgedf!!!!
head(verts)
dim(verts)

savepath <- paste0( path, "forConefor/") 
dir.exists(savepath)

## Before calculating IIC a binary passability has to be defined
###HIER TODO versteh das mal
E(river_graph)$pass_u_bin <- ifelse(is.na(E(river_graph)$pass_u), NA, 0)
E(river_graph)$pass_d_bin <- ifelse(is.na(E(river_graph)$pass_d), NA, 0)
# edgestest <- as.data.frame (get.edgelist(river_graph))
# head(edgestest)


write.table(fro_to, file = paste0(savepath, basinname, "_edges.txt"), sep = "\t",
            row.names = FALSE, col.names = FALSE) #Connections
write.table(verts, file = paste0(savepath, basinname, "_vertices.txt"), sep = "\t",
                         row.names = FALSE, col.names = FALSE)  #nodes
#savepath

#Catchment IIC
 igraph::get.edge.attribute(river_graph) %>% names
 igraph::get.vertex.attribute(river_graph) %>% names

 
catch_iic <-  index_calculation(graph = river_graph,
                                weight = "length",
                                pass_u = "pass_u_bin",
                                pass_d = "pass_d_bin",
                                param = 3,
                                disp_type = "threshold",
                                index_type = "full")

catch_iic
#PROB 0.08 (0.1 *0.8)


##################################################################################
# Initialize list where to store all the index calculation outputs
index <- list()
lab_index <- list()
letter_index <- list()

##Correct graph attributes!
# 1: Symmetric Dendritic Connectivity Index (no biotic effects)
lab_index[[1]] <- "Symmetric DCI"
letter_index[[1]] <- "A"
index[[1]] <- index_calculation(graph = river_graph,
                                weight = "length",
                                B_ij_flag = FALSE,
                                index_type = "reach",
                                index_mode = "from")

# 2: Asymmetric Dendritic Connectivity Index (no biotic effects)
lab_index[[2]] <- "Asymmetric DCI"
letter_index[[2]] <- "B"
index[[2]] <- index_calculation(graph = river_graph,
                                weight = "length",
                                B_ij_flag = FALSE,
                                dir_fragmentation_type = "asymmetric",
                                index_type = "reach",
                                index_mode = "from")

## Before calculating IIC a binary passability has to be defined
E(river_graph)$pass_u_bin <- ifelse(is.na(E(river_graph)$pass_u), NA, 0)
E(river_graph)$pass_d_bin <- ifelse(is.na(E(river_graph)$pass_d), NA, 0)

# 3: Symmetric Integral Index of Connectivity
lab_index[[3]] <- "IIC"
letter_index[[3]] <- "C"
index[[3]] <- index_calculation(graph = river_graph,
                                weight = "length",
                                pass_u = "pass_u_bin",
                                pass_d = "pass_d_bin",
                                param = 3,
                                disp_type = "threshold",
                                index_type = "reach",
                                index_mode = "from")

## Defining uniform weights for IIC
V(river_graph)$unif_w <- 1

# 4: Integral Index of Connectivity with uniform weights
lab_index[[4]] <- "IIC with uniform weights"
letter_index[[4]] <- "D"
index[[4]] <- index_calculation(graph = river_graph,
                                weight = "unif_w",
                                dir_fragmentation_type = "asymmetric",
                                pass_u = "pass_u_bin",
                                pass_d = "pass_d_bin",
                                param = 3,
                                disp_type = "threshold",
                                index_type = "reach",
                                index_mode = "from")

# 5: Population Connectivity Index for lowland fish
lab_index[[5]] <- "PCI lowland fish"
letter_index[[5]] <- "E"
index[[5]] <- index_calculation(graph = river_graph,
                                weight = "WUL_low",   ## weighted usable length
                                param = 0.8,
                                index_type = "reach",
                                index_mode = "from") 

# 6: Population Connectivity Index for upland fish
lab_index[[6]] <- "PCI upland fish"
letter_index[[6]] <- "F"
index[[6]] <- index_calculation(graph = river_graph,
                                weight = "WUL_high", 
                                param = 0.8,
                                index_type = "reach",
                                index_mode = "from") 


# 7: Population Connectivity Index for lowland passive drifter
## note that units of param_d must be consistent with length field (10s of km)
lab_index[[7]] <- "PCI lowland passive"
letter_index[[7]] <- "G"
index[[7]] <- index_calculation(graph = river_graph,
                                weight = "WUL_low", 
                                dir_distance_type  = "asymmetric",
                                disp_type = "threshold", 
                                param_u = 0, 
                                param_d = 3,
                                index_type = "reach",
                                index_mode = "from")

# 8: Population Connectivity Index for upland passive drifter
## note that units of param_d must be consistent with length field (10s of km)
lab_index[[8]] <- "PCI upland passive"
letter_index[[8]] <- "H"
index[[8]] <- index_calculation(graph = river_graph,
                                weight = "WUL_high", 
                                dir_distance_type  = "asymmetric",
                                disp_type = "threshold",
                                param_u = 0, 
                                param_d = 3,
                                index_type = "reach",
                                index_mode = "from")

#####
#plots
# Initialize empty list
plot_list <- list()

# iterate through list length
for (i in 1:length(index)) {
  
  # Join d_i information with dams shapefile
  river_plot <- river_net_simplified %>%
    mutate(name = as.character(EdgeID)) %>%
    left_join(index[[i]]) %>%
    mutate(rank = rank(desc(index)))
  
  # plot
  plot_iter <- ggplot() +
    coord_fixed() +
    ggspatial::layer_spatial(river_plot, aes(color = rank))+
    scale_color_viridis(name = "Ranking", direction = -1)+
    theme_void()+
    guides(size = FALSE) +
    theme(legend.position = "bottom")+
    ggtitle(paste0(letter_index[[i]], ") ", lab_index[[i]]))
  
  legend <- cowplot::get_legend(plot_iter)
  
  plot_list[[i]] <- plot_iter + 
    theme(legend.position = "none")
  
}

plot_list[[length(index)+1]] <- legend

plot_final <- ggpubr::ggarrange( plotlist = plot_list)
plot_final

#save
###
ggsave(plot = plot_final, paste0(savepath, "Figures/", basinname, "_habitat_prioritization.jpeg"),  
       width = 18, height = 12, units = "cm")
##

#correlation plot is not working because package is not working
#4.2 Barriers prioritization
##################################################################################
barriers_metadata <- data.frame("id_barrier" =  E(river_graph)$id_barrier[!is.na(E(river_graph)$id_barrier)],
                                "pass_u_updated" = 1,
                                "pass_d_updated" = 1)
head(barriers_metadata)

# Initialize list where to store all the index calculation outputs
d_index <- list()
lab_d_index <- list()
letter_d_index <- list()

# 1: Symmetric Dendritic Connectivity Index (no biotic effects)
lab_d_index[[1]] <- "Symmetric DCI"
letter_d_index[[1]] <- "A"
d_index[[1]] <- d_index_calculation(graph = river_graph, 
                                    barriers_metadata = barriers_metadata,
                                    B_ij_flag = FALSE, 
                                    parallel = TRUE,
                                    ncores = 7)

# 2: Asymmetric Dendritic Connectivity Index (no biotic effects)
lab_d_index[[2]] <- "Asymmetric DCI"
letter_d_index[[2]] <- "B"
d_index[[2]] <- d_index_calculation(graph = river_graph, 
                                    barriers_metadata = barriers_metadata, 
                                    B_ij_flag = FALSE, 
                                    parallel = TRUE,
                                    ncores = 7,
                                    dir_fragmentation_type = "asymmetric")

## Before calculating IIC a binary passability has to be defined
E(river_graph)$pass_u_bin <- ifelse(is.na(E(river_graph)$pass_u), NA, 0)
E(river_graph)$pass_d_bin <- ifelse(is.na(E(river_graph)$pass_d), NA, 0)

# 3: Symmetric Integral Index of Connectivity
lab_d_index[[3]] <- "Symmetric IIC"
letter_d_index[[3]] <- "C"
d_index[[3]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_low",
                                    pass_u = "pass_u_bin",
                                    pass_d = "pass_d_bin",
                                    param = 3,
                                    disp_type = "threshold",
                                    parallel = TRUE,
                                    ncores = 7)

## Defining uniform weights for IIC
V(river_graph)$unif_w <- 1

# 4: Integral Index of Connectivity with uniform weights
lab_d_index[[4]] <- "IIC with uniform weights"
letter_d_index[[4]] <- "D"
d_index[[4]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "unif_w",
                                    dir_fragmentation_type = "asymmetric",
                                    pass_u = "pass_u_bin",
                                    pass_d = "pass_d_bin",
                                    param = 3,
                                    disp_type = "threshold",
                                    parallel = TRUE,
                                    ncores = 7)

# 5: Population Connectivity Index for lowland fish
lab_d_index[[5]] <- "PCI lowland fish"
letter_d_index[[5]] <- "E"
d_index[[5]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_low", 
                                    param = 0.8,
                                    parallel = TRUE,
                                    ncores = 7) 

# 6: Population Connectivity Index for upland fish
lab_d_index[[6]] <- "PCI upland fish"
letter_d_index[[6]] <- "F"
d_index[[6]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_high", 
                                    param = 0.8,
                                    parallel = TRUE,
                                    ncores = 7)


# 7: Population Connectivity Index for lowland passive drifter
## note that units of param_d must be consistent with length field (10s of km)
lab_d_index[[7]] <- "PCI lowland passive"
letter_d_index[[7]] <- "G"
d_index[[7]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_low", 
                                    dir_distance_type  = "asymmetric",
                                    disp_type = "threshold", 
                                    param_u = 0, 
                                    param_d = 3,
                                    parallel = TRUE,
                                    ncores = 7)

# 8: Population Connectivity Index for upland passive drifter
## note that units of param_d must be consistent with length field (10s of km)
lab_d_index[[8]] <- "PCI upland passive"
letter_d_index[[8]] <- "H"
d_index[[8]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "WUL_high", 
                                    dir_distance_type  = "asymmetric",
                                    disp_type = "threshold",
                                    param_u = 0, 
                                    param_d = 3,
                                    parallel = TRUE,
                                    ncores = 7)

# 9: Catchment Area Fragmentation Index
## note that units of param_d must be consistent with length field (10s of km)
lab_d_index[[9]] <- "CAFI"
letter_d_index[[9]] <- "I"
d_index[[9]] <- d_index_calculation(graph = river_graph,
                                    barriers_metadata = barriers_metadata, 
                                    weight = "UP_CELLS", 
                                    dir_distance_type  = "symmetric",
                                    B_ij_flag = FALSE,
                                    parallel = TRUE,
                                    ncores = 7)

# Initialize empty list
plot_list <- list()

# iterate through list length
for (i in 1:length(d_index)) {
  
  # Join d_i information with dams shapefile
  network_links_plot <- network_links %>%
    filter(type == "dam") %>% 
    left_join(d_index[[i]]) %>%
    mutate(d_rank = rank(desc(d_index)))
  
  # plot
  plot_iter <- ggplot() +
    coord_fixed() +
    ggspatial::layer_spatial(river_net_simplified, color = "gray70")+
    ggspatial::layer_spatial(network_links_plot, 
                             aes(color = d_rank, size = 1/d_rank), alpha = 0.8)+
    scale_color_viridis(name = "Ranking", direction = -1)+
    theme_void()+
    guides(size = FALSE) +
    theme(legend.position = "bottom")+
    ggtitle(paste0(letter_d_index[[i]], ") ", lab_d_index[[i]]))
  
  # legend <- cowplot::get_legend(plot_iter)
  
  plot_list[[i]] <- plot_iter #+ theme(legend.position = "none")
  
}

# plot_list[[length(d_index)+1]] <- legend

plot_final <- ggpubr::ggarrange( plotlist = plot_list, 
                                 common.legend = TRUE, legend = "bottom")
plot_final
ggsave(plot = plot_final, paste0("Figures/" ,basinname,"barriers_prioritization.jpeg"),  
       width = 18, height = 12, units = "cm")

##Speed increase
#Rcpp
##Eigen and Armadillo
#// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>

# // [[Rcpp::export]]
# SEXP armaMatMult(arma::mat A, arma::mat B){
#   arma::mat C = A * B;
#   
#   return Rcpp::wrap(C);
# }
# 
# // [[Rcpp::export]]
# SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
#   Eigen::MatrixXd C = A * B;
#   
#   return Rcpp::wrap(C);
# }
# 
# // [[Rcpp::export]]
# SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
#   Eigen::MatrixXd C = A * B;
#   
#   return Rcpp::wrap(C);
# }


# library(Rfast)
# 
# help(Rfast)
