# 01_map_figure
# this script generates the map that is figure 1 in Min et al. 2025

# this is updated to use the sf package instead of rgdal (which is no longer available)

library(tidyverse)
library(readxl)
library(here)
# library(rgdal)
library(broom)
library(ggpubr)
library(sf)
library(lubridate)


##### Load shapefiles and prepare data ##### 

# choose a CRS to transform all shapefiles to
WGS_CRS <- 4326


### load USA outline
usa_spdf <- st_read(here("map_files", "USA_adm0.shp"))
usa_proj <- sf::st_transform(usa_spdf, crs = WGS_CRS)

## load USA states and territories
us_states <- st_read(here("map_files", "us_states_territories", "s_05mr24.shp"))
us_states_proj <- sf::st_transform(us_states, crs = WGS_CRS)

## load CAN outline (for Vancouver Island)

BC_spdf <- st_read(here("map_files", "canada", "lpr_000b16a_e.shp"))
BC_proj <- sf::st_transform(BC_spdf, crs = WGS_CRS)

## combine CAN into one so that we don't have the province boundaries
one_CAN_spdf <- st_union(BC_spdf)
# plot(st_union(one_CAN_spdf))
### load Major rivers
rivers_spdf <- st_read(here::here("map_files", "NA_Lakes_and_Rivers","hydrography_l_rivers_v2.shp"))
rivers_proj <- sf::st_transform(rivers_spdf, crs = WGS_CRS)

### load the Columbia river basin boundary
CRB_boundary_spdf <- st_read(here::here("map_files", "Columbia_Basin_Watershed_Boundary","Columbia_Basin_Watershed_Boundary.shp"))
CRB_boundary_proj <- sf::st_transform(CRB_boundary_spdf, crs = WGS_CRS)

# keep only the largest rivers for the inset map
rivers_proj
major_rivers_proj <- filter(rivers_proj, Shape_Leng >= 90092)

# keep only the columbia river basin area
CRB_rivers_proj <- st_crop(rivers_proj, xmin = -125.15, xmax = -113.23, ymin = 44.3, ymax = 48.45)

# subset only the columbia river
# The columbia River estuary and the McNary reservoir aren't considered parts of the Columbia or Snake
# in the shapefile, so we will also need to subset those individually

# let's make those manually, since they have this weird splitting thing
mcnary_reservoir_points <- data.frame(x = c(-119.324995, -119.253927, -119.166036, -119.064413,
                                            -119.014974, -118.946653, -118.937040, -118.985105,
                                            -119.031287, -119.104582, -119.214445, -119.242597),
                                      y = c(45.932081, 45.938767, 45.929216, 45.957864,
                                            45.976955, 46.029897, 46.093741, 46.137532,
                                            46.205528, 46.217883, 46.242583, 46.264898),
                                      NAMEEN = rep("McNary Reservoir", 12))

mcnary_reservoir_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> mcnary_reservoir_linestring

estuary_points <- data.frame(x = c(-122.876, -123.184092, -123.228724, -123.257906, -123.326914,
                                   -123.405878, -123.486216, -123.506472),
                                      y = c(46.0707, 46.183222, 46.162063, 46.146128, 46.149458,
                                            46.187501, 46.250451, 46.258522),
                                      NAMEEN = rep("Columbia River Estuary", 8))

estuary_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> estuary_linestring

# john day mouth
john_day_mouth_points <- data.frame(x = c(-120.513195, -120.530791, -120.547356, -120.553192, -120.564179,
                                          -120.573448, -120.602459, -120.606751, -120.622200,
                                          -120.643658, -120.652584),
                             y = c(45.665883, 45.671761, 45.669362, 45.677398, 45.683395,
                                   45.698024, 45.703059, 45.714447, 45.721878,
                                   45.722357, 45.736976),
                             NAMEEN = rep("John Day River Mouth", 11))

john_day_mouth_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> john_day_mouth_linestring

# walla walla mouth
walla_walla_mouth_points <- data.frame(x = c(-118.929555, -118.865462),
                                    y = c(46.063718, 46.058101),
                                    NAMEEN = rep("Walla Walla River Mouth", 2))

walla_walla_mouth_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> walla_walla_mouth_linestring



# create bbox objects with which to crop
estuary_bbox <- st_bbox(c(xmin=-124, xmax=-123, ymin=46, ymax=46.25), crs = "WGS84")
mcnary_reservoir_bbox <- st_bbox(c(xmin=-119.35, xmax=-118.79, ymin=45.87, ymax=46.27), crs = "WGS84")
rivers_proj %>% 
  st_crop(estuary_bbox) -> estuary_proj
rivers_proj %>% 
  st_crop(mcnary_reservoir_bbox) %>% 
  filter(!(NAMEEN == "Yakima River")) -> mcnary_reservoir_proj

st_geometry(mcnary_reservoir_proj)
rivers_proj %>% 
  filter(NAMEEN %in% c("Columbia River", "Snake River")) -> CS_rivers_proj

# get other major rivers from this DF
subset(rivers_proj, COUNTRY == "USA") -> USA_rivers_proj
unique(USA_rivers_proj$NAMEEN)

rivers_proj %>% 
  filter(NAMEEN %in% c("Fifteenmile Creek", "Deschutes River", 
                       "John Day River", "Umatilla River", "Yakima River", "Wenatchee River", "Entiat River",
                       "Okanogan River", "Methow River",
                       "Walla Walla River", "Tucannon River", "Clearwater River", "Grande Ronde River",
                       "Imnaha River", "Snake River", "Salmon River", "Asotin Creek")) -> major_tribs_proj

rivers_proj %>% 
  filter(NAMEEN %in% c("Clearwater River", "Salmon River",
                         "Deschutes River")) -> major_tribs_proj

unique(major_tribs_proj$NAMEEN)


### load Columbia River Basin streams
CRB_streams_spdf <- st_read(here::here("map_files", "crb_streams_over8m_100k","crb_streams_over8m_100k.shp"))
CRB_streams_proj <- sf::st_transform(CRB_streams_spdf, crs = WGS_CRS)

# # for testing: plot all of the streams
# CRB_streams_map <- ggplot(usa_spdf) +
#   geom_sf() +
#   ylab("Latitude")+
#   xlab("Longitude")+
#   # Lines for CRB streams
#   geom_sf(data = CRB_streams_proj,  color = "gray70", size = 0.2, fill = "gray70") +
#   geom_sf_text(data = CRB_streams_proj, aes(label =  HUCName), size = 1, check_overlap = TRUE) +
#   # Lines for North American Rivers (includes Columbia and Snake)
#   geom_sf(data = CS_rivers_proj, color = "gray70", size = 0.5, fill = "gray70") +
#   geom_sf(data = estuary_linestring, color = "gray70", size = 0.5, fill = "gray70") +
#   geom_sf(data = mcnary_reservoir_linestring, color = "gray70", size = 0.5, fill = "gray70") +
#   coord_sf(ylim = c(44.3,48.45),  xlim = c(-125.15,-113.23))
# 
# # save a test version
# ggsave(here::here("figures", "CRB_streams_map_for_testing.pdf"), plot = CRB_streams_map, height = 12, width = 20)
# 
# # ok, so there are just missing stream segments here (at least segments that don't connect to mainstem Columbia or Snake)
# # we'll have to manually code those

# trim the streams spdf
CRB_stream_names <-  data.frame(unique(CRB_streams_proj$GNIS_NAME))
colnames(CRB_stream_names) <- c("GNIS_NAME")
# Streams we need: Columbia River, Fifteenmile Creek, Deschutes River, 
# John Day River, Umatilla River, Yakima River, Wenatchee River, Entiat River,
# Okanogan River, Methow River,
# Walla Walla River, Tucannon River, Clearwater River, Grande Ronde River,
# Imnaha River, Snake River, Salmon River, Asotin Creek

subset_streams_GNIS_name <- c(CRB_stream_names$GNIS_NAME[grep("Columbia", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Fifteenmile", CRB_stream_names$GNIS_NAME)],
                    # CRB_stream_names$GNIS_NAME[grep("Deschutes", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("John Day", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Umatilla", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Yakima", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Wenatchee", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Entiat", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Walla Walla", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Tucannon", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Clearwater", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Grande Ronde", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Imnaha", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Snake", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Salmon", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Bridge", CRB_stream_names$GNIS_NAME)],
                    # CRB_stream_names$GNIS_NAME[grep("White Creek", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Catherine", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Looking", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Lapwai", CRB_stream_names$GNIS_NAME)],
                    # CRB_stream_names$GNIS_NAME[grep("Klickitat", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Asotin", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Okanogan", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Methow", CRB_stream_names$GNIS_NAME)])
                    # CRB_stream_names$GNIS_NAME[grep("Trout", CRB_stream_names$GNIS_NAME)])

CRB_stream_names <-  data.frame(unique(CRB_streams_proj$HUCName))
colnames(CRB_stream_names) <- c("HUCName")
subset_streams <- c(
  # CRB_stream_names$HUCName[grep("Columbia", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Fifteenmile", CRB_stream_names$HUCName)],
                            # CRB_stream_names$HUCName[grep("Deschutes", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("John Day", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Umatilla", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Yakima", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Wenatchee", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Entiat", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Walla Walla", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Tucannon", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Clearwater", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Grande Ronde", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Imnaha", CRB_stream_names$HUCName)],
                            # CRB_stream_names$HUCName[grep("Snake", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Salmon River", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Bridge", CRB_stream_names$HUCName)],
                            # CRB_stream_names$HUCName[grep("White Creek", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Catherine", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Looking", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Lapwai", CRB_stream_names$HUCName)],
                            # CRB_stream_names$HUCName[grep("Klickitat", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Asotin", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Okanogan", CRB_stream_names$HUCName)],
                            CRB_stream_names$HUCName[grep("Methow", CRB_stream_names$HUCName)])

# Remove the incorrectly selected ones
incorrect_streams <- c("Hoodoo Creek", "Little Wenatchee River", 
                       "Clearwater Creek", "Little Clearwater River", "Snake Creek", 
                       "No Snake Creek", "White Salmon River", "Salmon Falls Creek",
                       "Salmon Creek", "Little White Salmon River", "Little Salmon River",
                       "South Fork Salmon Falls Creek", "Little Salmon Creek",
                       "Big Salmon Creek", "Salmon la Sac Creek", "Trout Lake Creek")

subset_streams <- subset_streams[!(subset_streams %in% incorrect_streams)]

# Drop the ones that were incorrectly selected based on HUCName
# Also drop the north fork clearwater, since it's blocked by dworshak
HUCName_drop_streams <- c(CRB_stream_names$HUCName[grep("Frontal Columbia River", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Willamette River", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Lake Umatilla", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Little Wenatchee River", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Lake Entiat", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("North Fork Clearwater", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Hoodo", CRB_stream_names$HUCName)])

subset_streams <- subset_streams[!(subset_streams %in% HUCName_drop_streams)]


# CRB_streams_proj_subset <- subset(CRB_streams_proj, GNIS_NAME %in% subset_streams)
CRB_streams_proj_subset <- subset(CRB_streams_proj, HUCName %in% subset_streams)

# drop the Salmon River in Oregon
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(GNIS_NAME == "Salmon River" & ecoregion %in% c("Blue Mountains", "Cascades")))

# drop the Clearwater River in the Canadian Rockies
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(HUCName == "Clearwater River" & ecoregion == "Canadian Rockies"))

# drop the weird Clearwater river section that's way off in Montana (?)
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(HUCName == "Clearwater River" & UID > 2116900))

# drop weird fragments
# IDs keep changing so drop based on geography
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(HUCName == "Kachess River-Yakima River" & y_coord > 270000))

# create the missing river segments
# john day, connect to Columbia
# clearwater, connect to snake

# # Walla Walla missing segment
# wawa_stbbox <- st_bbox(c(xmin=-118.463298, xmax=-118.316699, ymin=45.937293, ymax=46.063931), crs = "WGS84")
# 
# CRB_streams_proj %>% 
#   st_crop(wawa_stbbox) -> wawa_missing

#### Make an inset map ####

inset_map <- ggplot(usa_spdf) +
  geom_sf(fill = "gray96") +
  geom_sf(data = one_CAN_spdf, fill = "gray96") +
  # add the CRB layer
  geom_sf(data = CRB_boundary_proj, color = "gray85", linewidth = 0.2, fill = "gray85") +
  # add the state outlines
  geom_sf(data = us_states_proj, fill = "transparent") + 
  # add the major rivers layer
  # geom_sf(data = major_rivers_proj, color = "gray70", linewidth = 0.2, fill = "gray70") +
  coord_sf(ylim = c(36,56),  xlim = c(-132, -90)) +
  annotate(geom = "text", x = -105, y = 38, label = "United States of America",
           fontface = "italic") + 
  annotate(geom = "text", x = -110, y = 52, label = "Canada",
           fontface = "italic") + 
  # remove axis text and labels
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # box for where our larger map is
  annotate(geom = "rect", ymin = 44.3, ymax = 48.45,  
           xmin = -125.15, xmax = -113.23, fill = NA, color = "black")

# ggsave(here::here("figures", "inset_map.png"), plot = inset_map, height = 5, width  = 8)
# ggsave(here::here("figures", "inset_map.pdf"), plot = inset_map, height = 5, width  = 8)


##### MAKE BASE MAP #####
CRB_map <- ggplot(usa_spdf) +
  geom_sf(fill = "gray96") +
  geom_sf(data = one_CAN_spdf, fill = "gray96") +
  # add the CRB layer
  geom_sf(data = CRB_boundary_proj, color = "gray90", linewidth = 0.2, fill = "gray90") +
  # add the state outlines
  geom_sf(data = us_states_proj, fill = "transparent") + 
  ylab("Latitude")+
  xlab("Longitude")+
  # Lines for CRB streams
  geom_sf(data = CRB_streams_proj_subset,  color = "gray20", linewidth = 0.5, fill = "gray20") +
  # Lines for North American Rivers (includes Columbia and Snake, and major tribs)
  geom_sf(data = CS_rivers_proj, color = "gray20", linewidth = 2, fill = "gray20") +
  geom_sf(data = estuary_linestring, color = "gray20", linewidth = 2, fill = "gray20") +
  geom_sf(data = mcnary_reservoir_linestring, color = "gray20", linewidth = 2, fill = "gray20") +
  geom_sf(data = major_tribs_proj, color = "gray20", linewidth = 0.5, fill = "gray20") +
  geom_sf(data = john_day_mouth_linestring, color = "gray20", linewidth = 0.5, fill = "gray20") +
  geom_sf(data = walla_walla_mouth_linestring, color = "gray20", linewidth = 0.5, fill = "gray20") +
  coord_sf(ylim = c(44.3,48.45),  xlim = c(-125.15,-113.23)) +
  # ADD LABELS FOR RIVERS
  # annotate("segment", x = -122.84, y = 47.68, xend = -123.18, yend = 47.68, size = 0.5, lty = 1) + # Columbia River
  annotate("text", x = -122.83, y = 46.15, label = "Columbia\nRiver", size = 4, fontface = 'italic', hjust = 0) + # Columbia River
  annotate("text", x = -121.2, y = 45.44, label = "Fifteenmile Cr.", size = 3, fontface = 'italic', hjust = 1, angle = 45) + # Fifteenmile Cr.
  annotate("text", x = -121.3, y = 44.8, label = "Deschutes R.", size = 3, fontface = 'italic', hjust = 1) + # Deschutes R.
  annotate("text", x = -119.63, y = 44.92, label = "John Day R.", size = 3, fontface = 'italic', hjust = 1) + # John Day R.
  annotate("text", x = -118.98, y = 45.60, label = "Umatilla R.", size = 3, fontface = 'italic', hjust = 1) + # Umatilla R.
  annotate("text", x = -117.95, y = 46.15, label = "Walla Walla R.", size = 3, fontface = 'italic', hjust = 1) + # Walla Walla R.
  annotate("text", x = -117.8, y = 46.366, label = "Tucannon R.", size = 3, fontface = 'italic', hjust = 1) + # Tucannon R.
  annotate("text", x = -118, y = 45.50, label = "Grande Ronde R.", size = 3, fontface = 'italic', hjust = 1) + # Grande Ronde R.
  annotate("text", x = -115.34, y = 46.475, label = "Clearwater R.", size = 3, fontface = 'italic', hjust = 1) + # Clearwater R.
  annotate("text", x = -116.92, y = 45.52, label = "Imnaha R.", size = 3, fontface = 'italic', hjust = 1) + # Imnaha R.
  annotate("text", x = -117.35, y = 44.6, label = "Snake\nRiver", size = 4, fontface = 'italic', hjust = 1) + # Snake R.
  annotate("text", x = -115.6, y = 45.62, label = "Salmon R.", size = 3, fontface = 'italic', hjust = 1) + # Salmon R.
  annotate("text", x = -120.44, y = 46.40, label = "Yakima R.", size = 3, fontface = 'italic', hjust = 1) + # Yakima R.
  annotate("text", x = -120.5, y = 47.42, label = "Wenatchee R.", size = 3, fontface = 'italic', hjust = 1) + # Wenatchee R.
  annotate("text", x = -120.6, y = 47.96, label = "Entiat R.", size = 3, fontface = 'italic', hjust = 1) + # Entiat R.
  annotate("text", x = -120.263715, y = 48.341111, label = "Methow R.", size = 3, fontface = 'italic', hjust = 1) + # Methow R.
  annotate("text", x = -119.434335, y = 48.468745, label = "Okanogan R.", size = 3, fontface = 'italic', hjust = 0) + # Okanogan R.
  annotate("text", x = -117.038870, y = 46.159808, label = "Asotin Cr.", size = 3, fontface = 'italic', hjust = 1) + # Asotin Creek
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  guides(fill = guide_legend(title = "Legend")) +
  # add a north arrow
  annotate("text", x = -124.9, y = 45.14, label = "N", fontface = "bold", size = 8) +
  annotate(geom = "polygon", x = c(-125.1, -124.9, -124.9), y = c(44.6, 45, 44.7), 
           color = "black", fill = "gray50") +
  annotate(geom = "polygon", x = c(-124.9, -124.9, -124.7), y = c(44.7, 45, 44.6), 
           color = "black", fill = "white")

# save this as the base map 
# ggsave(here::here("figures", "CRB_base_map.pdf"), plot = CRB_map, height = 7.5, width  = 12.5)

#### Create the base map, add the inset ####

CRB_map +
  annotation_custom(grob = ggplotGrob(inset_map),
                    xmin = -116.4, xmax = -113,
                    ymin = 46.7, ymax = 48.6) -> CRB_map_plus_inset

coord_sf(ylim = c(36,56),  xlim = c(-132, -90))


#### Add the dams on top of the base map ####
complete_event_det_counts <- read.csv(here::here("intermediate_outputs", "complete_event_det_counts.csv"))

# Make a change - add Lower Monumental and Little Goose to list of dams, highlight in grey.
complete_event_det_counts %>% 
  mutate(dam = ifelse(event_site_name %in% c("GOA - Little Goose Fish Ladder", "LMA - Lower Monumental Adult Ladders"), "dam", dam)) %>% 
  mutate(dam_abbr = ifelse(event_site_name == "GOA - Little Goose Fish Ladder", "LGO",
                           ifelse(event_site_name == "LMA - Lower Monumental Adult Ladders", "LMO", dam_abbr))) -> complete_event_det_counts


# JDA correction
complete_event_det_counts %>% 
  subset(!duplicated(event_site_name)) -> complete_event_det_counts

# Add the first dams that block fish passage: Chief Joseph on the Columbia and Hells Canyon on the Snake
blocking_dams <- data.frame(dam_abbr = c("CJO", "HEC"),
                            lat = c(48.00, 45.17),
                            lon = c(-119.6, -116.68),
                            dam_angle = c(0, 60),
                            label_lat = c(47.8224, 44.97),
                            label_lon = c(-119.424, -116.55))


dam_data_for_plot <- data.frame(dam_abbr = c("BON", "TDA", "JDA", "MCN",
                                               "PRA", "RIS", "RRE", "WEL",
                                               "ICH", "LMO", "LGO", "LGR"),
                                  dam_angle = c(0,0,0,0,
                                                -30,-30,80,-85,
                                                0,0,0,-30))

dam_data_for_plot %>% 
  left_join(., complete_event_det_counts, by = "dam_abbr") -> dam_data_for_plot

# add the label positions: By default place them 0.2 degrees above the dam, but manually adjust
dam_data_for_plot %>% 
  mutate(label_lat = event_site_latitude + 0.2) %>% 
  mutate(label_lon = event_site_longitude) -> dam_data_for_plot

dam_data_for_plot[dam_data_for_plot$dam_abbr == "RIS", "label_lat"] <- 47.362
dam_data_for_plot[dam_data_for_plot$dam_abbr == "RIS", "label_lon"] <- -119.773
dam_data_for_plot[dam_data_for_plot$dam_abbr == "PRA", "label_lon"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "PRA", "label_lon"]-0.05
dam_data_for_plot[dam_data_for_plot$dam_abbr == "PRA", "label_lat"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "PRA", "label_lat"]-0.38
dam_data_for_plot[dam_data_for_plot$dam_abbr == "RRE", "label_lat"] <- 47.541
dam_data_for_plot[dam_data_for_plot$dam_abbr == "RRE", "label_lon"] <- -119.9
dam_data_for_plot[dam_data_for_plot$dam_abbr == "WEL", "label_lat"] <- 47.964
dam_data_for_plot[dam_data_for_plot$dam_abbr == "WEL", "label_lon"] <- -120.23
dam_data_for_plot[dam_data_for_plot$dam_abbr == "LGR", "label_lon"] <- -117.3632
dam_data_for_plot[dam_data_for_plot$dam_abbr == "LMO", "label_lon"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "LMO", "label_lon"]-0.05
dam_data_for_plot[dam_data_for_plot$dam_abbr == "LGO", "label_lon"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "LGO", "label_lon"] + 0.05
dam_data_for_plot[dam_data_for_plot$dam_abbr == "TDA", "label_lon"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "TDA", "label_lon"]-0.05
dam_data_for_plot[dam_data_for_plot$dam_abbr == "JDA", "label_lon"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "JDA", "label_lon"] + 0.05
dam_data_for_plot[dam_data_for_plot$dam_abbr == "MCN", "label_lon"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "MCN", "label_lon"]-0.05
dam_data_for_plot[dam_data_for_plot$dam_abbr == "ICH", "label_lon"] <- dam_data_for_plot[dam_data_for_plot$dam_abbr == "ICH", "label_lon"] - 0.05

# Add Wanapum Dam (whoops!)
wanapum_dam <- data.frame(dam_abbr = c("WAN"),
                          event_site_latitude = c(46.85),
                          event_site_longitude = c(-119.98),
                            dam_angle = c(90),
                            label_lat = c(46.85),
                          label_lon = c(-119.5))

dam_data_for_plot %>% 
  bind_rows(., wanapum_dam) -> dam_data_for_plot



CRB_map_dams <- CRB_map_plus_inset +
  # Add dams - but grey out TDA, JDA, LMO, and LGO since they're not currently included as state delineators (and therefore we aren't estimating overshoot/fallback at them)
  geom_text(data = subset(dam_data_for_plot, !(dam_abbr %in% c("TDA", "JDA", "LGO", "LMO", "WAN"))), 
             aes(x = event_site_longitude, y = event_site_latitude, angle = dam_angle), 
             label = "|", fontface = "bold", size = 9, inherit.aes = FALSE) +
  geom_text(data = subset(dam_data_for_plot, !(dam_abbr %in% c("TDA", "JDA", "LGO", "LMO", "WAN"))), 
            aes(x = label_lon, y = label_lat, label = dam_abbr),
            size = 4.5, inherit.aes = FALSE) +
  # grey out TDA, JDA, LMO, and LGO
  geom_text(data = subset(dam_data_for_plot, dam_abbr %in% c("TDA", "JDA", "LGO", "LMO", "WAN")), 
             aes(x = event_site_longitude, y = event_site_latitude, angle = dam_angle), 
             label = "|", fontface = "bold", size = 9, color = "gray70", inherit.aes = FALSE) +
  # version of the plot without text labels for those dams
  # geom_text(data = subset(dam_data_for_plot, dam_abbr %in% c("TDA", "JDA", "LGO", "LMO", "WAN")), 
  #           aes(x = label_lon, y = label_lat, label = dam_abbr),
  #           size = 4.5, color = "gray70", inherit.aes = FALSE) +
  # add Chief Joseph and Hells Canyon as blocking dams
  geom_text(data = blocking_dams,
            label = "|", fontface = "bold",
             aes(x = lon, y = lat, angle = dam_angle), 
            hjust = 1,
             size = 9, color = "#800026", inherit.aes = FALSE) +
  geom_text(data = blocking_dams, 
            aes(x = label_lon, y = label_lat, label = dam_abbr),
            size = 4.5, color = "#800026", inherit.aes = FALSE)
  
  
  # geom_text(data = blocking_dams,
  #           aes(x = lon, y = lat-0.2, label = dam_abbr), 
  #           size = 6, color = "#800026", inherit.aes = FALSE)



# ggsave(here("figures", "CRB_map_dams.pdf"), CRB_map_dams, height = 7.5, width  = 12.5)

#### Final figure: Add the modeling states on top of the base map ####

# create a dataframe that has locations for all of the mainstem states (will be
# annotated as points)
mainstem_states <- data.frame(state = 1:9,
                              x = c(-122.75,
                                    -120.1,
                                    -119.3,
                                    -120,
                                    -120.4,
                                    -119.4,
                                    -119.7878,
                                    -118.3,
                                    -116.9),
                              y = c(45.65,
                                    45.75,
                                    46.6,
                                    47.1,
                                    47.2,
                                    47.55,
                                    48.5,
                                    46.55,
                                    46))

CRB_map_dams +
  # add segments to connect circles for states that aren't on mainstem
  annotate(geom = "segment", x = -120.24, y = 47.377, xend = -120.4, yend = 47.2) + # 5
  annotate(geom = "segment", x = -119.994, y = 47.7808, xend = -119.4, yend = 47.6) + #6
  annotate(geom = "segment", x = -119.7878, y = 48.0873, xend = -119.7878, yend = 48.5) + #7
  geom_point(data = mainstem_states, aes(x = x, y = y), 
             size = 12, shape = 21, fill = "white", color = "black") +
  geom_text(data = mainstem_states, aes(x = x, y = y, label = state),
            size = 5) -> paper_map_figure

ggsave(here("figures", "paper_figures", "Fig1.1_map.pdf"), paper_map_figure, height = 7.5, width  = 12.5)
ggsave(here("figures", "paper_figures", "Fig1.1_map.png"), paper_map_figure, height = 7.5, width  = 12.5)

# another map figure for presentation - no circles for states,
# and all dams the same color

CRB_map_plus_inset +
  # Add dams - but grey out TDA, JDA, LMO, and LGO since they're not currently included as state delineators (and therefore we aren't estimating overshoot/fallback at them)
  geom_text(data = dam_data_for_plot, 
            aes(x = event_site_longitude, y = event_site_latitude, angle = dam_angle), 
            label = "|", fontface = "bold", size = 9, inherit.aes = FALSE) +
  geom_text(data = dam_data_for_plot, 
            aes(x = label_lon, y = label_lat, label = dam_abbr),
            size = 4.5, inherit.aes = FALSE) +
  geom_text(data = blocking_dams,
            label = "|", fontface = "bold",
            aes(x = lon, y = lat, angle = dam_angle), 
            hjust = 1,
            size = 9, color = "#800026", inherit.aes = FALSE) +
  geom_text(data = blocking_dams, 
            aes(x = label_lon, y = label_lat, label = dam_abbr),
            size = 4.5, color = "#800026", inherit.aes = FALSE) -> presentation_map_figure

ggsave(here("figures", "pres_map.pdf"), presentation_map_figure, height = 7.5, width  = 12.5)


