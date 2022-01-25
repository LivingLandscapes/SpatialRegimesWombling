#=============================================================================
## Map GAM-predicted wombling values

# Create a color ramp for plotting the predicted wombling values
predBreaks <- seq(0, 
                  0.7, 
                  0.005)
predRamp <- rev(colorRampPalette(c("black", 
                                   "darkred",
                                   "red",
                                   "orange", 
                                   "palegoldenrod", 
                                   "greenyellow",
                                   "green",
                                   "green3"))(length(predBreaks) - 1))
predRamp <- c(predRamp, 
              rep(predRamp[length(predRamp)],
                  length(seq(max(predBreaks), max(predFit_dt$fit), 0.005))))

# Recode fit values to be between 0 - 1
new_predFit_dt <- copy(predFit_dt)
new_predFit_dt[fit > 1, "fit"] <- 1
new_predFit_dt[fit < 0, "fit"] <- 0

# Create a map to display the predicted wombling values across Ft Riley
gamPredMap_plot <-
  ggplot() +
  geom_tile(data = new_predFit_dt,
            mapping = aes(x = DMCE, y = DMCN, fill = fit)) +
  scale_fill_gradientn(colors = predRamp,
                       # values = predBreaks,
                       # limits = c(0,
                       #            1),
                       name = "Bird Boundary\nStrength\n(Wombling R^2\nValue)") +
  north(x.min = min(new_predFit_dt$DMCE),
        x.max = max(new_predFit_dt$DMCE),
        y.min = min(new_predFit_dt$DMCN),
        y.max = max(new_predFit_dt$DMCN),
        symbol = 12,
        location = "bottomleft",
        scale = 0.2) +
  facet_wrap(~ Year) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 13),
        plot.margin = unit(c(0,0,0,0.5), "cm")) +
  ylab("Northing") +
  xlab("Easting")
# gamPredMap_plot
# ggsave(plot = gamPredMap_plot,
#        filename = here(FigurePath, "FtRiley_Wombling_HGAM.tiff"),
#        dpi = 300,
#        width = 11, 
#        height = 7)

#=============================================================================
## Map spatial covariance

# Which years to plot?
plot_years <- c(1991, 2005, 2017)

# make spcov rasters into data.frames
spcov_df <- rbindlist(Map(X = c(1, 14, 26),
                          Y = plot_years,
                          function(X, Y) {
                            tempRast <- mask(projRasts_139[[X]],
                                             subset(riley, Maneuver != "IA"))
                            temp <- na.omit(as.data.frame(tempRast, xy = TRUE))
                            names(temp)[ncol(temp)] <- "spcov_139"
                            temp$Year <- Y
                            return(na.omit(temp))
                          }))

# Remove the Impact Area from the map
spcov_df <- inner_join(spcov_df, 
                       as.data.frame(rileyRast, xy = TRUE))[layer_Maneuver != "IA"]

# ## Map spatial covariance
vegBreaks <- seq(-1, 
                 0.02,
                 0.01)
vegRamp <- colorRampPalette(c("darkred", 
                              "red", 
                              "orange", 
                              "palegoldenrod", 
                              "lightskyblue3"))(length(vegBreaks) - 1)
vegRamp <- c(rep(vegRamp[1], 
                 length(seq(-2.93,
                            -1.01, 
                            0.01))),
             vegRamp)

# Map spatial covariance
spcov_maps <-
  ggplot() +
  geom_tile(data = spcov_df,
            aes(x = x,
                y = y,
                fill = spcov_139)) +
  scale_fill_gradientn(colors = vegRamp,
                       # values = rescale(vegBreaks),
                       limits = c(min(spcov_df$spcov_139),
                                  max(spcov_df$spcov_139)),
                       name = "Grass:Woody\nBoundary\nStrength",
                       breaks = c(0, -1, -2),
                       guide = guide_colorbar(reverse = TRUE)) +
  facet_wrap(~ Year, nrow = 1) +
  coord_equal() +
  theme_bw() +
  ggsn::scalebar(x.min = min(spcov_df$x),
                 x.max = max(spcov_df$x),
                 y.min = min(spcov_df$y),
                 y.max = max(spcov_df$y),
                 dist = 5, st.size=2.5, height=0.03,
                 st.dist = .05,
                 dist_unit = "km",
                 transform = FALSE,
                 location="bottomleft" ) +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 13),
        plot.margin = unit(c(0,0,0,0.5), "cm")) +
  ylab("Northing") +
  xlab("Easting")
# spcov_maps
# ggsave(plot = spcov_maps,
#        filename = here(FigurePath,
#                        "weirdTest.tiff"),
#        dpi = 300,
#        height = 7,
#        width = 11)

#=============================================================================
## Map cumulative number of fires

# Map cumulative number of fires
fire_maps <- 
  ggplot() +
  geom_tile(data = fire_dt,
            mapping = aes(x = DMCE, y = DMCN, fill = n_fires)) +
  geom_polygon(fortify(spTransform(riley_boundary, CRS(proj4string(riley)))), 
               mapping = aes(x = long, y = lat),
               fill = NA,
               color = "black") +
  scale_fill_continuous(breaks = seq(1, 9, 2),
                        low = "orange", high = "darkred",
                        name = "Cumulative\nNumber of\nFires") +
  facet_wrap(~ Year, nrow = 1) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 13),
        plot.margin = unit(c(0,0,0,0.5), "cm")) +
  ylab("Northing") +
  xlab("Easting")
# fire_maps
# ggsave(plot = fire_maps,
#        filename = here(FigurePath,
#                        "FireHistory_Maps.tiff"),
#        dpi = 300,
#        height = 7,
#        width = 11)

#=============================================================================
## Make plots for hypothesis test effects plots

# Draw effects plot for number of fires
fireEffects_Plot <-
  mydraw.gam(fit_womble, ci_level = 0.99)[[2]] +
  geom_hline(data = data.frame(y = 0),
             mapping = aes(yintercept = 0), 
             linetype = 3,
             color = "darkred") +
  scale_x_continuous(breaks = seq(0, 8, 2)) +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 10)) +
  ylab(expression("Bird Boundary Strength (wombling)" %->% "")) +
  xlab("Number of Fires")

# Draw effects plot for spatial covariance
spcovEffects_Plot <-
  mydraw.gam(fit_womble, ci_level = 0.99)[[1]] +
  geom_hline(data = data.frame(y = 0),
             mapping = aes(yintercept = 0), 
             linetype = 3,
             color = "darkred") +
  scale_y_continuous(breaks = seq(-0.1, 0.4, 0.1)) +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 10)) +
  ylab(expression("Bird Boundary Strength (wombling)" %->% "")) +
  xlab(expression("" %<-% "Grass:Woody Boundary Strength (spatial covariance)"))

# Combine effects plots
effects_plots <- 
  plot_grid(spcovEffects_Plot, fireEffects_Plot,
            ncol = 2,
            align = "hv",
            labels = c("A", "B"))
# ggsave(plot = effects_plots,
#        filename = here(FigurePath,
#                        "Effects_Plots_jpeg.jpeg"),
#        dpi = 600)

#=============================================================================
# Combine Veg transitions, wombling, and fire history

spcov_womble_fire <- 
  plot_grid(gamPredMap_plot, spcov_maps,  fire_maps, 
            ncol = 1, 
            labels = c("A", "B", "C"),
            label_size = 18,
            align = "hv",
            hjust = -1)
# ggsave(plot = spcov_womble_fire,
#        filename = here(FigurePath,
#                        "SpCov_Womble_Fire_Plot.jpeg"),
#        dpi = 600, 
#        height = 8)

#============================================================
## Plot bird point count locations

# 
pointCount_map <- 
  ggplot() +
  geom_polygon(fortify(spTransform(riley_boundary, CRS(proj4string(riley)))), 
               mapping = aes(x = long, y = lat),
               fill = NA,
               color = "black") +
  geom_polygon(data = fortify(subset(riley, Maneuver == "IA")),
               mapping = aes(x = long, y = lat, group = group),
               fill = "grey50", 
               color = NA) + 
  geom_point(data = outputs[Year == 2017],
             aes(x = DMCE, y = DMCN),
             color = "orange",
             size = 1.5) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 13),
        plot.margin = unit(c(0,0,0,0.5), "cm")) +
  ylab("Northing") +
  xlab("Easting")
# ggsave(plot = pointCount_map,
#        filename = here(FigurePath,
#                        "Supplemental_Map_Plot.tiff"),
#        dpi = 300)

#=============================================================================
## Figure 1: workflow description

flowchart1 <- 
  ggplot() +
  geom_polygon(fortify(spTransform(riley_boundary, CRS(proj4string(riley)))), 
               mapping = aes(x = long, y = lat),
               fill = NA,
               color = "black") +
  geom_point(data = ALL[Year == 2017],
             mapping = aes(x = DMCE, y = DMCN)) +
  coord_equal() +
  ggsn::scalebar(x.min = min(ALL$DMCE),
                 x.max = max(ALL$DMCE),
                 y.min = min(ALL$DMCN),
                 y.max = max(ALL$DMCN),
                 dist = 2.5, st.size=3, height=0.02,
                 st.dist = .05,
                 dist_unit = "km",
                 transform = FALSE,
                 location="bottomleft" ) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.margin = unit(c(0,0,0,0.5), "cm")) +
  ylab("Northing") +
  xlab("Easting")
flowchart2 <- 
  ggplot() +
  geom_polygon(fortify(spTransform(riley_boundary, CRS(proj4string(riley)))), 
               mapping = aes(x = long, y = lat),
               fill = NA,
               color = "black") +
  geom_point(data = ALL[Year == 2017],
             mapping = aes(x = DMCE, y = DMCN, color = GWRr2),
             size = 4) +
  scale_color_gradientn(colors = predRamp,
                        # values = predBreaks,
                        # limits = c(0,
                        #            1),
                        name = "Bird Boundary\nStrength\n(Wombling R^2\nValue)") +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.margin = unit(c(0,0,0,0.5), "cm")) +
  ylab("Northing") +
  xlab("Easting")
flowchart3  <-
  ggplot() +
  geom_tile(data = new_predFit_dt[Year == 2017],
            mapping = aes(x = DMCE, y = DMCN, fill = fit)) +
  scale_fill_gradientn(colors = predRamp,
                       # values = predBreaks,
                       # limits = c(0,
                       #            1),
                       name = "Bird Boundary\nStrength\n(Wombling R^2\nValue)") +
  north(x.min = min(new_predFit_dt$DMCE),
        x.max = max(new_predFit_dt$DMCE),
        y.min = min(new_predFit_dt$DMCN),
        y.max = max(new_predFit_dt$DMCN),
        symbol = 12,
        location = "bottomleft",
        scale = 0.2) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 9),
        plot.margin = unit(c(0,0,0,0.5), "cm")) +
  ylab("Northing") +
  xlab("Easting")

# ggsave(plot = flowchart1,
#        filename = here(FigurePath, "FtRiley_Wombling_Flowchart1.jpeg"),
#        dpi = 600)
# ggsave(plot = flowchart2,
#        filename = here(FigurePath, "FtRiley_Wombling_Flowchart2.jpeg"),
#        dpi = 600)
# ggsave(plot = flowchart3,
#        filename = here(FigurePath, "FtRiley_Wombling_Flowchart3.jpeg"),
#        dpi = 600)