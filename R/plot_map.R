#' Plot data on Ecuador map
#'
#' @param df A data frame with the following columns:
#'               col_1: the name of the region of Ecuador OR the id
#'               col_2: the data value to be plotted
#' @param geo2 Object created in "SETUP: Import Ecuador GIS data" section
#' @param title Map title
#' @param limits Limits for the color scale
#' @return A ggplot2 map object

plot_map <- function(df, geo2, title, limits=NULL, sfg="purple") {
  
  if (sfg=="purple") {
    sfg2 <- scale_fill_gradient(
      low = "#FFFFFF",
      high = "#280AA1",
      limits = limits
    )
  }
  if (sfg=="rwg") {
    sfg2 <- scale_fill_gradient2(
      low = "#8F2E2E",
      mid = "#FFFFFF",
      high = "#2E8F44",
      midpoint = 0,
      limits = limits
    )
  }
  
  names(df) <- c("id", "Value")
  
  if (df$id[1]==1) {
    df %<>% mutate( id=reg_names[id] )
  }
  
  geo3 = inner_join(geo2, df, by=c("id"))
  
  map_plot <- ggplot(data=geo3) +
    geom_polygon(aes(x=long, y=lat, group=group, fill=Value),
                 color="black") +
    theme_void() +
    coord_map() +
    sfg2 +
    labs(title=title)
  
  return (map_plot)
  
}
