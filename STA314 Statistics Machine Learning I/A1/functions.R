# check existence of the required packages
needed_packages = c("imager", "tidyverse", "tidymodels", "sp", "scales", "cowplot", "dmc")
if(!require(imager)) {
  stop("The imager packages must be installed. Run install.packages(\"imager\") and then try again.")
}
if(!require(tidyverse)) {
  stop("The tidyverse packages must be installed. Run install.packages(\"tidyverse\") and then try again.")
}
if(!require(tidymodels)) {
  stop("The tidymodels packages must be installed. Run install.packages(\"tidymodels\") and then try again.")
}
if(!require(sp)) {
  stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
}
if(!require(scales)) {
  stop("The scales packages must be installed. Run install.packages(\"scales\") and then try again.")
}
if(!require(cowplot)) {
  stop("The cowplot packages must be installed. Run install.packages(\"cowplot\") and then try again.")
}
if(!require(dmc)) {
  stop("The dmc packages must be installed. Run devtools::install_github(\"sharlagelfand/dmc\") and then try again.")
}

library(imager)
library(tidyverse)
library(tidymodels)
library(sp)
library(scales)
library(cowplot)
library(dmc)

#### Given function, change_resolution ####
change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}

#### function1 ####
process_image <- function(image_file_name, k_list)
{
  ## process_image(image_file_name, k_list) takes a image <image_file_name> and number of clusters <k_list> 
  ## to get information requires for cross-stictch
  ## 
  ## Input:
  ## - image_file_name: a string of the file name of a PNG or JPEG image
  ## - k_list: a list of potential number of centres that you want for the clustering, for each k in k_list, 
  ## the program will produce a corresponding kmeans clusterings
  ##
  ## Output:
  ##   – cluster_info: returns a list of tibbles of the information derived from the k_means. This is the only function that
  ## computes a clustering. The list has two tibbles for each k_clusterings at cluster_info[[k]] as the followings:
  ##   * (cluster_info[[k]]$centres) the tidied clusters, their associated RGB values and their nearest DMC thread colour, information to construct colour strips
  ##   * (cluster_info[[k]]$dat) the tidied data after clustered (x, y, R, G, B, cluster, dmc), information for make_patterns function
  ##   * (the last entry, cluster_info$clusterings) the original output of the kclust calls, information for scree plots
  ##
  ## Example:
  ## > cluster_info <- process_image("cat.jpg", c(1:10))
  ## > cluster_info$clusterings
  ## # A tibble: 10 x 6
  ##       k kclust   totss tot.withinss betweenss  iter
  ##   <int> <list>   <dbl>        <dbl>     <dbl> <int>
  ## 1     1 <kmeans> 5493.        5493. -2.71e-10     1
  ## 2     2 <kmeans> 5493.        1932.  3.56e+ 3     1
  ## 3     3 <kmeans> 5493.        1120.  4.37e+ 3     3
  ## 4     4 <kmeans> 5493.         756.  4.74e+ 3     3
  ## 5     5 <kmeans> 5493.         600.  4.89e+ 3     3
  ## 6     6 <kmeans> 5493.         500.  4.99e+ 3     4
  ## 7     7 <kmeans> 5493.         406.  5.09e+ 3     3
  ## 8     8 <kmeans> 5493.         336.  5.16e+ 3     5
  ## 9     9 <kmeans> 5493.         295.  5.20e+ 3     5
  ## 10    10 <kmeans> 5493.         262.  5.23e+ 3     6
  ##
  ## > cluster_info[[5]]$centres
  ## # A tibble: 5 x 9
  ##       R     G      B  size withinss cluster col     hex     dmc                           
  ##    <dbl> <dbl>  <dbl> <int>    <dbl> <fct>   <chr>   <chr>   <chr>                         
  ## 1 0.746 0.549 0.418  13105    125.  1       #BE8C6B #C48E70 Desert Sand (3064)            
  ## 2 0.294 0.106 0.0494  6229    102.  2       #4B1B0D #492A13 Coffee Brown - Very Dark (898)
  ## 3 0.613 0.398 0.284  12642    144.  3       #9C6548 #A06C50 Desert Sand - Very Dark (3772)
  ## 4 0.893 0.708 0.580   8422     98.3 4       #E4B494 #E4BB8E Tan - Light (437)             
  ## 5 0.489 0.261 0.158  10227    130.  5       #7D4328 #7A451F Brown - Medium (433)          
  
  im <- imager::load.image(image_file_name)
  tidy_dat <-as.data.frame(im, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)  # each pixel is stored with their RGB value

  kclusts <-
    tibble(k = k_list) %>%
    mutate(
      kclust = map(k, ~kmeans(x = select(tidy_dat,c(-x,-y)) , centers = .x, nstart=4)),  # set nstart to be something not 1
      glanced = map(kclust, glance), # stores information of the kcluster
    )
  clusterings <- kclusts %>% unnest(cols = c(glanced))

  result <- list() # preparing for storing the big list of data
  for (k in k_list){
    # reconstruct data <dat> to preform kmean again for each iteration
    dat <- as.data.frame(im, wide = "c") %>% rename(R = c.1, G = c.2, B = c.3)
    # construct clusterings for <k> clusters
    kclust <- kmeans(select(dat, -x, -y), centers = k, nstart = 4)
    # save the centres data, its tot.withiness, etc.
    centres <- tidy(kclust)
    # We can also add a column to the tidied centres to add the colour in a way that we can use for plots. The rgb
    # function will do this and display the colour as a hex string.
    centres <- centres %>% mutate(col = rgb(R,G,B))  #creating a new column called "col" of hex string
    # convert the RGB colours to usable colours with DMC library, attach name + dmc and its hex
    dmc_col <- c(1:length(centres$col))
    hex_col <- c(1:length(centres$col))
    for (j in seq(1:length(centres$col))){
      dmc_col[j] <- paste(dmc(centres$col[j], visualize = FALSE)$name, 
                          " (", dmc(centres$col[j], visualize = FALSE)$dmc, ")", sep = "")
      hex_col[j] <- dmc(centres$col[j], visualize = FALSE)$hex
    }
    centres <- centres %>% add_column(hex = hex_col, dmc = dmc_col)
    # Then we need to do the cluster centre replacement. 
    dat <- augment(kclust, dat) %>% rename(cluster = .cluster)
    # Add DMC to each data point based on the cluster they are at
    temp <- centres %>% select(cluster, dmc, hex)
    dat <- dat %>% full_join(temp, by="cluster")
    # put all the data we need for a <k> clustering trial into a list then add it to the 
    # big lst <result>
    result[[k]] = lst(dat, centres)
  }
  result["clusterings"] = lst(clusterings) # information for the scree_plot
  return(result)
}

#### function2 ####
scree_plot <- function(cluster_info)
{
  ## scree_plot takes cluster_info generated earlier for the data and returns a scree plot of the clusterings 
  ## Input:
  ##   – cluster_info: a list of tibbles of the information derived from the k_means. This is the only function that
  ## computes a clustering. The list has two tibbles for each k_clusterings at cluster_info[[k]] as the followings:
  ##   * (cluster_info[[k]]$centres) the tidied clusters, their associated RGB values and their nearest DMC thread colour, information to construct colour strips
  ##   * (cluster_info[[k]]$dat) the tidied data (after clustered), information for make_patterns function
  ##   * (the last entry, cluster_info$clusterings) the original output of the kclust calls, information for scree plots
  ##
  ## Output:
  ##   – scree_plot: plots a ggplot (the scree plot for the clusterings) that has x-axis to be the compoent number (cluster number), y axis to be the 
  ## total within-variances for the number of clustering
  ##
  ## Example:
  ## > cluster_info <- process_image("Liaz.jpg", c(5:8))  # the same picture used for tutorial 2
  ## > scree_plot(cluster_info)
  clusterings <- cluster_info$clusterings
  scree_plot_result <- ggplot(clusterings, aes(k, tot.withinss)) + scale_x_continuous(breaks = pretty_breaks()) + 
    geom_line() + geom_point() + ggtitle("Scree Plot")
  return(scree_plot_result)
}


#### function3 ####
colour_strips <- function(cluster_info)
{
  ## colour_strips(cluster_info) produces colour strips with the DMC colour closest to the cluster centre colour
  ##
  ## Input:
  ##   – cluster_info: a list of tibbles of the information derived from the k_means. This is the only function that
  ## computes a clustering. The list has two tibbles for each k_clusterings at cluster_info[[k]] as the followings:
  ##   * (cluster_info[[k]]$centres) the tidied clusters, their associated RGB values and their nearest DMC thread colour, information to construct colour strips
  ##   * (cluster_info[[k]]$dat) the tidied data (after clustered), information for make_patterns function
  ##   * (the last entry, cluster_info$clusterings) the original output of the kclust calls, information for scree plots
  ##
  ## Output:
  ##   – colour_strips: plots a visualization (diagram) of the chosen clusters colours (dmc colours)
  ##
  ## Example:
  ## > cluster_info <- process_image("Liaz.jpg", c(5:8))
  ## > colour_strips(cluster_info)
  
  # helper function
  square <- function(x, label_size) { 
    scalar = 1
    if( floor(20/label_size) <= 2){
      scalar = 0.5
    }
    ggplot()  + 
      coord_fixed(xlim=c(0,1), ylim = c(0,1)) + theme_void() + 
      theme(plot.background = element_rect(fill = x)) + 
      geom_text(aes(0.5,0.5), label = x , size = scalar*label_size)
  }
  
  for (k in c(1:(length(cluster_info)-1))){
    if (!is.null(cluster_info[[k]])){
      centres <- cluster_info[[k]]$centres
      # visualize the DMCed colours of the cluster centres
      t <- tibble(colours = centres$hex,
                  squares = purrr::map(colours, ~ square(.x, 24/length(colours))))
      title <- ggdraw() +
        draw_label(paste("Colour Strips for k = ", k, " clusters", sep = ""), fontface = 'bold', x = 0, hjust = 0)
      p <- plot_grid(nrow = 1, plotlist = t$squares)
      print(plot_grid(ncol = 1, title, p))
    }
  }
}

#### function4 ####
make_pattern <- function(cluster_info, k, x_size, black_white = FALSE, background_colour = NULL)
{
  ## make_pattern(cluster_info, k, x_size, black_white = FALSE, background_colour = NULL) plots the pattern of cross-stitich
  ## pattern based on the cluster_info
  ##
  ## Input:
  ## – cluster_info: a list of tibbles of information derived from the chosen picture, the output of process_image
  ## – k: a int, the chosen number of clusters after going through the scree plot and colour strips
  ## – x_size: a int, the (approximate) total number of possible stitches in the horizontal direction
  ## – black_white: (logical) print the pattern in black and white (TRUE) or colour (FALSE, default)
  ## – background_colour: a string of dmc number to identify the colour of the background, which should not be stitched in the
  ##                      pattern. (Default is NULL, to not have a colour)
  ##
  ## Output:
  ##   – a cross-stitches pattern ggplot: It produces a cross-stitch pattern that can be followed, complete with a legend that has
  ## thread colour, and a guide grid.
  ##
  ## Note: This function is the only function to use change_resolution(image_df, x_size) from change_resolution.R
  ##
  ## Example:
  ## > cluster_info <- process_image("Liaz.jpg", c(5:8))
  ## > cross_stitches <- make_pattern(cluster_info, k = 5, x_size = 50)

  dat <- change_resolution(cluster_info[[k]]$dat, 50)  # use change_resolution to decrease the quality of the picture
  
  col_frame <- tibble(cluster = cluster_info[[k]]$cluster, col = dat$hex, name = dat$dmc)
  shape_stop <- length(cluster_info[[k]]$centres$cluster) # the number of clusters we will need to plot
  
  # plot the grid, how exciting!
  if (black_white){
    # make plot black_white
    cross_stitch_plot <- dat %>% ggplot(aes(x, y)) + geom_point(aes(col = factor(cluster), shape = factor(cluster)), size = 1.5)  + 
      scale_y_reverse() + scale_shape_manual(values = c(1:shape_stop)) + scale_colour_grey() + ggtitle("Cross-Stitch Pattern")
    return(cross_stitch_plot)
    } else {
      if (is.null(background_colour)){
        # background_colour is NULL
        cross_stitch_plot <- dat %>% ggplot(aes(x, y)) + geom_point(aes(col = factor(cluster), shape = factor(cluster)), size = 1.5)  + 
          scale_colour_manual(values = dat %>% select(cluster, hex) %>% deframe,
                              label =  cluster_info[[k]]$centres %>% select(cluster, dmc) %>% deframe) + 
          scale_y_reverse() + scale_shape_manual(values = c(1:shape_stop)) + ggtitle("Cross-Stitch Pattern")
        return(cross_stitch_plot)
        } else {
          background_colour = paste("(", background_colour, ")", sep = "")
          flag = TRUE
          actual_dmc = NULL
          for (possible_dmc in cluster_info[[k]]$centres$dmc){
            if (grepl(background_colour, possible_dmc, fixed = TRUE)){
              flag = FALSE # valid input
              actual_dmc = possible_dmc
              }
            }
          if (flag){
            # invalid input for the background_colour, tell user to input a correct one again
            print("Wrong input for backgroud_colour! This colour is not one of the colours of the clusters. Please try again.") 
          } else {
            hex_delete = filter(cluster_info[[k]]$centres, dmc == actual_dmc)$hex
            cluster_delete = filter(cluster_info[[k]]$centres, hex == hex_delete)
            cluster_delete_index = cluster_delete$cluster
            new_dat <- dat %>% filter(cluster != cluster_delete_index)
            cross_stitch_plot <- new_dat %>% ggplot(aes(x, y)) + geom_point(aes(col = factor(cluster), shape = factor(cluster)), size = 1.5)  + 
              scale_colour_manual(values = new_dat %>% select(cluster, hex) %>% deframe,
                                  label =  new_dat %>% select(cluster, dmc) %>% deframe) + 
              scale_y_reverse() + scale_shape_manual(values = c(1:shape_stop)) + ggtitle("Cross-Stitch Pattern")
            return(cross_stitch_plot)
            }
        }
    }
}