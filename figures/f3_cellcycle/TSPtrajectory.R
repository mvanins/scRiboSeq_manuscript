#' R inplementation of wanderlust
#'
#' @param data Input data matrix.
#' @param s Starting point ID.
#' @param l l nearest neighbours.
#' @param k k nearest neighbours, k < l.
#' @param num_waypoints Number of waypoints to guide the trajectory detection.
#' @param flock_waypoints The number of times for flocking the waypoints, default is 2.
#' @param num_graphs Number of repreated graphs.
#' @param waypoints_seed The seed for reproducing the results.
#' @param metric Distance calculation metric for nearest neighbour detection.
#' @param voting_scheme The scheme of voting.
#' @param band_sample Boolean, if band the sample
#' @param partial_order default NULL
#' @param verbose Boolean, if print the details
#'
#' @return a list containing Trajectory, Order, Waypoints
#' @author Hao Chen
#' @importFrom RANN nn2
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples
#' set.seed(15)
#' shuffled_iris <- iris[sample(150, 150, replace = FALSE), ]
#' data <- shuffled_iris[,1:4]
#' data_label <- shuffled_iris[,5]
#' wishbone <- Rwanderlust(data = data, num_waypoints = 100, waypoints_seed = 2)
#' pd1 <- data.frame(id = wishbone$Trajectory, label=data_label, stringsAsFactors = FALSE)
#' pd2 <- data.frame(id = seq_along(row.names(data)), label=data_label, stringsAsFactors = FALSE)
#' #ggplot(pd1, aes(x=id, y=id, colour = label)) + geom_point() + theme_bw()
#' #ggplot(pd2, aes(x=id, y=id, colour = label)) + geom_point() + theme_bw()


#' downsample points to achieve similar spatial densities
#'
#' @param tbl data table
#' @param x name of column with x coordinate
#' @param y name of column with y coordinate
#' @param rounds number of rounds to downsample
#' @param final final number of points to return
#' @param noise_fraction fraction of mean nearest neighbour distance to add to calculated distances
#' @param seed seed (optional)
#' 
#' @return a downsampled tbl
#' @author Michael VanInsberghe
#' @importFrom RANN nn2
#' @import dplyr
#' @import tidyr
#' @examples
#' set.seed(112)
#' data <- rbind(data.frame(x = rnorm(10000, mean = 0, sd = 0.5),
#'                          y = rnorm(10000, mean = 0, sd = 0.5),
#'                          group = "1"),
#'               data.frame(x = rnorm(1000, mean = 0, sd = 1),
#'                          y = rnorm(1000, mean = 0, sd = 1),
#'                          group = "2"))
#' downsampled <- data %>%
#'   downsample_nn()
#' #ggplot(data,aes(x=x,y=y,colour=group))+geom_point()
#' #ggplot(downsampled,aes(x=x,y=y,colour=group))+geom_point()
downsample_nn <- function(tbl, x=x, y=y, rounds = 10, nn_final = 500, final = 50, noise_fraction = 0.25, seed = NULL) {
  if(!is.null(seed)){
    set.seed(seed)
  }
  downsample_fraction = (nn_final/nrow(tbl))^(1/(rounds))
  tblds <- tbl
  for(i in 1:rounds){
    nn_dist <- nn2(tblds %>% select(.data[[x]],.data[[y]]), k=2)$nn.dists[,2] 
    tblds$nn_dist <- nn_dist + runif(nrow(tblds))*(noise_fraction*mean(nn_dist))
    tblds <- tblds %>%
      slice_sample(n=round(downsample_fraction * nrow(tblds)), weight_by = nn_dist)
  }

  return(tblds %>% slice_sample(n=final) %>% select(-nn_dist))
}


#' check orientation of path
#' 
#' @param a dataframe with x,y coordinates in the first, and second columns
#'
#' @return TRUE if clockwise, FALSE if counterclockwise
#' @author Michael VanInsberghe
#' @import dplyr
is_clockwise <- function(a){
  cw <- a %>%
    mutate(direction_product = (.[[1]] - lag(.[[1]]))*(.[[2]] + lag(.[[2]])))

  if(sum(cw$direction_product, na.rm=TRUE) < 0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

#' Perform TSP simulations
#' Adds point distance and path distance columns to input tibble
#' 
#' @param tbl with x and y coordinates
#' @param x name of column with x coordinate
#' @param y name of column with y coordinate
#' @param start_list list of possible start cells
#' @param downsample_rounds number of rounds to downsample
#' @param downsample_points final number of points to return
#' @param downsample_noise_fraction fraction of mean nearest neighbour distance to add to calculated distances
#' 
#' @return tibble with added point_distance and path_distance columns
#' @author Michael VanInsberghe
#' @import dplyr
#' @import TSP
#' @import tidyr
solveTour <- function(tbl, iteration, x=x, y=y, start_list, downsample_rounds = 10, downsample_nn = 500, downsample_final = 50, downsample_noise_fraction = 0.25){
  includes_start <- FALSE
  while(!includes_start){
    ds <- tbl %>%
      downsample_nn(x = x, y = y, rounds = downsample_rounds, nn_final = downsample_nn, final = downsample_final, noise_fraction = downsample_noise_fraction)
    
    start_cells <- pull(ds %>% filter(CB %in% start_list) %>% select(CB))
    if(length(start_cells) < 1){
      next
    }else{
      includes_start <- TRUE
    }
  }

  tsp <- TSP(dist(as.matrix(ds[,c(x, y)])))

  #tsp <- TSP(dist(as.matrix(ds %>% select(.data[[x]], .data[[y]]))))
  tour <- solve_TSP(tsp, method = "concorde", verbose = FALSE)
  
  # find cut cell
  if(is_clockwise(ds[names(tour),c(x, y)])){
    # clockwise
    cut_index <- min(which(names(tour) %in% start_cells))
    if(all(c(first(names(tour)), last(names(tour))) %in% start_cells)){
      # start cells are split between start and end of tour, only consider last half
      cut_name <- names(tail(tour, round(length(tour)/2)))[min(which(names(tail(tour, round(length(tour)/2))) %in% start_cells))]
    }else{
      cut_name <- names(tour)[min(which(names(tour) %in% start_cells))]
    }
    tour_order <- names(cut_tour(tour, names(tour)[cut_index], exclude_cut = FALSE))
  }else{
    # counterclockwise
    cut_index <- max(which(names(tour) %in% start_cells))+1
    if(cut_index > length(tour)){
      # start cells are split between start and end of tour, only condisder first half
      cut_index = max( c(0, which(names(tour[1:round(length(tour)/2)]) %in% start_cells)) ) +1
    }
    tour_order <- rev(names(cut_tour(tour, names(tour)[cut_index], exclude_cut = FALSE)))
  }

  ds$tour_order <- 0
  ds[tour_order, "tour_order"] <- 1:nrow(ds)

  ds <- ds %>%
    arrange(tour_order) %>%
    mutate(point_distance = sqrt((.data[[x]]-lag(.data[[x]]))^2+(.data[[y]]-lag(.data[[y]]))^2)) %>%
    mutate(point_distance = replace_na(point_distance,0)) %>%
    mutate(path_distance = cumsum(point_distance),
           iteration = iteration)

  # return(ds %>% select(-.data[[x]], -.data[[y]]))
  return(ds %>% select(CB, tour_order, point_distance, path_distance, iteration))

}

# from https://stackoverflow.com/questions/37750937/doparallel-package-foreach-does-not-work-for-big-iterations-in-r
idivix <- function(n, chunks) {
  i <- 1
  it <- idiv(n, chunks=chunks)
  nextEl <- function() {
    m <- nextElem(it)  # may throw 'StopIterator'
    value <- list(i=i, m=m)
    i <<- i + m
    value
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}

