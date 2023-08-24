library(XML)

extract_data_from_SVG <- function(filename, x_tick_values, y_tick_values,
                                  x_axis_layer_name = "x-axis",
                                  y_axis_layer_name = "y-axis",
                                  point_layer_name = "points",
                                  loess_span = 0.1)
{
  xml <- xmlTreeParse(filename)
  xml_variables <- xmlRoot(xml)["g", all=T]
  
  x_ticks <- NULL
  y_ticks <- NULL
  points_xy <- data.frame(x=NULL, y=NULL)
  
  for (i in 1:length(xml_variables)) {
    layer_id <- xmlAttrs(xml_variables[[i]])["id"]
    if (layer_id == x_axis_layer_name) {
      # X-axis tick marks
      ticks <- xmlChildren(xml_variables[[i]])
      x_ticks <- rep(NA, length(ticks))
      for (j in 1:length(ticks)) {
        x_ticks[j] <- as.numeric(xmlAttrs(ticks[[j]])["x1"])
      }
    }
    if (layer_id == y_axis_layer_name) {
      # Y-axis tick marks
      ticks <- xmlChildren(xml_variables[[i]])
      y_ticks <- rep(NA, length(ticks))
      for (j in 1:length(ticks)) {
        y_ticks[j] <- as.numeric(xmlAttrs(ticks[[j]])["y1"])
      }
    }
  }
  
  for (i in 1:length(xml_variables)) {
    layer_id <- xmlAttrs(xml_variables[[i]])["id"]
    if (layer_id == point_layer_name) {
      points <- xmlChildren(xml_variables[[i]])
      if (names(points)[1] == "path") {
        # Its one or more paths, get all the path coordinates
        merged_paths_x <- NULL
        merged_paths_y <- NULL
        for (j in 1:length(points)) {
          path_data <- as.character(xmlAttrs(points[[j]])["d"])
          all_path_points_x <- NULL
          all_path_points_y <- NULL
          
          command_locations <- gregexpr("[A-Za-z]", path_data)[[1]]
          for (k in 1:(length(command_locations)-1)) {
            pos <- command_locations[k]
            command <- substr(path_data, pos, pos)
            coords <- substr(path_data, pos+1, command_locations[k+1]-1)
            
            # Add commas
            coords <- paste(substr(coords, 1, 1), gsub("-", ",-", substr(coords, 2, nchar(coords))), sep="")
            xy <- strsplit(coords, ",")[[1]]
            if (command == "c" || command == "C" || command == "S" || command == "s") {
              xy <- tail(xy, 2)
            }
            if (command == "H") {
              xy <- c(xy, tail(all_path_points_y, 1))
            }
            if (command == "h") {
              xy <- c(xy, 0)
            }
            if (command == "V") {
              xy <- c(tail(all_path_points_x, 1), xy)
            }
            if (command == "v") {
              xy <- c(0, xy)
            }
            stopifnot(length(xy) == 2)
            
            if (length(grep("[[:upper:]]", command)) > 0) {
              all_path_points_x <- c(all_path_points_x, as.numeric(xy[1]))
              all_path_points_y <- c(all_path_points_y, as.numeric(xy[2]))
            } else {
              all_path_points_x <- c(all_path_points_x, tail(all_path_points_x, 1) + as.numeric(xy[1]))
              all_path_points_y <- c(all_path_points_y, tail(all_path_points_y, 1) + as.numeric(xy[2]))
            }
          }
          merged_paths_x <- c(merged_paths_x, all_path_points_x)
          merged_paths_y <- c(merged_paths_y, all_path_points_y)
        }
        
        # Get evenly distributed points
        predict_x <- seq(min(x_ticks), max(x_ticks), len=100)
        #use_ix <- which(predict_x > min(merged_paths_x) & predict_x < max(merged_paths_x))
        predict_y <- predict(loess(merged_paths_y ~ merged_paths_x, span=loess_span), predict_x)
        
        #plot(merged_paths_x, merged_paths_y)
        #lines(predict_x, predict_y, col='red', lwd=3)
        
        points_xy <- data.frame(x = predict_x, y = predict_y)
      } else if (names(points)[1] == "polyline") {
        polyline_data <- as.character(xmlAttrs(points$polyline)["points"])
        coords <- strsplit(polyline_data, " ")[[1]]
        points_xy <- data.frame(x=rep(NA, length(coords)), y=rep(NA, length(coords)))
        for (j in 1:length(coords)) {
          if (coords[j] != "") {
            xy <- strsplit(coords[j], ",")[[1]]
            points_xy$x[j] <- as.numeric(xy[1])
            points_xy$y[j] <- as.numeric(xy[2])
          }
        }
        points_xy <- points_xy[!is.na(points_xy$x),]
      } else {
        # They should be primitives - just get the coordinates
        points_xy <- data.frame(x=rep(NA, length(points)), y=rep(NA, length(points)))
        for (j in 1:length(points)) {
          points_xy$x[j] <- as.numeric(xmlAttrs(points[[j]])["cx"])
          points_xy$y[j] <- as.numeric(xmlAttrs(points[[j]])["cy"])
        }
      }
    }
  }
  
  stopifnot(length(x_tick_values) == length(x_ticks))
  stopifnot(length(y_tick_values) == length(y_ticks))
  
  x_regression <- lm(x_tick_values ~ sort(as.numeric(x_ticks)))
  y_regression <- lm(y_tick_values ~ sort(as.numeric(y_ticks)))
  
  ordering <- order(points_xy$x)
  
  data <- data.frame(row.names = 1:nrow(points_xy))
  data$x <- points_xy$x[ordering] * coef(x_regression)[2] + coef(x_regression)[1]
  data$y <- points_xy$y[ordering] * coef(y_regression)[2] + coef(y_regression)[1]
  return(data)
}


extract_barplot_data_from_SVG <- function(filename, y_tick_values,
                                          y_axis_layer_name = "y-axis",
                                          point_layer_name = "points",
                                          loess_span = 0.1)
{
  xml <- xmlTreeParse(filename)
  xml_variables <- xmlRoot(xml)["g", all=T]
  
  y_ticks <- NULL
  points_xy <- data.frame(x=NULL, y=NULL)
  
  for (i in 1:length(xml_variables)) {
    layer_id <- xmlAttrs(xml_variables[[i]])["id"]
    if (layer_id == y_axis_layer_name) {
      # Y-axis tick marks
      ticks <- xmlChildren(xml_variables[[i]])
      y_ticks <- rep(NA, length(ticks))
      for (j in 1:length(ticks)) {
        y_ticks[j] <- as.numeric(xmlAttrs(ticks[[j]])["y1"])
      }
    }
  }
  
  for (i in 1:length(xml_variables)) {
    layer_id <- xmlAttrs(xml_variables[[i]])["id"]
    if (layer_id == point_layer_name) {
      points <- xmlChildren(xml_variables[[i]])
      
      if (names(points)[1] == "line") {
        # They should be lines - just get the coordinates
        points_xy <- data.frame(x=rep(NA, length(points)), y=rep(NA, length(points)))
        for (j in 1:length(points)) {
          points_xy$x[j] <- mean(c(as.numeric(xmlAttrs(points[[j]])["x1"]), as.numeric(xmlAttrs(points[[j]])["x2"])))
          points_xy$y[j] <- as.numeric(xmlAttrs(points[[j]])["y1"])
        }
      } else {
        stop("Unsupported primitive for barplots - tops of the bars should be marked with a line")
      }
    }
  }
  
  stopifnot(length(y_tick_values) == length(y_ticks))
  y_regression <- lm(y_tick_values ~ sort(as.numeric(y_ticks)))
  
  ordering <- order(points_xy$x)
  
  data <- data.frame(row.names = 1:nrow(points_xy))
  data$y <- points_xy$y[ordering] * coef(y_regression)[2] + coef(y_regression)[1]
  return(data)
}