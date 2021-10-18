ausplots <- function(overall_minimum_abundance = 100, # Minimum overall number of individuals for the species not to be considered 'rare'
                     plot_minimum_abundance = 5, # Minimum number of individuals in a given plot for the species not to be considered 'rare'
                     thinning_probability = 0, # Remove this percent of individuals
                     condition = "Alive", # Select only live or dead individuals
                     use_marks = FALSE, # Use individuals' diameters as marks?
                     to_consider = 48, # How many of the plots to consider?
                     jitter = 1e-5) { # Jitter locations with this std error
  # Constants related to the dataset
  square_width <- 100
  long_square_width <- 160
  
  # to_consider should not be larger than the total number of plots
  to_consider <- max(1, min(to_consider, 48))
  
  # Load the required data
  data <- utils::read.csv(file = "https://shared.tern.org.au/attachment/0e503109-2fb6-4182-969f-2d570abdbabd/data_large_tree_survey.csv", 
                          stringsAsFactors = FALSE, 
                          allowEscapes = TRUE, 
                          flush = TRUE)
  
  # Remove individuals with NA values for coordinates
  data <- data[!is.na(data$Ausplot_X) & !is.na(data$Ausplot_Y), ]
  
  # Filter only trees which are alive/dead depending on what we're doing
  data <- data[data$Tree_Condition == condition, ]
  
  # Plots we will be considering
  plots <- unique(data$Site_Name)[seq_len(to_consider)]
  data <- data[data$Site_Name %in% plots, ]
  
  # Remove some of the plants to help with model stability
  nr <- nrow(data)
  if(thinning_probability > 0) {
    data <- data[-base::sample(seq_len(nr), thinning_probability * nr), ]
  }
  
  # Jitter 
  data$Ausplot_X <- rnorm(length(data$Ausplot_X), mean = data$Ausplot_X, sd = jitter)
  data$Ausplot_Y <- rnorm(length(data$Ausplot_Y), mean = data$Ausplot_Y, sd = jitter)
  
  # Group species with small abundance together into a new category
  species <- factor(data$Genus_Species)
  species_abundance <- sapply(levels(species), function(sp) {
    sum(species == sp)
  })
  data$Genus_Species[species_abundance[data$Genus_Species] < overall_minimum_abundance] <- "Rare species"
  
  # Put everything into Configuration format
  configurations <- setNames(lapply(plots, function(plot) {
    d <- data[data$Site_Name == plot, ]
    types <- factor(d$Genus_Species)
    types_names <- levels(types)
    types_names[sapply(types_names, function(sp) sum(types == sp)) < plot_minimum_abundance] <- "Rare species"
    levels(types) <- types_names
    if(use_marks) {
      ppjsdm::Configuration(x = d$Ausplot_X, y = d$Ausplot_Y, types = types, marks = d$Diameter / 100)
    } else {
      ppjsdm::Configuration(x = d$Ausplot_X, y = d$Ausplot_Y, types = types)
    }
  }), nm = plots)
  
  # Construct the corresponding windows
  windows <- setNames(lapply(plots, function(plot) {
    if(plot == "Supersite") { # Supersite with longer x-range
      ppjsdm::Rectangle_window(c(0, long_square_width), c(0, square_width))
    } else {
      ppjsdm::Rectangle_window(c(0, square_width), c(0, square_width))
    }
  }), nm = plots)
  
  # Convert to AU coordinate system
  require(rgdal)
  d <- data.frame(lon = data$Longitude, lat = data$Latitude)
  coordinates(d) <- c("lon", "lat")
  proj4string(d) <- CRS("+init=epsg:4326") # WGS 84
  CRS_new <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs") # GDA2020
  d <- sp::spTransform(d, CRS_new)
  data$GDA2020_X <- d@coords[, 1]
  data$GDA2020_Y <- d@coords[, 2]
  
  # Construct complete configuration
  xs <- data$GDA2020_X - min(data$GDA2020_X) + data$Ausplot_X
  ys <- data$GDA2020_Y - min(data$GDA2020_Y) + data$Ausplot_Y
  configuration <- if(use_marks) {
    ppjsdm::Configuration(x = xs, y = ys, types = factor(factor(data$Genus_Species)), marks = data$Diameter / 100)
  } else {
    ppjsdm::Configuration(x = xs, y = ys, types = factor(factor(data$Genus_Species)))
  }
  
  xs <- unique(data$GDA2020_X - min(data$GDA2020_X))
  ys <- unique(data$GDA2020_Y - min(data$GDA2020_Y))
  window <- ppjsdm::Rectangle_window_union(lapply(seq_len(length(xs)), function(n) {
    if(xs[n] == 471794) { # Supersite with longer x-range
      c(xs[n], long_square_width + xs[n])
    } else {
      c(xs[n], square_width + xs[n])
    }
  }), lapply(seq_len(length(xs)), function(n) {
    c(ys[n], square_width + ys[n])
  }))
  
  file_path <- file.path(tempdir(), "sv.RData")
  save(list = c("configuration",
                "configurations",
                "window",
                "windows"), file = file_path, envir = environment())
  file_path
}