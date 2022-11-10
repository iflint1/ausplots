ausplots <- function(overall_minimum_abundance = 100, # Minimum overall number of individuals for the species not to be considered 'rare'
                     plot_minimum_abundance = 5, # Minimum number of individuals in a given plot for the species not to be considered 'rare'
                     thinning_probability = 0, # Remove this percent of individuals
                     condition = "Alive", # Select only live or dead individuals
                     use_marks = FALSE, # Use individuals' diameters as marks?
                     jitter = 1e-5, # Jitter locations with this std error
                     download_path = getwd(), # Where should we download the Ausplots CSV?
                     only_species_with_per_plot = 0, # Set this to >0 value to force all considered species to have a min amount of individuals on each plot
                     plots_to_consider) { # Which plots to consider?
  # Constants related to the dataset
  square_width <- 100
  long_square_width <- 160
  
  # If the data does not yet exist, download it
  if(!file.exists(file.path(download_path, "ausplots.csv"))) {
    download.file(url = "https://shared.tern.org.au/attachment/0e503109-2fb6-4182-969f-2d570abdbabd/data_large_tree_survey.csv",
                  destfile = file.path(download_path, "ausplots.csv"))
  }
  
  # Load the required data
  data <- utils::read.csv(file = file.path(download_path, "ausplots.csv"), 
                          stringsAsFactors = FALSE, 
                          allowEscapes = TRUE, 
                          flush = TRUE)
  
  # Remove individuals with NA values for coordinates
  data <- data[!is.na(data$Ausplot_X) & !is.na(data$Ausplot_Y), ]
  
  # Filter only trees which are alive/dead depending on what we're doing
  data <- data[data$Tree_Condition == condition, ]
  
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
  min_x <- min(data$GDA2020_X)
  min_y <- min(data$GDA2020_Y)
  
  # Plots we will be considering
  plots <- unique(data$Site_Name)
  
  # Fill in default value for plot_to_consider
  if(!missing(plots_to_consider)) {
    plots <- plots[plots %in% plots_to_consider]
  }
  
  # Subset data to the plots we are interested in
  data <- data[data$Site_Name %in% plots, ]
  
  # Remove some of the plants to help with model stability
  nr <- nrow(data)
  if(thinning_probability > 0) {
    data <- data[-base::sample(seq_len(nr), thinning_probability * nr), ]
  }
  
  # Jitter 
  data$Ausplot_X <- rnorm(length(data$Ausplot_X), mean = data$Ausplot_X, sd = jitter)
  data$Ausplot_Y <- rnorm(length(data$Ausplot_Y), mean = data$Ausplot_Y, sd = jitter)
  
  # Restrict to proper range
  data$Ausplot_X[data$Site_Name == "Supersite"] <- pmax(pmin(data$Ausplot_X[data$Site_Name == "Supersite"], long_square_width), 0)
  data$Ausplot_X[data$Site_Name != "Supersite"] <- pmax(pmin(data$Ausplot_X[data$Site_Name != "Supersite"], square_width), 0)
  data$Ausplot_Y <- pmax(pmin(data$Ausplot_Y, square_width), 0)
  
  # Group species with small abundance together into a new category
  species <- factor(data$Genus_Species)
  species_abundance <- sapply(levels(species), function(sp) {
    sum(species == sp)
  })
  data$Genus_Species[species_abundance[data$Genus_Species] < overall_minimum_abundance] <- "Rare species"
  
  # Ensure that we only consider species with a sufficient number of individuals per plot
  new_names <- sapply(unique(data$Genus_Species), function(sp) {
    for(pl in plots) {
      if(sum(data$Genus_Species == sp & data$Site_Name == pl) < only_species_with_per_plot) {
        return("Restricted species")
      }
    }
    return(sp)
  })
  data$Genus_Species <- sapply(data$Genus_Species, function(sp) new_names[sp])
  
  # Construct complete configuration
  xs <- data$GDA2020_X - min_x + data$Ausplot_X
  ys <- data$GDA2020_Y - min_y + data$Ausplot_Y
  configuration <- if(use_marks) {
    xs <- xs[!is.na(data$Diameter)]
    ys <- ys[!is.na(data$Diameter)]
    ppjsdm::Configuration(x = xs, y = ys, types = factor(data$Genus_Species[!is.na(data$Diameter)]), marks = data$Diameter[!is.na(data$Diameter)] / 100)
  } else {
    ppjsdm::Configuration(x = xs, y = ys, types = factor(data$Genus_Species))
  }
  
  # Construct the big window
  xs <- unique(data$GDA2020_X - min_x)
  ys <- unique(data$GDA2020_Y - min_y)
  super_site_x <- which(xs == data$GDA2020_X[data$Site_Name == "Supersite"][1] - min_x)
  if(length(super_site_x) == 0) {
    super_site_x <- Inf
  }
  window <- ppjsdm::Rectangle_window_union(lapply(seq_len(length(xs)), function(n) {
    if(n == super_site_x) { # Supersite with longer x-range
      c(xs[n], long_square_width + xs[n])
    } else {
      c(xs[n], square_width + xs[n])
    }
  }), lapply(seq_len(length(xs)), function(n) {
    c(ys[n], square_width + ys[n])
  }))
  
  # Put everything into Configuration format
  configurations <- setNames(lapply(plots, function(plot) {
    d <- data[data$Site_Name == plot, ]
    types <- factor(d$Genus_Species)
    types_names <- levels(types)
    types_names[sapply(types_names, function(sp) sum(types == sp)) < plot_minimum_abundance] <- "Rare species"
    levels(types) <- types_names
    if(use_marks) {
      ppjsdm::Configuration(x = d$Ausplot_X + d$GDA2020_X - min_x, y = d$Ausplot_Y + d$GDA2020_Y - min_y, types = types, marks = d$Diameter / 100)
    } else {
      ppjsdm::Configuration(x = d$Ausplot_X + d$GDA2020_X - min_x, y = d$Ausplot_Y + d$GDA2020_Y - min_y, types = types)
    }
  }), nm = plots)
  
  # Construct the corresponding windows
  xs <- unique(data$GDA2020_X - min_x)
  ys <- unique(data$GDA2020_Y - min_y)
  super_site_x <- which(xs == data$GDA2020_X[data$Site_Name == "Supersite"][1] - min_x)
  if(length(super_site_x) == 0) {
    super_site_x <- Inf
  }
  windows <- setNames(lapply(seq_len(length(xs)), function(n) {
    if(n == super_site_x) { # Supersite with longer x-range
      ppjsdm::Rectangle_window(c(xs[n], long_square_width + xs[n]), c(ys[n], square_width + ys[n]))
    } else {
      ppjsdm::Rectangle_window(c(xs[n], square_width + xs[n]), c(ys[n], square_width + ys[n]))
    }
  }), nm = plots)
  
  file_path <- file.path(tempdir(), "sv.RData")
  save(list = c("configuration",
                "configurations",
                "window",
                "windows"), file = file_path, envir = environment())
  file_path
}
