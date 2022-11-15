suppressMessages(library(ggplot2))
suppressMessages(library(ppjsdm))
suppressMessages(library(rgdal))
suppressMessages(library(spatstat))

################################################################################
### The code starts by plotting the coefficients by region, while at the same time
### constructing a data.frame containing all coefficients to further analysis,
### as well as a list containing all the fits. This latter object is used further down
### to plot the potentials and the papangelou conditional intensities.
################################################################################

source("R/ausplots.R")

short <- 10 # Short-range distances
medium <- 40 # Medium-range distances
# medium <- 150 # Use this to put more weight on between plot interactions
long <- 120 # Long-range distances
nspecies <- 7 # This is the number of species we target on each plot, putting all others in the "Other species" class

no_medium <- FALSE # Set this to TRUE if you want there to not be ANY medium-range interactions

save_dir <- "~/AusplotsPicsSaturation50MediumWithinPlot/"

dir.create(save_dir, showWarnings = FALSE)

R_dataframe_to_latex <- function(df) {
  df <- as.data.frame(df)
  str <- paste0("\\begin{tabular}{|", paste0(rep("c", ncol(df)), collapse = "|"), "|}\n\\hline\n")
  for(i in seq_len(nrow(df))) {
    for(j in seq_len(ncol(df))) {
      if(j < ncol(df)) {
        str <- paste0(str, df[i, j], " & ")
      } else {
        str <- paste0(str, df[i, j], "\\\\")
      }
    }
    str <- paste0(str, "\n\\hline\n")
  }
  str <- paste0(str, "\\end{tabular}")
  cat(str)
}

update_dataframe <- function(dat,
                             estimates,
                             low,
                             high,
                             labels,
                             region,
                             range) {
  species1 <- ifelse(grepl("↔", labels, fixed = TRUE),
                     gsub("(.*) ↔ (.*)", "\\1", labels),
                     labels)
  species2 <- ifelse(grepl("↔", labels, fixed = TRUE),
                     gsub("(.*) ↔ (.*)", "\\2", labels),
                     labels)
  dat <- rbind(dat, data.frame(species1 = species1, species2 = species2, 
                               estimates = estimates, low = low, high = high, 
                               region = region, range = range))
  dat[!duplicated(dat), ]
}

plot_coef <- function(title,
                      estimates,
                      CI95_lo,
                      CI95_hi,
                      labels = names(CI95_lo),
                      base_size = 20,
                      only_significant = FALSE,
                      range = trunc(c(min(CI95_lo), max(CI95_hi)) * 10) / 10,
                      emphasize_zero = TRUE) {
  nc <- length(estimates)
  x <- factor(1:nc)
  levels(x) <- labels
  df <- data.frame(x = x, E = estimates, L = CI95_lo, U = CI95_hi, 
                   is_significant = CI95_lo > 0 | CI95_hi < 0)
  
  if(only_significant) {
    df <- df[df$L > 0 | df$U < 0, ]
  }
  
  g <- ggplot(df, aes(y = x, x = E))
  
  if(emphasize_zero) {
    g <- g + geom_vline(xintercept = 0, color = "red", size = 2) # Vertical line at 0
  }
  
  # Make sure the plot does not have too many ticks on the x axis
  # breaks <- seq(from = range[1], to = range[2], by = 0.1)
  # if(length(breaks) > 10) {
  #   breaks <- seq(from = range[1], to = range[2], by = 0.25) 
  #   if(length(breaks) > 10) {
  #     breaks <- seq(from = range[1], to = range[2], by = 0.5)
  #     if(length(breaks) > 10) {
  #       breaks <- seq(from = range[1], to = range[2], by = 1)
  #       minor_breaks <- seq(from = range[1], to = range[2], by = 0.5)
  #     } else {
  #       minor_breaks <- seq(from = range[1], to = range[2], by = 0.25)
  #     }
  #   } else {
  #     minor_breaks <- seq(from = range[1], to = range[2], by = 0.125)
  #   }
  # } else {
  #   minor_breaks <- seq(from = range[1], to = range[2], by = 0.05)
  # }
  
  g <- g +
    geom_point(aes(colour = is_significant), size = 5) + # Size of point estimates
    geom_errorbarh(aes(xmax = U, xmin = L, colour = is_significant), size = 1.5, height = 0.1) + # Error bars
    xlab(NULL) + # Remove x labels
    ylab(NULL) + # Remove y labels
    ggtitle(title) + # Title
    # scale_x_continuous(minor_breaks = minor_breaks, 
    #                    breaks = breaks) +
    scale_colour_manual(values = c("grey", "black"), breaks = c(FALSE, TRUE)) +
    guides(colour = "none") +
    theme_minimal(base_size = base_size) + # Theme
    theme(panel.grid.major = element_line(size = 2),
          axis.text = element_text(colour = "black", face = "bold", size = 0.8 * base_size),
          plot.title = element_text(size = 0.8 * base_size, face = "bold", hjust = 0.5)) # Title size & position
  g
}

plot_for_region <- function(region,
                            show = TRUE,
                            split_threshold = NA,  # Threshold at which to split Eucalypts into thick/thin
                            use_plot_long = FALSE, # Set this to true to use the plot-specific long-range distances
                            dat = dat) { # Data frame to be updated when executing this function
  set.seed(1) # Reproducibility wrt dummy point distribution
  
  # Code below loads up the data successive times with varying abundance, until
  # we have reached the specified number of species. Definitely not the prettiest, but
  # works for the time being...
  abundance <- 30
  while(TRUE) {
    load(ausplots(plots_to_consider = region$plots,
                  overall_minimum_abundance = abundance, 
                  use_marks = FALSE,
                  split_threshold = split_threshold,
                  plot_minimum_abundance = 0))
    if(nlevels(configuration$types) <= nspecies) {
      break
    } else {
      abundance <- abundance + 10
    }
  }
  
  if(no_medium) {
    medium <- long <- 0
  }
  
  if(use_plot_long) {
    long <- region$long
  }
  
  ntypes <- nlevels(configuration$types)
  tm <- Sys.time()
  fit <- gibbsm(configuration,
                window = window,
                model = "square_exponential", # Short-range interaction
                medium_range_model = "tanh", # Medium-range interaction
                short_range = matrix(short, ntypes, ntypes), # Matrix of short-range interaction distances
                medium_range = matrix(medium, ntypes, ntypes), # Matrix of medium-range interaction distances
                long_range = matrix(long, ntypes, ntypes), # Matrix of long-range interaction distances
                saturation = 50, # Saturation parameter
                fitting_package = "glm", # Fitting package 
                use_regularization = FALSE, # Do we want to run a regularised fit? Useful if fitting psckage is glmnet
                min_dummy = 1e5, # Force at least this many dummy points per species
                max_dummy = 1e5, # Force at maximum this many dummy points per species
                nthreads = 4, # Number of CPU cores to use
                dummy_distribution = "binomial") # Distribution of dummy points
  if(any(abs(fit$coefficients_vector) > 1e5)) {
    print(coef(fit))
    stop("Absurd fitted coefficients, there is probably an error in the model specification.")
  }
  summary_fit <- summary(fit)
  print(Sys.time() - tm)
  
  if(show) {
    print(fit$coefficients$alpha[[1]])
    print(fit$coefficients$gamma)
    print(summary_fit)
    print(abundance)
  }
  
  # Plot the coefficients
  nspecies <- nlevels(configuration$types)
  parameters <- fit$coefficients
  alpha <- parameters$alpha[[1]]
  gamma <- parameters$gamma
  
  ci_fit <- cbind(summary_fit$coefficients$CI95_lo, summary_fit$coefficients$CI95_hi)
  rownames(ci_fit) <- rownames(summary_fit$coefficients)
  
  alpha_diagonal_str <- paste0("alpha1_", 1:nspecies, "_", 1:nspecies)
  gamma_diagonal_str <- paste0("gamma_", 1:nspecies, "_", 1:nspecies)
  
  alpha_off_diagonal_str <- c()
  gamma_off_diagonal_str <- c()
  interaction_description <- c()
  for(i in seq_len(nspecies - 1)) {
    for(j in (i + 1):nspecies) {
      alpha_off_diagonal_str <- c(alpha_off_diagonal_str, paste0("alpha1_", i, "_", j))
      gamma_off_diagonal_str <- c(gamma_off_diagonal_str, paste0("gamma_", i, "_", j))
      interaction_description <- c(interaction_description, paste0(levels(ppjsdm::types(configuration))[i],
                                                                   " ↔ ",
                                                                   levels(ppjsdm::types(configuration))[j]))
    }
  }
  
  estimates_alpha_diagonal <- summary_fit$coefficients[rownames(summary_fit$coefficients) %in% alpha_diagonal_str, 1]
  estimates_gamma_diagonal <- summary_fit$coefficients[rownames(summary_fit$coefficients) %in% gamma_diagonal_str, 1]
  i <- order(estimates_alpha_diagonal)
  im <- order(estimates_gamma_diagonal)
  estimates_alpha_diagonal <- estimates_alpha_diagonal[i]
  estimates_gamma_diagonal <- estimates_gamma_diagonal[im]
  
  CI95_lo_alpha_diagonal <- ci_fit[rownames(ci_fit) %in% alpha_diagonal_str, 1][i]
  CI95_hi_alpha_diagonal <- ci_fit[rownames(ci_fit) %in% alpha_diagonal_str, 2][i]
  CI95_lo_gamma_diagonal <- ci_fit[rownames(ci_fit) %in% gamma_diagonal_str, 1][im]
  CI95_hi_gamma_diagonal <- ci_fit[rownames(ci_fit) %in% gamma_diagonal_str, 2][im]
  
  estimates_alpha_off_diagonal <- summary_fit$coefficients[startsWith(rownames(summary_fit$coefficients), "alpha1_") & !(rownames(summary_fit$coefficients) %in% alpha_diagonal_str), 1]
  estimates_gamma_off_diagonal <- summary_fit$coefficients[startsWith(rownames(summary_fit$coefficients), "gamma_") & !(rownames(summary_fit$coefficients) %in% gamma_diagonal_str), 1]
  
  j <- order(estimates_alpha_off_diagonal)
  jm <- order(estimates_gamma_off_diagonal)
  estimates_alpha_off_diagonal <- estimates_alpha_off_diagonal[j]
  estimates_gamma_off_diagonal <- estimates_gamma_off_diagonal[jm]
  
  CI95_lo_alpha_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "alpha1_") & !(rownames(ci_fit) %in% alpha_diagonal_str), 1][j]
  CI95_hi_alpha_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "alpha1_") & !(rownames(ci_fit) %in% alpha_diagonal_str), 2][j]
  CI95_lo_gamma_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "gamma_") & !(rownames(ci_fit) %in% gamma_diagonal_str), 1][jm]
  CI95_hi_gamma_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "gamma_") & !(rownames(ci_fit) %in% gamma_diagonal_str), 2][jm]
  
  estimates_beta_intercept <- summary_fit$coefficients[startsWith(rownames(summary_fit$coefficients), "beta0"), 1][i]
  CI95_lo_beta_intercept <- ci_fit[startsWith(rownames(ci_fit), "beta0"), 1][i]
  CI95_hi_beta_intercept <- ci_fit[startsWith(rownames(ci_fit), "beta0"), 2][i]
  
  dat <- update_dataframe(dat = dat,
                          estimates = estimates_alpha_diagonal,
                          low = CI95_lo_alpha_diagonal,
                          high = CI95_hi_alpha_diagonal,
                          labels = levels(ppjsdm::types(configuration))[i],
                          region = region$name,
                          range = "Short")
  dat <- update_dataframe(dat = dat,
                          estimates = estimates_alpha_off_diagonal,
                          low = CI95_lo_alpha_off_diagonal,
                          high = CI95_hi_alpha_off_diagonal,
                          labels = interaction_description[j],
                          region = region$name,
                          range = "Short")
  if(!all(estimates_gamma_diagonal == 0)) {
    dat <- update_dataframe(dat = dat,
                            estimates = estimates_gamma_diagonal,
                            low = CI95_lo_gamma_diagonal,
                            high = CI95_hi_gamma_diagonal,
                            labels = levels(ppjsdm::types(configuration))[im],
                            region = region$name,
                            range = "Medium")
  }
  if(!all(estimates_gamma_off_diagonal == 0)) {
    dat <- update_dataframe(dat = dat,
                            estimates = estimates_gamma_off_diagonal,
                            low = CI95_lo_gamma_off_diagonal,
                            high = CI95_hi_gamma_off_diagonal,
                            labels = interaction_description[jm],
                            region = region$name,
                            range = "Medium")
  }
  
  png(file = paste0(save_dir, region$name, "_short_within.png"), bg = "white", width = 1600, height = 1000)
  print(plot_coef(title = paste0("Short range interaction within species in ", region$name),
                  estimates = estimates_alpha_diagonal,
                  CI95_lo = CI95_lo_alpha_diagonal,
                  CI95_hi = CI95_hi_alpha_diagonal,
                  labels = levels(ppjsdm::types(configuration))[i],
                  base_size = 35))
  dev.off()
  png(file = paste0(save_dir, region$name, "_short_between.png"), bg = "white", width = 1600, height = 1000)
  print(plot_coef(title = paste0("Short range interaction between species in ", region$name),
                  estimates = estimates_alpha_off_diagonal,
                  CI95_lo = CI95_lo_alpha_off_diagonal,
                  CI95_hi = CI95_hi_alpha_off_diagonal,
                  labels = interaction_description[j],
                  only_significant = FALSE,
                  base_size = 35))
  dev.off()
  if(!all(estimates_gamma_diagonal == 0)) {
    png(file = paste0(save_dir, region$name, "_medium_within.png"), bg = "white", width = 1600, height = 1000)
    print(plot_coef(title = paste0("Medium range interaction within species in ", region$name),
                    estimates = estimates_gamma_diagonal,
                    CI95_lo = CI95_lo_gamma_diagonal,
                    CI95_hi = CI95_hi_gamma_diagonal,
                    labels = levels(ppjsdm::types(configuration))[im],
                    base_size = 35))
    dev.off()
  }
  if(!all(estimates_gamma_off_diagonal == 0)) {
    png(file = paste0(save_dir, region$name, "_medium_between.png"), bg = "white", width = 1600, height = 1000)
    print(plot_coef(title = paste0("Medium range interaction between species in ", region$name),
                    estimates = estimates_gamma_off_diagonal,
                    CI95_lo = CI95_lo_gamma_off_diagonal,
                    CI95_hi = CI95_hi_gamma_off_diagonal,
                    labels = interaction_description[jm],
                    only_significant = FALSE,
                    base_size = 35))
    dev.off()
  }
  
  list(dat = dat,
       fit = fit,
       abundance = abundance)
}

fits <- list() # Initial list that will contain all the fit objects

result <- plot_for_region(list(name = "NVic",
                               plots = c("ANU101", "Ada Tree", "ANU363", "ANU589", "HardyCreek"),
                               long = 1000),
                          split_threshold = 80,
                          dat = data.frame())
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$NVic <- result$fit

result <- plot_for_region(list(name = "SVic",
                               plots = c("Weeaproinah", "Turtons", "Lardner"),
                               long = 1000),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$SVic <- result$fit

result <- plot_for_region(list(name = "TASNoDelega",
                               plots = c("Bird", "Supersite", "Weld", "ZigZag", "BlackRiver", "BondTier", 
                                         "Flowerdale", "Dip"),
                               long = 5000),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$TASNoDelega <- result$fit

result <- plot_for_region(list(name = "TASDelega",
                               plots = c("NorthStyx", "MtField", "MtMaurice", "BenRidge", "Mackenzie", "Caveside"),
                               abundance = 200,
                               long = 5000),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$TASDelega <- result$fit

result <- plot_for_region(list(name = "SENSW",
                               plots = c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"),
                               long = 1000),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$SENSW <- result$fit

result <- plot_for_region(list(name = "NNSW",
                               plots = c("MinesRd", "Tinebank", "Lorne", "BirdTree", "A-Tree", 
                                         "BlackBull", "Bruxner", "Osullivans"),
                               long = 1000),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$NNSW <- result$fit

result <- plot_for_region(list(name = "QLD",
                               plots = c("Herberton", "Lamb Range", "Baldy", "Koombooloomba"),
                               long = 500),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$QLD <- result$fit

result <- plot_for_region(list(name = "SWA",
                               plots = c("Frankland", "Clare", "Giants", "Dawson"),
                               long = 500),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$SWA <- result$fit

result <- plot_for_region(list(name = "NWA",
                               plots = c("Collins", "Sutton", "Warren", "Carey", "Dombakup"),
                               long = 1000),
                          split_threshold = 80,
                          dat = result$dat)
result$fit$complete <- NULL # Used to reduce RAM usage, this object is quite large
fits$NWA <- result$fit

dat <- result$dat # Data.frame containing all the estimated coefficients

################################################################################
### At this point in the script, we have a dat object which contains all the fitted
### coefficients, and a fits object which contains the fit objects for each region.
### Below, we show how to extract useful information from these objects.
################################################################################

# Plotting the potentials of one of the fits
region <- "SVic"
print(names(coef(fits[[region]])$beta0)) # This is a quick way to get a vector of all species in this fit
sp <- "Big Eucalyptus regnans" # Choose one of the names below for the analysis
plot(potentials(fits[[region]], sp, sp)) # Plot the potential of this species with itself, for this fit
plot(potentials(fits[[region]], "Big Eucalyptus regnans", "Small Eucalyptus regnans")) # Another interesting potential

# Plotting the Papangelou conditional intensity on one of the plots, given one of our fits
region <- "SVic"
print(names(coef(fits[[region]])$beta0)) # This is a quick way to get a vector of all species in this fit
sp <- "Small Eucalyptus regnans"
plot <- "Turtons" # Turtons is one of the plots in the South Vic region
# IMPORTANT: In the code below, it is crucial to use the SAME list of plots
# as in the target region, AND the overall_minimum_abundance used to generate the fits.
# The reason for this is that if we change the parameters, then which species
# fall into the "Other species" class will change, and thus the fit cannot be applied.
load(ausplots(plots = c("Weeaproinah", "Turtons", "Lardner"), 
              plot_minimum_abundance = 0, 
              split_threshold = 80,
              overall_minimum_abundance = 100))
# The plot below shows the conditional intensity of the target species `sp`,
# on a single one of the plots, conditioning on all species other than `sp`.
# This allows us, by inspection, to check how well the model is performing on this
# plot, for this species.
# DISCLAIMER: On other plots, and for other species, the results do not look as good.
plot_papangelou(fits[[region]],
                window = windows[[plot]],
                configuration = configurations[[plot]],
                type = sp,
                use_log = TRUE,
                drop_type_from_configuration = TRUE,
                show_only_type = TRUE,
                grid_steps = c(500, 500))

# The code above can be executed for all species on a given plot, and we may, e.g., compute
# AUCs to check how well the model is performing on this plot.
region <- "SVic"
print(names(coef(fits[[region]])$beta0)) # This is a quick way to get a vector of all species in this fit
sp <- "Small Eucalyptus regnans"
plot <- "Turtons" # Turtons is one of the plots in the South Vic region
papangelous <- setNames(lapply(names(coef(fits[[region]])$beta0), function(sp) {
  plot_papangelou(fits[[region]],
                  window = windows[[plot]],
                  configuration = configurations[[plot]],
                  type = sp,
                  use_log = FALSE,
                  drop_type_from_configuration = TRUE,
                  return_papangelou = TRUE,
                  grid_steps = c(100, 100))
}), nm = names(coef(fits[[region]])$beta0))

aucs <- sapply(names(coef(fits[[region]])$beta0), function(sp) {
  intensity <- papangelous[[sp]]
  conf <- ppp(x = configurations[[plot]]$x[configurations[[plot]]$types == sp],
              y = configurations[[plot]]$y[configurations[[plot]]$types == sp],
              window = owin(windows[[plot]]$x_range, windows[[plot]]$y_range))
  auc(X = conf, covariate = intensity)
})

# Performance (measured as conditional AUC) on this plot:
print(aucs)

# Example plot by using the data.frame containing all the coefficients
eucalypts <- dat[dat$range == "Short" & # Short-range
                   dat$species1 == dat$species2 & # Within species
                   grepl("Eucalyptus", dat$species1, fixed = TRUE), ] # Is a Eucalypt
eucalypts$is_big <- grepl("Big", eucalypts$species1, fixed = TRUE)
ggplot() + 
  geom_point(aes(x = eucalypts$estimates[eucalypts$is_big],
                 y = eucalypts$estimates[!eucalypts$is_big])) +
  theme_minimal() +
  geom_vline(aes(xintercept = 0), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "red") +
  xlab("Big Eucalyptus within-species alpha") +
  ylab("Small Eucalyptus within-species alpha")

# Example data that can be extracted from the df
# vector below contains all Big Eucalyptus vs Small Eucalyptus short-range interactions
big_small_short <- sort(dat$estimates[startsWith(as.character(dat$species1), "Big") & startsWith(as.character(dat$species2), "Small") & dat$range == "Short"])
ggplot() +
  geom_histogram(aes(x = big_small_short), binwidth = 0.05) +
  theme_minimal() +
  geom_vline(aes(xintercept = mean(big_small_short)), colour = "red") +
  xlab("Big Eucalyptus VS Small Eucalyptus short-range interaction coefficient")

################################################################################
### Code below generates complementary analysis of the Bruxner plot used in a presentation
################################################################################

# Bruxner plot
load(ausplots(plots = "Bruxner", plot_minimum_abundance = 15, overall_minimum_abundance = 15))
window <- windows$Bruxner

ntypes <- nlevels(configuration$types)
tm <- Sys.time()
fit <- gibbsm(configuration,
              window = window,
              model = "square_exponential", # Short-range interaction
              medium_range_model = "tanh", # Medium-range interaction
              short_range = matrix(5, ntypes, ntypes), # Matrix of short-range interaction distances
              medium_range = matrix(10, ntypes, ntypes), # Matrix of short-range interaction distances
              long_range = matrix(20, ntypes, ntypes), # Matrix of short-range interaction distances
              saturation = 10, # Saturation parameter
              fitting_package = "glm", # Fitting package 
              min_dummy = 1e4, # Force at least this number of dummy points per species
              max_dummy = 1e4, # Force at maximum this number of dummy points per species
              nthreads = 4, # Number of threads to use
              dummy_distribution = "stratified") # Distribution of dummy points
Sys.time() - tm

summary_fit <- summary(fit)

# Plot the coefficients
nspecies <- nlevels(configuration$types)
parameters <- fit$coefficients
alpha <- parameters$alpha[[1]]
gamma <- parameters$gamma

ci_fit <- cbind(summary_fit$coefficients$CI95_lo, summary_fit$coefficients$CI95_hi)
rownames(ci_fit) <- rownames(summary_fit$coefficients)

alpha_diagonal_str <- paste0("alpha1_", 1:nspecies, "_", 1:nspecies)
gamma_diagonal_str <- paste0("gamma_", 1:nspecies, "_", 1:nspecies)

alpha_off_diagonal_str <- c()
gamma_off_diagonal_str <- c()
interaction_description <- c()
for(i in seq_len(nspecies - 1)) {
  for(j in (i + 1):nspecies) {
    alpha_off_diagonal_str <- c(alpha_off_diagonal_str, paste0("alpha1_", i, "_", j))
    gamma_off_diagonal_str <- c(gamma_off_diagonal_str, paste0("gamma_", i, "_", j))
    interaction_description <- c(interaction_description, paste0(levels(ppjsdm::types(configuration))[i],
                                                                 " ↔ ",
                                                                 levels(ppjsdm::types(configuration))[j]))
  }
}

estimates_alpha_diagonal <- summary_fit$coefficients[rownames(summary_fit$coefficients) %in% alpha_diagonal_str, 1]
estimates_gamma_diagonal <- summary_fit$coefficients[rownames(summary_fit$coefficients) %in% gamma_diagonal_str, 1]
i <- order(estimates_alpha_diagonal)
im <- order(estimates_gamma_diagonal)
estimates_alpha_diagonal <- estimates_alpha_diagonal[i]
estimates_gamma_diagonal <- estimates_gamma_diagonal[im]

CI95_lo_alpha_diagonal <- ci_fit[rownames(ci_fit) %in% alpha_diagonal_str, 1][i]
CI95_hi_alpha_diagonal <- ci_fit[rownames(ci_fit) %in% alpha_diagonal_str, 2][i]
CI95_lo_gamma_diagonal <- ci_fit[rownames(ci_fit) %in% gamma_diagonal_str, 1][im]
CI95_hi_gamma_diagonal <- ci_fit[rownames(ci_fit) %in% gamma_diagonal_str, 2][im]

estimates_alpha_off_diagonal <- summary_fit$coefficients[startsWith(rownames(summary_fit$coefficients), "alpha1_") & !(rownames(summary_fit$coefficients) %in% alpha_diagonal_str), 1]
estimates_gamma_off_diagonal <- summary_fit$coefficients[startsWith(rownames(summary_fit$coefficients), "gamma_") & !(rownames(summary_fit$coefficients) %in% gamma_diagonal_str), 1]

j <- order(estimates_alpha_off_diagonal)
jm <- order(estimates_gamma_off_diagonal)
estimates_alpha_off_diagonal <- estimates_alpha_off_diagonal[j]
estimates_gamma_off_diagonal <- estimates_gamma_off_diagonal[jm]

CI95_lo_alpha_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "alpha1_") & !(rownames(ci_fit) %in% alpha_diagonal_str), 1][j]
CI95_hi_alpha_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "alpha1_") & !(rownames(ci_fit) %in% alpha_diagonal_str), 2][j]
CI95_lo_gamma_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "gamma_") & !(rownames(ci_fit) %in% gamma_diagonal_str), 1][jm]
CI95_hi_gamma_off_diagonal <- ci_fit[startsWith(rownames(ci_fit), "gamma_") & !(rownames(ci_fit) %in% gamma_diagonal_str), 2][jm]

estimates_beta_intercept <- summary_fit$coefficients[startsWith(rownames(summary_fit$coefficients), "beta0"), 1][i]
CI95_lo_beta_intercept <- ci_fit[startsWith(rownames(ci_fit), "beta0"), 1][i]
CI95_hi_beta_intercept <- ci_fit[startsWith(rownames(ci_fit), "beta0"), 2][i]

png(file = paste0(save_dir, "Bruxner_short_within.png"), bg = "white", width = 1600, height = 1000)
print(plot_coef(title = paste0("Short range interaction within species in Bruxner"),
                estimates = estimates_alpha_diagonal,
                CI95_lo = CI95_lo_alpha_diagonal,
                CI95_hi = CI95_hi_alpha_diagonal,
                labels = levels(ppjsdm::types(configuration))[i],
                base_size = 35))
dev.off()
png(file = paste0(save_dir, "Bruxner_short_between.png"), bg = "white", width = 1600, height = 1000)
print(plot_coef(title = paste0("Short range interaction between species in Bruxner"),
                estimates = estimates_alpha_off_diagonal,
                CI95_lo = CI95_lo_alpha_off_diagonal,
                CI95_hi = CI95_hi_alpha_off_diagonal,
                labels = interaction_description[j],
                only_significant = TRUE,
                base_size = 35))
dev.off()
if(!all(estimates_gamma_diagonal == 0)) {
  png(file = paste0(save_dir, "Bruxner_medium_within.png"), bg = "white", width = 1600, height = 1000)
  print(plot_coef(title = paste0("Medium range interaction within species in Bruxner"),
                  estimates = estimates_gamma_diagonal,
                  CI95_lo = CI95_lo_gamma_diagonal,
                  CI95_hi = CI95_hi_gamma_diagonal,
                  labels = levels(ppjsdm::types(configuration))[im],
                  base_size = 35))
  dev.off()
}
if(!all(estimates_gamma_off_diagonal == 0)) {
  png(file = paste0(save_dir, "Bruxner_medium_between.png"), bg = "white", width = 1600, height = 1000)
  print(plot_coef(title = paste0("Medium range interaction between species in Bruxner"),
                  estimates = estimates_gamma_off_diagonal,
                  CI95_lo = CI95_lo_gamma_off_diagonal,
                  CI95_hi = CI95_hi_gamma_off_diagonal,
                  labels = interaction_description[jm],
                  only_significant = TRUE,
                  base_size = 35))
  dev.off()
}

png(file = paste0(save_dir, "Bruxner_Allocasuarina.png"), bg = "white", width = 500, height = 300)
plot_papangelou(fit,
                window = window,
                configuration = configuration,
                type = "Allocasuarina torulosa",
                use_log = TRUE,
                drop_type_from_configuration = TRUE,
                show_only_type = TRUE,
                grid_steps = c(100, 100))
dev.off()

png(file = paste0(save_dir, "Bruxner_grandis.png"), bg = "white", width = 500, height = 300)
plot_papangelou(fit,
                window = window,
                configuration = configuration,
                type = "Eucalyptus grandis",
                use_log = TRUE,
                drop_type_from_configuration = TRUE,
                show_only_type = TRUE,
                grid_steps = c(100, 100))
dev.off()

species <- levels(configuration$types)

papangelous <- setNames(lapply(species, function(sp) {
  plot_papangelou(fit,
                  window = window,
                  configuration = configuration,
                  type = sp,
                  use_log = FALSE,
                  drop_type_from_configuration = TRUE,
                  return_papangelou = TRUE,
                  grid_steps = c(100, 100))
}), nm = species)

aucs <- sapply(species, function(sp) {
  intensity <- papangelous[[sp]]
  conf <- ppp(x = configuration$x[configuration$types == sp],
              y = configuration$y[configuration$types == sp],
              window = owin(window$x_range, window$y_range))
  auc(X = conf, covariate = intensity)
})

R_dataframe_to_latex(data.frame(Species = names(sort(aucs)),
                                AUC = round(sort(aucs), digits = 2)))

sp <- "Allocasuarina torulosa"
plot_papangelou(fit,
                window = window,
                configuration = configuration,
                type = sp,
                use_log = TRUE,
                drop_type_from_configuration = TRUE,
                show_only_type = TRUE,
                grid_steps = c(500, 500))

sp <- "Lophostemon sp."
plot_papangelou(fit,
                window = window,
                configuration = configuration,
                type = sp,
                use_log = TRUE,
                drop_type_from_configuration = TRUE,
                show_only_type = TRUE,
                grid_steps = c(500, 500))

sp <- "Eucalyptus grandis"
plot_papangelou(fit,
                window = window,
                configuration = configuration,
                type = sp,
                use_log = TRUE,
                drop_type_from_configuration = TRUE,
                show_only_type = TRUE,
                grid_steps = c(500, 500))

sp <- "Allocasuarina torulosa"
R_dataframe_to_latex(data.frame(Species = names(sort(fit$coefficients$alpha[[1]][sp, ])),
                                Alpha = round(sort(fit$coefficients$alpha[[1]][sp, ]), digits = 2),
                                Gamma = round(fit$coefficients$gamma[sp, ][names(sort(fit$coefficients$alpha[[1]][sp, ]))], digits = 2)))
plot_papangelou(fit,
                window = window,
                configuration = configuration[unique(c(sp, 
                                                       names(sort(fit$coefficients$alpha[[1]][sp, ]))[1],
                                                       names(sort(fit$coefficients$alpha[[1]][sp, ], decreasing = TRUE))[1],
                                                       names(sort(fit$coefficients$gamma[sp, ]))[1],
                                                       names(sort(fit$coefficients$gamma[sp, ], decreasing = TRUE))[1]))],
                type = sp,
                use_log = TRUE,
                drop_type_from_configuration = TRUE,
                show_only_type = FALSE,
                grid_steps = c(500, 500))
