suppressMessages(library(ggplot2))
suppressMessages(library(ppjsdm))
suppressMessages(library(rgdal))
suppressMessages(library(spatstat))

source("R/ausplots.R")

short <- 10 # Short-range distances
medium <- 0 # Medium-range distances, larger than the plot size
long <- 0 # Long-range distances, chosen to forbid the plots located far away from each other from "interacting" much

save_dir <- "~/AusplotsPics/"

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

plot_coef <- function(title,
                      estimates,
                      CI95_lo,
                      CI95_hi,
                      labels = names(CI95_lo),
                      base_size = 20,
                      only_significant = FALSE,
                      range = trunc(c(min(CI95_lo), max(CI95_hi)) * 2) / 2 + c(-0.5, 0.5),
                      emphasize_zero = TRUE) {
  nc <- length(estimates)
  x <- factor(1:nc)
  levels(x) <- labels
  df <- data.frame(x = x, E = estimates, L = CI95_lo, U = CI95_hi)
  
  if(only_significant) {
    df <- df[df$L > 0 | df$U < 0, ]
  }
  
  g <- ggplot(df, aes(y = x, x = E))
  
  if(emphasize_zero) {
    g <- g + geom_vline(xintercept = 0, color = "red", size = 2) # Vertical line at 0
  }
  
  g <- g +
    geom_point(size = 5) + # Size of point estimates
    geom_errorbarh(aes(xmax = U, xmin = L), size = 1.5, height = 0.1) + # Error bars
    xlab(NULL) + # Remove x labels
    ylab(NULL) + # Remove y labels
    ggtitle(title) + # Title
    scale_x_continuous(minor_breaks = seq(from = range[1], to = range[2], by = 0.125), breaks = seq(from = range[1], to = range[2], by = 0.25)) +
    theme_minimal(base_size = base_size) + # Theme
    theme(panel.grid.major = element_line(size = 2),
          axis.text = element_text(colour = "black", face = "bold", size = 0.8 * base_size),
          plot.title = element_text(size = 0.8 * base_size, face = "bold", hjust = 0.5)) # Title size & position
  g
}

plot_for_region <- function(region,
                            show = FALSE) {
  load(ausplots(plots_to_consider = region$plots,
                overall_minimum_abundance = region$abundance, 
                plot_minimum_abundance = 1))
  
  ntypes <- nlevels(configuration$types)
  tm <- Sys.time()
  fit <- gibbsm(configuration,
                window = window,
                model = "square_exponential", # Short-range interaction
                medium_range_model = "tanh", # Medium-range interaction
                short_range = matrix(short, ntypes, ntypes), # Matrix of short-range interaction distances
                medium_range = matrix(medium, ntypes, ntypes), # Matrix of medium-range interaction distances
                long_range = matrix(long, ntypes, ntypes), # Matrix of long-range interaction distances
                saturation = 4, # Saturation parameter
                fitting_package = "glm", # Fitting package 
                min_dummy = 1e5, # Force at least this many dummy points per species
                max_dummy = 1e5, # Force at maximum this many dummy points per species
                nthreads = 4, # Number of CPU cores to use
                dummy_distribution = "binomial") # Distribution of dummy points
  print(Sys.time() - tm)
  summary_fit <- summary(fit)
  
  if(show) {
    print(fit$coefficients$alpha[[1]])
    print(fit$coefficients$gamma)
    print(summary_fit)
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
                  only_significant = TRUE,
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
                    only_significant = TRUE,
                    base_size = 35))
    dev.off()
  }
}

plot_for_region(list(name = "NVic",
                     plots = c("ANU101", "Ada Tree", "ANU363", "ANU589", "HardyCreek"),
                     abundance = 50))

plot_for_region(list(name = "SVic",
                     plots = c("Weeaproinah", "Turtons", "Lardner"),
                     abundance = 100))

plot_for_region(list(name = "SENSW",
                     plots = c("Newline", "WaratahMix", "WogWay", "Goodenia", "Candelo"),
                     abundance = 100))

plot_for_region(list(name = "NNSW",
                     plots = c("MinesRd", "Tinebank", "Lorne", "BirdTree", "A-Tree", "BlackBull", "Bruxner", "Osullivans"),
                     abundance = 200))

plot_for_region(list(name = "QLD",
                     plots = c("Herberton", "Lamb Range", "Baldy", "Koombooloomba"),
                     abundance = 50))

plot_for_region(list(name = "SWA",
                     plots = c("Frankland", "Clare", "Giants", "Dawson"),
                     abundance = 100))

plot_for_region(list(name = "NWA",
                     plots = c("Collins", "Sutton", "Warren", "Carey", "Dombakup"),
                     abundance = 100))

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
