# KL Divergence Sample Size Optimization for cLHS
# Based on Malone et al. (2019)

# Required libraries
library(clhs)        # For conditioned Latin hypercube sampling
library(terra)       # For raster operations
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(minpack.lm)  # For non-linear fitting

# Silence R CMD check notes for non-standard evaluation
utils::globalVariables(c(
  "sample_size", "kl_divergence", "mean_kl", "sd_kl", "fitted_kl", "cdf"
))

# ============================================================================
# RASTER LOADING UTILITIES
# ============================================================================

#' Load all predictor rasters from a directory into a SpatRaster
#' 
#' @param dir character, directory containing .tif predictor rasters
#' @param pattern filename pattern (default "*.tif")
#' @param recursive logical, search subdirectories
#' @return SpatRaster stack
load_predictor_rasters <- function(dir = "data", pattern = "*.tif", recursive = FALSE) {
  if (!dir.exists(dir)) {
    stop("Directory does not exist: ", dir)
  }
  
  files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (length(files) == 0) {
    stop("No raster files found in ", dir, " matching pattern ", pattern)
  }
  
  rast(files)
}

#' Convenience function: load either a multilayer .tif or a folder of single-layer rasters
#' 
#' @param path path to a .tif file OR a directory containing .tif files
#' @return SpatRaster
load_predictors <- function(path = "data/predictors.tif") {
  if (file.exists(path) && grepl("\\.tif$", path, ignore.case = TRUE)) {
    return(rast(path))
  } else if (dir.exists(path)) {
    return(load_predictor_rasters(path))
  } else {
    stop("Path is neither an existing .tif file nor a directory: ", path)
  }
}

# ============================================================================
# KL DIVERGENCE FUNCTIONS
# ============================================================================

#' Calculate KL divergence between sample and population distributions
#' 
#' @param population_data data frame of population ancillary data
#' @param sample_data data frame of sample ancillary data
#' @param n_bins number of histogram bins (default 25)
#' @return mean KL divergence across all variables
calculate_kl_divergence <- function(population_data, sample_data, n_bins = 25) {
  kl_values <- numeric()
  
  # Calculate KL divergence for each continuous variable
  for (col in names(population_data)) {
    if (is.numeric(population_data[[col]])) {
      
      # Create histogram bins based on population data range
      pop_range <- range(population_data[[col]], na.rm = TRUE)
      breaks <- seq(pop_range[1], pop_range[2], length.out = n_bins + 1)
      
      # Calculate population distribution (Ei)
      pop_hist <- hist(population_data[[col]], breaks = breaks, plot = FALSE)
      pop_density <- pop_hist$counts / sum(pop_hist$counts)
      
      # Calculate sample distribution (Oi)
      sample_hist <- hist(sample_data[[col]], breaks = breaks, plot = FALSE)
      sample_density <- sample_hist$counts / sum(sample_hist$counts)
      
      # Avoid log(0) by adding small constant
      pop_density <- pmax(pop_density, 1e-10)
      sample_density <- pmax(sample_density, 1e-10)
      
      # Calculate KL divergence: KL = Σ Oi * log(Oi) - log(Ei)
      kl <- sum(sample_density * (log(sample_density) - log(pop_density)))
      kl_values <- c(kl_values, kl)
    }
  }
  
  # Return mean KL divergence across all variables
  return(mean(kl_values, na.rm = TRUE))
}

#' Extract ancillary data at sample points
#' 
#' @param population_rasters SpatRaster stack of ancillary variables
#' @param sample_points spatial points (sf object or matrix with x,y coordinates)
#' @return data frame of extracted values
extract_sample_data <- function(population_rasters, sample_points) {
  # Extract values at sample points
  if (inherits(sample_points, "sf")) {
    coords <- sf::st_coordinates(sample_points)
  } else {
    coords <- sample_points
  }
  
  extracted_values <- terra::extract(population_rasters, coords, ID = FALSE)
  return(as.data.frame(extracted_values))
}

# ============================================================================
# SAMPLE SIZE OPTIMIZATION
# ============================================================================

#' Optimize sample size using KL divergence approach (Malone et al. 2019)
#' 
#' @param population_data data frame of population ancillary data
#' @param min_samples minimum sample size to test
#' @param max_samples maximum sample size to test
#' @param step_size increment between sample sizes
#' @param n_replicates number of replicates per sample size
#' @param n_bins number of histogram bins for KL calculation
#' @param probability_threshold CDF threshold for optimal size (default 0.95)
#' @return list with results and plots
optimize_sample_size <- function(population_data, 
                               min_samples = 10, 
                               max_samples = 500, 
                               step_size = 10,
                               n_replicates = 10,
                               n_bins = 25,
                               probability_threshold = 0.95) {
  
  # Generate sequence of sample sizes
  sample_sizes <- seq(min_samples, max_samples, by = step_size)
  
  # Initialize results storage
  results <- data.frame(
    sample_size = integer(),
    replicate = integer(),
    kl_divergence = numeric(),
    stringsAsFactors = FALSE
  )
  
  cat("Running cLHS optimization...\n")
  
  # Main optimization loop
  for (n_samples in sample_sizes) {
    cat(paste("Testing sample size:", n_samples, "\n"))
    
    for (rep in 1:n_replicates) {
      # Run cLHS sampling with error handling
      tryCatch({
        clhs_result <- clhs(population_data, size = n_samples, iter = 10000)
        sample_data <- population_data[clhs_result, ]
        
        # Calculate KL divergence
        kl_div <- calculate_kl_divergence(population_data, sample_data, n_bins)
        
        # Store results
        results <- rbind(results, data.frame(
          sample_size = n_samples,
          replicate = rep,
          kl_divergence = kl_div,
          stringsAsFactors = FALSE
        ))
        
      }, error = function(e) {
        cat(paste("Error at sample size", n_samples, "replicate", rep, ":", e$message, "\n"))
      })
    }
  }
  
  # Check if we have any results
  if (nrow(results) == 0) {
    cat("No successful results obtained\n")
    return(list(
      raw_results = results,
      summary_results = data.frame(),
      fitted_model = NULL,
      fitted_curve = NULL,
      optimal_sample_size = NA,
      plot_kl = NULL,
      plot_cdf = NULL
    ))
  }
  
  # Calculate mean KL divergence for each sample size
  summary_results <- results %>%
    dplyr::group_by(.data$sample_size) %>%
    dplyr::summarise(
      mean_kl = mean(.data$kl_divergence, na.rm = TRUE),
      sd_kl = sd(.data$kl_divergence, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Initialize variables for model fitting
  exp_model <- NULL
  fitted_curve <- NULL
  cdf_values <- NULL
  optimal_sample_size <- NA
  
  # Fit exponential decay function: y = b1 * exp(-k*x) + b0
  cat("Fitting exponential decay function...\n")
  
  if (nrow(summary_results) > 3) {
    tryCatch({
      # Initial parameter estimates
      start_params <- list(
        b0 = min(summary_results$mean_kl), 
        b1 = max(summary_results$mean_kl) - min(summary_results$mean_kl),
        k = 0.01
      )
      
      # Fit exponential model
      exp_model <- minpack.lm::nlsLM(
        mean_kl ~ b1 * exp(-k * sample_size) + b0,
        data = summary_results,
        start = start_params
      )
      
      # Generate fitted curve for plotting
      fitted_curve <- data.frame(
        sample_size = sample_sizes,
        fitted_kl = predict(exp_model, newdata = data.frame(sample_size = sample_sizes))
      )
      
      # Calculate cumulative density function of (1 - KL divergence)
      # This represents the proportion of maximum possible improvement achieved
      max_improvement <- max(fitted_curve$fitted_kl) - min(fitted_curve$fitted_kl)
      if (max_improvement > 1e-10) {
        cdf_values <- (max(fitted_curve$fitted_kl) - fitted_curve$fitted_kl) / max_improvement
        
        # Find optimal sample size where CDF reaches threshold
        optimal_idx <- which(cdf_values >= probability_threshold)[1]
        optimal_sample_size <- ifelse(is.na(optimal_idx), max_samples, sample_sizes[optimal_idx])
      } else {
        optimal_sample_size <- max_samples
      }
      
      cat(paste("Optimal sample size:", optimal_sample_size, "\n"))
      
    }, error = function(e) {
      cat("Error fitting exponential model:", e$message, "\n")
      # Use the largest sample size as fallback
      optimal_sample_size <- max_samples
    })
  } else {
    cat("Not enough data points to fit exponential model\n")
    optimal_sample_size <- max_samples
  }
  
  # Create visualizations
  cat("Creating visualizations...\n")
  
  # Plot 1: KL divergence vs sample size
  p1 <- NULL
  tryCatch({
    p1 <- ggplot(summary_results, aes(x = .data$sample_size, y = .data$mean_kl)) +
      geom_point(size = 2, color = "blue") +
      geom_errorbar(aes(ymin = pmax(.data$mean_kl - .data$sd_kl, 0), 
                        ymax = .data$mean_kl + .data$sd_kl), 
                    width = step_size/2, alpha = 0.7) +
      labs(title = "KL Divergence vs Sample Size",
           x = "Number of Samples",
           y = "KL Divergence") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Add fitted curve if available
    if (!is.null(fitted_curve)) {
      p1 <- p1 + geom_line(data = fitted_curve, 
                           aes(x = .data$sample_size, y = .data$fitted_kl), 
                           color = "red", size = 1, linetype = "solid")
    }
    
    cat("   ✓ KL divergence plot created\n")
  }, error = function(e) {
    cat("   ⚠ Error creating KL plot:", e$message, "\n")
    p1 <- NULL
  })
  
  # Plot 2: Cumulative density function
  p2 <- NULL
  if (!is.null(cdf_values) && !is.na(optimal_sample_size)) {
    tryCatch({
      cdf_data <- data.frame(sample_size = sample_sizes, cdf = cdf_values)
      
      p2 <- ggplot(cdf_data, aes(x = .data$sample_size, y = .data$cdf)) +
        geom_line(size = 1, color = "darkgreen") +
        geom_hline(yintercept = probability_threshold, color = "red", linetype = "dashed") +
        geom_vline(xintercept = optimal_sample_size, color = "red", linetype = "dashed") +
        labs(title = "Cumulative Density Function of (1 - KL Divergence)",
             x = "Number of Samples", 
             y = "CDF of (1 - KL divergence)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)) +
        annotate("text", x = optimal_sample_size + step_size * 2, y = probability_threshold - 0.05,
                 label = paste("Optimal size:", optimal_sample_size), color = "red")
      
      cat("   ✓ CDF plot created\n")
    }, error = function(e) {
      cat("   ⚠ Error creating CDF plot:", e$message, "\n")
      p2 <- NULL
    })
  } else {
    cat("   ⚠ CDF plot not created (insufficient data for model fitting)\n")
  }
  
  # Return results
  return(list(
    raw_results = results,
    summary_results = summary_results,
    fitted_model = exp_model,
    fitted_curve = fitted_curve,
    optimal_sample_size = optimal_sample_size,
    plot_kl = p1,
    plot_cdf = p2
  ))
}

# ============================================================================
# SPATIAL SAMPLING OPTIMIZATION
# ============================================================================

#' Optimize sample size for spatial raster data
#' 
#' @param raster_stack SpatRaster stack of ancillary variables
#' @param min_samples minimum sample size to test
#' @param max_samples maximum sample size to test  
#' @param step_size increment between sample sizes
#' @param n_replicates number of replicates per sample size
#' @return optimization results
optimize_spatial_sampling <- function(raster_stack, 
                                     min_samples = 10,
                                     max_samples = 500,
                                     step_size = 10,
                                     n_replicates = 10) {
  
  # Convert raster to data frame (sample if too large)
  cat("Converting raster data to data frame...\n")
  
  # Get all cell values
  population_values <- terra::values(raster_stack, na.rm = TRUE)
  population_data <- as.data.frame(population_values)
  
  # Remove rows with any NA values
  population_data <- population_data[complete.cases(population_data), ]
  
  cat("Population size after removing NA:", nrow(population_data), "\n")
  
  # If population is very large, take a representative sample
  if (nrow(population_data) > 50000) {
    cat("Large population detected, sampling 50,000 points for analysis...\n")
    sample_idx <- sample(nrow(population_data), 50000)
    population_data <- population_data[sample_idx, ]
  }
  
  # Run optimization
  results <- optimize_sample_size(
    population_data = population_data,
    min_samples = min_samples,
    max_samples = max_samples,
    step_size = step_size,
    n_replicates = n_replicates
  )
  
  return(results)
}

# ============================================================================
# END-TO-END WORKFLOW FUNCTIONS
# ============================================================================

#' End-to-end KL-based sample size optimization using project rasters
#'
#' @param predictor_path path to multi-layer .tif OR directory of single-layer rasters
#' @param output_dir directory to write outputs (CSV + PNG)
#' @param min_samples minimum sample size to test
#' @param max_samples maximum sample size to test
#' @param step_size increment between sample sizes
#' @param n_replicates number of replicates per sample size
#' @param n_bins number of histogram bins for KL calculation
#' @param probability_threshold CDF threshold for optimal size (default 0.95)
#' @param recursive search subdirectories when predictor_path is a directory
#' @return list results as from optimize_sample_size plus file paths
run_kl_optimization <- function(predictor_path = "data",
                                output_dir = "outputs",
                                min_samples = 10,
                                max_samples = 500,
                                step_size = 10,
                                n_replicates = 10,
                                n_bins = 25,
                                probability_threshold = 0.95,
                                recursive = FALSE) {
  
  cat("\n=== KL Sample Size Optimization ===\n")
  cat("Loading predictors from: ", predictor_path, "\n", sep = "")
  
  # Load raster data with error handling
  tryCatch({
    rasters <- load_predictors(predictor_path)
    cat("Predictor layers (", nlyr(rasters), "): ", paste(names(rasters), collapse = ", "), "\n", sep = "")
  }, error = function(e) {
    stop("Error loading raster data: ", e$message)
  })
  
  # Extract values
  population_values <- terra::values(rasters, na.rm = TRUE)
  population_data <- as.data.frame(population_values)
  population_data <- population_data[complete.cases(population_data), ]
  
  cat("Population (complete cases): ", nrow(population_data), " rows\n", sep = "")
  
  if (nrow(population_data) == 0) {
    stop("No valid data found in rasters")
  }
  
  if (nrow(population_data) > 100000) {
    cat("Large dataset detected; sampling 100,000 rows for KL processing to keep runtime reasonable.\n")
    idx <- sample(nrow(population_data), 100000)
    population_data <- population_data[idx, ]
  }
  
  # Run optimization
  results <- optimize_sample_size(
    population_data = population_data,
    min_samples = min_samples,
    max_samples = max_samples,
    step_size = step_size,
    n_replicates = n_replicates,
    n_bins = n_bins,
    probability_threshold = probability_threshold
  )
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write CSV outputs
  raw_path <- file.path(output_dir, "kl_raw_results.csv")
  summary_path <- file.path(output_dir, "kl_summary_results.csv")
  
  write.csv(results$raw_results, raw_path, row.names = FALSE)
  write.csv(results$summary_results, summary_path, row.names = FALSE)
  
  # Write fitted curve if available
  curve_path <- NA
  if (!is.null(results$fitted_curve)) {
    curve_path <- file.path(output_dir, "kl_fitted_curve.csv")
    write.csv(results$fitted_curve, curve_path, row.names = FALSE)
  }
  
  # Save plots with error handling
  kl_plot_path <- NA
  cdf_plot_path <- NA
  
  # Save KL divergence plot
  if (!is.null(results$plot_kl)) {
    kl_plot_path <- file.path(output_dir, "kl_divergence_vs_sample_size.png")
    tryCatch({
      ggplot2::ggsave(kl_plot_path, plot = results$plot_kl, width = 7, height = 5, dpi = 300)
      cat("   ✓ KL divergence plot saved\n")
    }, error = function(e) {
      cat("   ⚠ Warning: Could not save KL plot:", e$message, "\n")
      kl_plot_path <- NA
    })
  } else {
    cat("   ⚠ No KL plot generated\n")
  }
  
  # Save CDF plot
  if (!is.null(results$plot_cdf)) {
    cdf_plot_path <- file.path(output_dir, "kl_cdf_threshold.png")
    tryCatch({
      ggplot2::ggsave(cdf_plot_path, plot = results$plot_cdf, width = 7, height = 5, dpi = 300)
      cat("   ✓ CDF plot saved\n")
    }, error = function(e) {
      cat("   ⚠ Warning: Could not save CDF plot:", e$message, "\n")
      cdf_plot_path <- NA
    })
  } else {
    cat("   ⚠ No CDF plot generated (normal if model fitting failed)\n")
  }
  
  # Report outputs
  cat("\nOutputs written to: ", normalizePath(output_dir, winslash = "/"), "\n", sep = "")
  cat(" - Raw results: ", raw_path, "\n", sep = "")
  cat(" - Summary: ", summary_path, "\n", sep = "")
  if (!is.na(curve_path)) cat(" - Fitted curve: ", curve_path, "\n", sep = "")
  if (!is.na(kl_plot_path)) cat(" - KL plot: ", kl_plot_path, "\n", sep = "")
  if (!is.na(cdf_plot_path)) cat(" - CDF plot: ", cdf_plot_path, "\n", sep = "")
  
  # Add file paths to results
  results$file_paths <- list(
    raw = raw_path,
    summary = summary_path,
    curve = curve_path,
    kl_plot = kl_plot_path,
    cdf_plot = cdf_plot_path
  )
  
  return(results)
}

# ============================================================================
# PLOTTING UTILITIES
# ============================================================================

#' Save optimization plots manually
#' 
#' @param results results from optimize_sample_size()
#' @param output_dir directory to save plots
#' @param prefix file name prefix (default "kl")
#' @return vector of saved file paths
save_optimization_plots <- function(results, output_dir = "outputs", prefix = "kl") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saved_files <- character()
  
  # Save KL divergence plot
  if (!is.null(results$plot_kl)) {
    kl_path <- file.path(output_dir, paste0(prefix, "_divergence_plot.png"))
    tryCatch({
      ggplot2::ggsave(kl_path, plot = results$plot_kl, width = 8, height = 6, dpi = 300)
      cat("✓ KL divergence plot saved to:", kl_path, "\n")
      saved_files <- c(saved_files, kl_path)
    }, error = function(e) {
      cat("✗ Error saving KL plot:", e$message, "\n")
    })
  }
  
  # Save CDF plot
  if (!is.null(results$plot_cdf)) {
    cdf_path <- file.path(output_dir, paste0(prefix, "_cdf_plot.png"))
    tryCatch({
      ggplot2::ggsave(cdf_path, plot = results$plot_cdf, width = 8, height = 6, dpi = 300)
      cat("✓ CDF plot saved to:", cdf_path, "\n")
      saved_files <- c(saved_files, cdf_path)
    }, error = function(e) {
      cat("✗ Error saving CDF plot:", e$message, "\n")
    })
  }
  
  if (length(saved_files) == 0) {
    cat("⚠ No plots were saved\n")
  }
  
  return(saved_files)
}


# Full end-to-end workflow with outputs: uncomment to run
# results <- run_kl_optimization()


