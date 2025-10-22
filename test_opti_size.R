source("kl_optimization.r")

# Load predictor rasters from directory or multi-layer file
rasters <- load_predictors('data/predictors.tif')

# Load predictor rasters from directory
#rasters_layers <- load_predictor_rasters('data', pattern = "*.tif", recursive = FALSE)

population_values <- terra::values(rasters, na.rm = TRUE)

# Convert extracted values to a data frame
population_data <- as.data.frame(population_values)

# Remove any rows with missing values
population_data <- population_data[complete.cases(population_data), ]

# cLHS implementation with error handling
clhs_result <- clhs(population_data, size = 100, iter = 10000)
sample_data <- population_data[clhs_result, ]

# Calculate KL divergence
kl_div <- calculate_kl_divergence(population_data, sample_data, n_bins = 25)
print(kl_div)

# Optimize sample size - Option 1
opt_size <- optimize_sample_size(population_data, step_size = 20, min_samples = 30, max_samples = 300, n_bins = 25)
print(opt_size)

# Optimize spatial sampling - Option 2
opt_spatial <- optimize_spatial_sampling(rasters, step_size = 10, min_samples = 30, max_samples = 200)
print(opt_spatial)

opt_kl <- run_kl_optimization('data', step_size = 10, min_samples = 30, max_samples = 200)
print(opt_kl)
