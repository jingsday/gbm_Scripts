# Define the path to the gem folder
gem_folder <- "~/project_GBM/gbm_DATA/gbm_DATA_single_atlas"

# Get the list of files in the gem folder
files <- list.files(gem_folder, pattern = "_batch2", full.names = TRUE)

# Loop through each file and rename it
for (file in files) {
  # Get the filename without the path
  filename <- basename(file)
  
  # Replace '_batch2' with 'v2'
  new_filename <- gsub("_batch2", "v2", filename)
  
  # Construct the full path for the new file
  new_file <- file.path(gem_folder, new_filename)
  
  # Rename the file
  file.rename(file, new_file)
  
  # Print the result
  cat("Renamed:", filename, "->", new_filename, "\n")
