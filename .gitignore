#' 3D Spatiotemporal Clustering for Compound Drought-Heat Events
#' 
#' This package provides functions to identify, track, and analyze compound drought-heat events
#' using 3D spatiotemporal clustering (space-time DBSCAN) on climate data.
#' 
#' @import ncdf4 ggplot2 dbscan factoextra zoo raster sf abind
#' @importFrom stats sd na.omit
NULL

#' Apply 3D DBSCAN clustering to identify compound events
#'
#' @param cdhi_data 3D array [lon, lat, time] of CDHI (Compound Drought-Heat Index)
#' @param eps DBSCAN epsilon parameter (spatial distance threshold)
#' @param minPts Minimum points to form a cluster (default is 4)
#' @return 3D array with cluster labels (same dimensions as input, NA where no event)
#' @export
#'
#' @examples
#' \dontrun{
#' cdhi_data <- ncvar_get(nc, "cdhi")
#' clustered <- cluster_3d_events(cdhi_data, eps = 1.2)
#' }
cluster_3d_events <- function(cdhi_data, eps = 1.2, minPts = 4) {
  
  # Input validation
  if (!is.array(cdhi_data) || length(dim(cdhi_data)) != 3) {
    stop("cdhi_data must be a 3D array [lon, lat, time]")
  }
  
  if (eps <= 0) stop("eps must be positive")
  if (minPts < 1) stop("minPts must be at least 1")
  
  # Get dimensions
  nlon <- dim(cdhi_data)[1]
  nlat <- dim(cdhi_data)[2]
  ntime <- dim(cdhi_data)[3]
  
  # Initialize output array
  clustered <- array(NA, dim = c(nlon, nlat, ntime))
  
  # Process each time step
  for (t in 1:ntime) {
    # Extract current time slice
    slice <- cdhi_data[, , t]
    
    # Replace NA with sentinel value (100) for clustering
    slice_for_cluster <- slice
    na_mask <- is.na(slice_for_cluster)
    slice_for_cluster[na_mask] <- 100
    
    # Create coordinate matrices for DBSCAN
    # X coordinates (longitude index)
    coord_lon <- matrix(NA, nrow = nlon, ncol = nlat)
    coord_lon[, 1:nlat] <- c(1:nlon)
    
    # Y coordinates (latitude index)
    coord_lat <- matrix(rep(c(1:nlat), nlon), 
                        nrow = nlon, ncol = nlat, byrow = TRUE)
    
    # Prepare input matrix: [coord_lon, coord_lat, value]
    input_matrix <- cbind(matrix(coord_lon, ncol = 1, byrow = TRUE),
                          matrix(coord_lat, ncol = 1, byrow = TRUE),
                          matrix(slice_for_cluster, ncol = 1, byrow = TRUE))
    
    # Apply DBSCAN
    db_result <- dbscan::dbscan(input_matrix, eps = eps, minPts = minPts)
    
    # Reshape cluster labels back to 2D
    cluster_2d <- matrix(db_result$cluster, nrow = nlon, ncol = nlat)
    
    # Mark original NA positions (restore NA values)
    cluster_2d[na_mask] <- 100
    cluster_2d[cluster_2d == 100] <- NA
    
    clustered[, , t] <- cluster_2d
  }
  
  return(clustered)
}

#' Calculate event centroid migration trajectories
#'
#' @param event_labels 3D array of event labels (from cluster_3d_events)
#' @return List containing event trajectories and statistics
#' @export
#'
#' @examples
#' \dontrun{
#' trajectories <- calculate_centroid_migration(event_labels)
#' }
calculate_centroid_migration <- function(event_labels) {
  
  if (!is.array(event_labels) || length(dim(event_labels)) != 3) {
    stop("event_labels must be a 3D array [lon, lat, time]")
  }
  
  nlon <- dim(event_labels)[1]
  nlat <- dim(event_labels)[2]
  ntime <- dim(event_labels)[3]
  
  # Initialize storage
  centroids_lon <- matrix(NA, nrow = ntime, ncol = 1)
  centroids_lat <- matrix(NA, nrow = ntime, ncol = 1)
  total_area <- matrix(NA, nrow = ntime, ncol = 1)
  
  # Calculate centroid for each time step
  for (t in 1:ntime) {
    slice <- event_labels[, , t]
    
    # Skip if no events
    if (all(is.na(slice))) next
    
    # Weighted centroid calculation (using event presence = 1)
    event_mask <- !is.na(slice)
    
    if (sum(event_mask) > 0) {
      # Create coordinate grids
      lon_idx <- matrix(1:nlon, nrow = nlon, ncol = nlat)
      lat_idx <- matrix(1:nlat, nrow = nlon, ncol = nlat, byrow = TRUE)
      
      # Calculate weighted centroids
      total_weight <- sum(event_mask)
      centroids_lon[t] <- sum(lon_idx[event_mask]) / total_weight
      centroids_lat[t] <- sum(lat_idx[event_mask]) / total_weight
      total_area[t] <- total_weight
    }
  }
  
  # Convert to geographical coordinates (assuming 0.5 degree grid)
  # Adjust based on your actual grid definition
  lon_geo <- (centroids_lon / 2) - 180
  lat_geo <- (centroids_lat / 2) - 90
  
  # Calculate displacement between consecutive time steps
  displacement_lon <- diff(centroids_lon)
  displacement_lat <- diff(centroids_lat)
  
  # Calculate cumulative displacement
  cum_displacement <- sqrt(cumsum(displacement_lon^2, na.rm = TRUE) + 
                            cumsum(displacement_lat^2, na.rm = TRUE))
  
  return(list(
    centroid_lon = centroids_lon,
    centroid_lat = centroids_lat,
    centroid_lon_geo = lon_geo,
    centroid_lat_geo = lat_geo,
    displacement_lon = displacement_lon,
    displacement_lat = displacement_lat,
    total_area = total_area,
    cumulative_displacement = cum_displacement,
    total_migration_lon = sum(displacement_lon, na.rm = TRUE),
    total_migration_lat = sum(displacement_lat, na.rm = TRUE)
  ))
}

#' Create standard NetCDF output file
#'
#' @param data 3D array [lon, lat, time] to write
#' @param filename Output file path
#' @param var_name Variable name (e.g., "cdhi", "cdti")
#' @param var_units Variable units
#' @param var_description Variable description
#' @param lon_vals Longitude values (default: -179.75 to 179.75 at 0.5 degree)
#' @param lat_vals Latitude values (default: -89.75 to 89.75 at 0.5 degree)
#' @export
#'
#' @examples
#' \dontrun{
#' write_netcdf(data, "output.nc", "cdhi", "1", "CDHI index")
#' }
write_netcdf <- function(data, filename, var_name, var_units, var_description,
                         lon_vals = seq(-179.75, 179.75, 0.5),
                         lat_vals = seq(-89.75, 89.75, 0.5)) {
  
  # Load required package
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required. Please install it.")
  }
  
  # Input validation
  if (!is.array(data) || length(dim(data)) != 3) {
    stop("data must be a 3D array [lon, lat, time]")
  }
  
  nlon <- dim(data)[1]
  nlat <- dim(data)[2]
  ntime <- dim(data)[3]
  
  if (length(lon_vals) != nlon) {
    stop("Length of lon_vals must match first dimension of data")
  }
  if (length(lat_vals) != nlat) {
    stop("Length of lat_vals must match second dimension of data")
  }
  
  # Define dimensions
  lon_dim <- ncdf4::ncdim_def(name = 'longitude', units = 'degrees_east', 
                              vals = lon_vals)
  lat_dim <- ncdf4::ncdim_def(name = 'latitude', units = 'degrees_north', 
                              vals = lat_vals)
  time_dim <- ncdf4::ncdim_def(name = 'time', units = 'months', 
                               vals = 1:ntime)
  
  # Define variable
  var_def <- ncdf4::ncvar_def(name = var_name, units = var_units, 
                              dim = list(lon_dim, lat_dim, time_dim), 
                              missval = NA, prec = 'double')
  
  # Create file and write data
  nc_new <- ncdf4::nc_create(filename = filename, vars = var_def)
  ncdf4::ncvar_put(nc = nc_new, varid = var_def, vals = data)
  ncdf4::ncatt_put(nc = nc_new, varid = 0, attname = 'description', 
                   attval = var_description)
  ncdf4::nc_close(nc_new)
}

#' Process basin-level event trajectories
#'
#' @param clustered_events 3D array of clustered event labels
#' @param basin_shapefile Path to basin shapefile or sf object
#' @param lon_vals Longitude values
#' @param lat_vals Latitude values
#' @return Data frame with basin trajectories
#' @export
#'
#' @examples
#' \dontrun{
#' basin_results <- process_basin_trajectories(event_labels, "basins.shp")
#' }
process_basin_trajectories <- function(clustered_events, basin_shapefile,
                                       lon_vals = seq(-179.75, 179.75, 0.5),
                                       lat_vals = seq(-89.75, 89.75, 0.5)) {
  
  # Load required packages
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Package 'raster' is required")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required")
  }
  
  # Load basins
  if (is.character(basin_shapefile)) {
    basins <- sf::read_sf(basin_shapefile)
  } else {
    basins <- basin_shapefile
  }
  
  # Create raster template
  nlon <- dim(clustered_events)[1]
  nlat <- dim(clustered_events)[2]
  ntime <- dim(clustered_events)[3]
  
  raster_template <- raster::raster(nrows = nlat, ncols = nlon,
                                    xmn = min(lon_vals), xmx = max(lon_vals),
                                    ymn = min(lat_vals), ymx = max(lat_vals))
  
  # Initialize results storage
  n_basins <- nrow(basins)
  results <- data.frame(basin_id = basins$GRIDCODE,
                        total_displacement_lon = NA,
                        total_displacement_lat = NA)
  
  # Process each basin
  for (b in 1:n_basins) {
    basin <- basins[b, ]
    
    # Mask events to basin
    basin_mask <- raster::mask(raster_template, basin)
    
    # Calculate basin statistics
    total_pixels <- sum(!is.na(basin_mask@data@values))
    
    if (total_pixels > 0) {
      # Initialize tracking for this basin
      centroids_lon <- rep(NA, ntime)
      centroids_lat <- rep(NA, ntime)
      
      for (t in 1:ntime) {
        event_slice <- clustered_events[, , t]
        event_slice_masked <- event_slice
        event_slice_masked[is.na(basin_mask@data@values)] <- NA
        
        if (!all(is.na(event_slice_masked))) {
          # Calculate centroid for this time step within basin
          event_mask <- !is.na(event_slice_masked)
          lon_idx <- matrix(1:nlon, nrow = nlon, ncol = nlat)
          lat_idx <- matrix(1:nlat, nrow = nlon, ncol = nlat, byrow = TRUE)
          
          centroids_lon[t] <- sum(lon_idx[event_mask], na.rm = TRUE) / sum(event_mask, na.rm = TRUE)
          centroids_lat[t] <- sum(lat_idx[event_mask], na.rm = TRUE) / sum(event_mask, na.rm = TRUE)
        }
      }
      
      # Calculate total displacement
      results$total_displacement_lon[b] <- sum(diff(centroids_lon), na.rm = TRUE)
      results$total_displacement_lat[b] <- sum(diff(centroids_lat), na.rm = TRUE)
    }
  }
  
  return(results)
}

#' Main pipeline for 3D event detection and tracking
#'
#' @param cdhi_file Path to CDHI (Compound Drought-Heat Index) NetCDF file
#' @param output_dir Directory for output files
#' @param eps DBSCAN epsilon parameter
#' @param time_range Optional time range (indices) to subset
#' @return List of results including clustered events and trajectories
#' @export
#'
#' @examples
#' \dontrun{
#' results <- run_3d_event_detection("cdhi_data.nc", "output/", eps = 1.2)
#' }
run_3d_event_detection <- function(cdhi_file, output_dir = ".", eps = 1.2,
                                   time_range = NULL) {
  
  # Load required package
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required. Please install it.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read CDHI data
  nc_conn <- ncdf4::nc_open(cdhi_file)
  cdhi_data <- ncdf4::ncvar_get(nc_conn, varid = 'cdhi')
  ncdf4::nc_close(nc_conn)
  
  # Subset time if specified
  if (!is.null(time_range)) {
    cdhi_data <- cdhi_data[, , time_range]
  }
  
  # Apply 3D clustering
  cat("Applying 3D DBSCAN clustering to CDHI data...\n")
  clustered_events <- cluster_3d_events(cdhi_data, eps = eps)
  
  # Save clustered events
  output_file <- file.path(output_dir, "clustered_cdhe_events.nc")
  write_netcdf(clustered_events, output_file, "event_label", "1", 
               "3D clustered CDHE event labels")
  cat("Saved clustered events to:", output_file, "\n")
  
  # Calculate centroid trajectories
  cat("Calculating centroid trajectories for CDHEs...\n")
  trajectories <- calculate_centroid_migration(clustered_events)
  
  # Save trajectories as CSV
  traj_file <- file.path(output_dir, "cdhe_centroid_trajectories.csv")
  traj_df <- data.frame(
    time = 1:length(trajectories$centroid_lon),
    centroid_lon = trajectories$centroid_lon,
    centroid_lat = trajectories$centroid_lat,
    centroid_lon_geo = trajectories$centroid_lon_geo,
    centroid_lat_geo = trajectories$centroid_lat_geo,
    area = trajectories$total_area
  )
  write.csv(traj_df, traj_file, row.names = FALSE)
  cat("Saved trajectories to:", traj_file, "\n")
  
  return(list(
    clustered_events = clustered_events,
    trajectories = trajectories,
    output_dir = output_dir
  ))
}

# Example usage (commented out)
# 
# # Example: Run the complete pipeline for CDHI data
# cdhi_file <- "path/to/your/cdhi_data.nc"
# results <- run_3d_event_detection(cdhi_file, output_dir = "results/", eps = 1.2)
# 
# # Example: Read existing CDHI data and apply clustering
# library(ncdf4)
# nc_data <- nc_open("cdhi_data.nc")
# cdhi <- ncvar_get(nc_data, "cdhi")
# nc_close(nc_data)
# clustered <- cluster_3d_events(cdhi, eps = 1.2)
# trajectories <- calculate_centroid_migration(clustered)
