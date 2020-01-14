save_stat_to_nii <- function(stat_index, all_stats, allmask_vertices, mask, prefix, statname){
  suppressPackageStartupMessages(library(data.table))
  
  mask.dims <- dim(mask)
  len <- prod(mask.dims)
  
  the_stats_l <- lapply(all_stats, function(x){unlist(x[stat_index,])})
  stat_1 <- the_stats_l[[1]]
  mask_vert_1 <- allmask_vertices[[1]]
  
  # make sure that we are filling in data of the right size
  # there are two options - either there is one value per vertex
  # or the output is four dimensional
  is3d <-length(mask_vert_1)==length(stat_1)
  is4d <- length(mask_vert_1) < length(stat_1)
  
  #they need to have the same number of chunks
  stopifnot(length(the_stats_l) == length(allmask_vertices))
  n_chunks <- length(the_stats_l)
  
  stopifnot(is3d || is4d)
  if (is3d) {
    the_stats_l_ul <- unlist(the_stats_l, recursive = F)
    y <- vector(mode="numeric", length=len)
    y[mask.vertices] <- the_stats_l_ul
    nim <- nifti.image.copy.info(mask)
    nifti.image.setdatatype(nim, "NIFTI_TYPE_FLOAT32")
    nim$dim <- mask.dims
    nifti.image.alloc.data(nim)
    nim[,,] <- y
  } else {# a 4d image
    nvolumes <- length(stat_1)/length(mask_vert_1)
    nvolumes.trunc <- trunc(nvolumes)
    if (nvolumes != nvolumes.trunc) {
      cat("processVoxel returned a 4d volume that is not an even multiple of the mask.\n")
      stop("Check your processVoxel function.")
    }
    
    the_stats_l_ul <- unlist(the_stats_l, recursive = F)
    the_stats_dt <- data.table(Z = the_stats_l_ul, 
                               volume = names(the_stats_l_ul), 
                               vertex = rep(allmask_vertices_l, each = nvolumes))
    setorder(the_stats_dt, volume, vertex)
    
    ####YOU NEED TO GET THE CONTIGUOUS 1028 voxels together from each bootstrap (eg. t0001) before rotating through the bootstraps!!
    y <- vector(mode="numeric", length=len)
    y[allmask_vertices_l] <- 1
    y <- rep(y, nvolumes)
    y[y==1] <- the_stats_dt$Z
    
    nim <- nifti.image.copy.info(mask)
    nifti.image.setdatatype(nim, "NIFTI_TYPE_FLOAT32")
    dims <- c(mask.dims, nvolumes)
    nim$dim <- dims
    nifti.image.alloc.data(nim)
    #test_array <- array(y, dim = dims)
    nim[,,,] <- y
  }
  outputfilename <- paste0(prefix, '.', statname, '.nii.gz')
  if(!nifti.set.filenames(nim, outputfilename)){
    message('filename set correctly')
  } else {
    message('error setting filename')
  }
  if(is.null(nifti.image.write(nim))){
    message('file written successfully: ', outputfilename)
  } else {
    message('failed to write: ', outputfilename)
  }
}

setwd("~/NewNeuropoint/eptot.gad_lead.fear/eptot.gad_lead.fear.amy_mask")

suppressPackageStartupMessages(library(Rniftilib))

prefix <- 'eptot.gad_lead.fear'
mask <- paste0(prefix, '.0001.nii.gz')
designmat_file <- paste0(prefix, '.designmat.rds')
voxeldata_file <- paste0(prefix, '.0001.rds')
raw_results_file <- paste0(prefix, '.0001.rawresults.rds')

all_raw_rds <- sort(dir(path = '.', pattern = '.*rawresults.rds'))
allrez <- lapply(all_raw_rds, readRDS)

allmask_files <- sort(dir(path = '.', pattern = '.*[0-9]{4}.nii.gz'))
allmask_vertices <- lapply(allmask_files, function(maskfname){
  mask <- nifti.image.read(maskfname)
  mask.vector <- as.vector(mask[,,])
  mask.vertices <- which(mask.vector >0)
  return(mask.vertices)
})

allmask_vertices_l <- unlist(allmask_vertices)

#first load the first resultsfile.
designmat <- readRDS(designmat_file)
voxeldata <- readRDS(voxeldata_file)
stopifnot(dim(voxeldata)[1] == dim(designmat)[1])

results <- readRDS(raw_results_file)
mask <- nifti.image.read(mask)


stat_names <- attributes(results[,1])$names


