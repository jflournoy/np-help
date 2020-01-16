suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rniftilib))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(parallel))

unpack_stat_from_rds <- function(file, stat_index, subchunk = NULL){
  print(paste0('Unpacking stat ', stat_index, ' from ', file))
  if(is.null(subchunk)){
    rez <- readRDS(file)[stat_index,]
  } else {
    rez <- lapply(readRDS(file)[stat_index,], `[`, subchunk)
  }
  return(rez)
}

save_stat_to_nii <- function(stat_index, all_stats, allmask_vertices, mask, prefix, postfix, statname){
  if(is.null(postfix)){
    outputfilename <- paste0(prefix, '.', statname, '.nii.gz')
  } else {
    outputfilename <- paste0(prefix, '.', statname, '.', postfix, '.nii.gz')
  }
  message('Working on: ', outputfilename)
  
  allmask_vertices_l <- unlist(allmask_vertices)
  
  mask.dims <- dim(mask)
  len <- prod(mask.dims)
  
  the_stats_l <- lapply(all_stats, function(x){unlist(x)})
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
    y[allmask_vertices_l] <- the_stats_l_ul
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
  
  if(!nifti.set.filenames(nim, outputfilename)){
  } else {
    message('error setting filename: ', outputfilename)
  }
  if(is.null(nifti.image.write(nim))){
    message('file written successfully: ', outputfilename)
  } else {
    message('failed to write: ', outputfilename)
  }
  gc_report <- gc(verbose = getOption("verbose"))
  # print(gc_report)
  return(NULL)
}

parser <- ArgumentParser(description='Collect data from a raw dump from neuropointillist (from a very specific model)')
parser$add_argument('--prefix', type="character",
                    help = 'Filename prefix. For example, whatever precedes .0001.nii.gz',
                    required = TRUE)
parser$add_argument('--ncores', type="integer",
                    help = 'Number of cores to use in parallelizing multi-volume images.',
                    required = FALSE, default = 1)
parser$add_argument('--skipchecks', action='store_true',
                    help = 'Skip some file existence checks? (can be useful for testing)')
args <- parser$parse_args()

ncores <- args$ncores

prefix <- args$prefix
mask <- paste0(prefix, '.0001.nii.gz')
designmat_file <- paste0(prefix, '.designmat.rds')
voxeldata_file <- paste0(prefix, '.0001.rds')
raw_results_file <- paste0(prefix, '.0001.rawresults.rds')

if(!args$skipchecks){
  #check if files exist
  message('Checking ', mask)
  stopifnot(file.exists(mask))
  message('Checking ', designmat_file)
  stopifnot(file.exists(designmat_file))
  message('Checking ', voxeldata_file)
  stopifnot(file.exists(voxeldata_file))
  message('Checking ', raw_results_file)
  stopifnot(file.exists(raw_results_file))
  #first load the first resultsfile.
  designmat <- readRDS(designmat_file)
  voxeldata <- readRDS(voxeldata_file)
  stopifnot(dim(voxeldata)[1] == dim(designmat)[1])
}

mask <- nifti.image.read(mask)

all_raw_rds <- sort(dir(path = '.', pattern = '.*rawresults.rds$'))
allmask_files <- unlist(lapply(all_raw_rds, gsub, pattern = 'rawresults.rds$', replacement = 'nii.gz'))

message('Found ', length(all_raw_rds), ' rds files and assumed ', length(allmask_files), ' mask files...')

#check number of stats in the RDS file
a_rez <- readRDS(all_raw_rds[[1]])

n_stats <- length(a_rez[,1])
stat_names <- attributes(a_rez[,1])$names
n_vox <- length(a_rez[3,])
message('Found ', n_stats, ' outputs for ', n_vox, ' voxels in first raw data file: ', all_raw_rds[[1]], '.')
stats_per_voxel <- lapply(1:n_stats, function(i){
  n <- length(a_rez[i,1][[1]])
  message(n, ' values per voxel for ', stat_names[i])
  return(n)
})

message('Reading in all voxel indices... (takes about a minute)')
allmask_vertices <- lapply(allmask_files, function(maskfname){
  mask <- nifti.image.read(maskfname)
  mask.vector <- as.vector(mask[,,])
  mask.vertices <- which(mask.vector >0)
  return(mask.vertices)
})

#MAX_VOLUMES is hardcoded at 100, but this could change in the future
MAX_VOLUMES=100
for(i in 1:length(stat_names)){
system.time({
  n_per_voxel <- stats_per_voxel[[i]]
  chunks <- NULL
  if(n_per_voxel > MAX_VOLUMES){
    chunks <- split(1:n_per_voxel, ceiling(1:n_per_voxel/MAX_VOLUMES))
  }
  message('Unpacking statistic ', i, ': ', stat_names[[i]])
  if(is.null(chunks)){
    allrez <- mclapply(all_raw_rds, unpack_stat_from_rds, stat_index = i, mc.cores = ncores)
    message(sprintf('Size of R data object: %s', format(object.size(allrez), units = 'Mb')))
    nada <- save_stat_to_nii(stat_index = i, 
                             all_stats = allrez, 
                             allmask_vertices = allmask_vertices, 
                             mask = mask, 
                             prefix = prefix, 
                             postfix = NULL,
                             statname = stat_names[i])  
    print(gc())
  } else {
    for(k in 1:length(chunks)){
      achunk <- chunks[[k]]
      allrez <- mclapply(all_raw_rds, unpack_stat_from_rds, stat_index = i, subchunk = achunk, mc.cores = ncores, mc.silent = FALSE)
      message(sprintf('Size of R data object (chunk %d of %d): %s', k, length(chunks), format(object.size(allrez), units = 'Mb')))
      nada <- save_stat_to_nii(stat_index = i, 
                               all_stats = allrez, 
                               allmask_vertices = allmask_vertices, 
                               mask = mask, 
                               prefix = prefix, 
                               postfix = sprintf('k%04d', k),
                               statname = stat_names[i])  
    }
  }
})
}
