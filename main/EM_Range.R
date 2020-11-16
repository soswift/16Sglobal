# Global 16S - Sample Spread Distance Corellation - R Script
# Document Summary:
# This script corellates the number of samples an OTU is found in with it's global geographic range
# Written as a script to run on HPC.

library(data.table)
library(ggplot2)
library(parallel)

# set ggplot theme
theme_set(theme_bw())

# Define Functions ----------

## set_standardID(): Identify the metadata column that contains sample names
# In all other functions, OTUs and metadata will always be linked using column standardID
set_standardID <- function(metadata, ID_column){
  ids <- metadata[[ID_column]]
  metadata[, standardID := ids]
}

## make_pa(): converts a list to presence absence.
# Use lapply to convert a full OTU abundance table to presence absence data

make_pa <- function(Col) {
  Col[Col != 0] <- 1
  return(Col)
}

## gev_value(): conditionally replaces numeric value in list with numeric value from another list

get_value <- function(Val1, Val2){
  
  New_Val <- Val2[Val1 !=0]
  
  return(New_Val)
}

## sub_by_samples(): subsets an OTU table to match a vector of sample names.
# Assumes that OTU table columns are samples, rows are otus.
# If a column of OTU names should be kept, indicate the column name. 

sub_by_samples <- function(samples, otu_table, OTU_name_column = NA){
  
  if(is.na(OTU_name_column)){
    # No OTU name column present
    
    # subset columns by selected sample names
    new_otu  <- otu_table[ , ..samples]
    
    # remove OTUs with rowSums that = 0
    new_otu[, Sums := rowSums(.SD)]
    
  } else {
    # OTU name column present, add it to the vector of sample names
    samples <- c(OTU_name_column, samples)
    
    # subset columns by selected sample names, preserving OTU names
    new_otu  <- otu_table[ , ..samples]
    
    # remove OTUs with rowSums that = 0, exlcuding OTU name column
    new_otu[, Sums := rowSums(.SD), .SDcols = !OTU_name_column]
  }
  
  new_otu <- new_otu[Sums != 0, .SD, .SDcols = !"Sums"]
  return(new_otu)
}

## count_spread(): Based on a grouping column in the metadata, count the total number of groups an OTU is in.

count_spread <-function(otu_table, metadata_table, group_column, key_column = "standardID", OTU_name_column = "OTU_ID"){
  
  # copy data.table
  otus_merge <- copy(otu_table)
  
  # create vector of metadata columns that will be used for aggregating
  meta_cols <- c(group_column, key_column)
  
  # transpose OTU table so samples are rows and OTUs are columns
  otus_merge <- transpose(otus_merge, keep.names = key_column, make.names = OTU_name_column)
  
  # add grouping column to OTU table 
  otus_merge <- merge(otus_merge, metadata_table[, ..meta_cols], by = key_column)
  setkeyv(otus_merge, cols = group_column)
  
  # generate vector of OTU data columns, excluding the metadata columns
  otu_cols <- colnames(otus_merge)
  otu_cols <- otu_cols[!(otu_cols %in% meta_cols)]
  
  # sum OTU counts by group
  otus_merge <- otus_merge[, lapply(.SD, function(x) any(x > 1)), by = group_column, .SDcols = otu_cols]

  # count how many groups each OTU is in
  counts <- colSums(otus_merge[, ..otu_cols])
  
  # return OTU names and counts as a data.table
  result <- data.table(OTU_ID = names(counts), group_count = counts)
  return(result)
  }

## get_geographic_range(): Calculate the geographic range of each OTU based on latitude or longitude.

get_geographic_range <-function(otu_table, metadata_table, coordinate_column, key_column = "standardID", OTU_name_column = "OTU_ID"){

  # copy metadata table prior to editing
  meta_merge <- copy(metadata_table)

  # copy coordinates to new column 'coords', convert to numeric, drop NA values
  suppressWarnings(meta_merge[ , coords := lapply( .SD, as.numeric ), .SDcols = coordinate_column])
  meta_merge <- na.omit(meta_merge, cols = "coords")
  
  # transform coordinate values to absolute value
  meta_merge[ , coords := abs(coords)]
  
  # save a vector of sample names
  samples <- colnames(otu_table)[colnames(otu_table) != OTU_name_column]
  
  # transform OTU counts to logical presence-absence
  otu_table[ , (samples):= lapply(.SD, as.logical), .SDcols = samples]

  # convert to matrix
  otus_mat <- as.matrix(otu_table, rownames = OTU_name_column)
  
  # arrange to match metadata
  otus_mat <- otus_mat[, meta_merge[[key_column]]]
  
  # keep only OTUs that are present in a sample
  keep <- which( apply( otus_mat,1, any ) )
  otus_mat <- otus_mat[ keep,]
  

  # get range of coordinates for each row of the OTU matrix
  coord_ranges <- apply(otus_mat, 1, function(x){
                    range( meta_merge$coords[x], na.rm = T )
  })
  
  # make a nice data.table output with consistent OTU name column
  output <- data.table(coord_min = coord_ranges[1,], coord_max = coord_ranges[2,])
  output[ , (OTU_name_column) := colnames(coord_ranges)]
  
  # calculate total range
  output[ , range_coord := coord_max - coord_min]
  }

# Define Inputs ----------

# Toy data version
otu_file  <- "../../toy_data/toy_data/toy_emp_data.csv"
meta_file <- "../../toy_data/toy_data/toy_emp_metadata.csv"

# full data version

# Read in ----------

# full data version
# all_otus <- fread(otu_file)
# colnames(all_otus)[1] <- "OTU_ID"
# 
# all_meta <- fread(meta_file)
# colnames(emp_meta)[1] <- "SampleID"

all_otus <- fread(otu_file)
all_meta <- fread(meta_file)

# set standard ID for meta
all_meta <- set_standardID(all_meta, "SampleID")

# Subset waimea samples
wai_meta <- all_meta[site_type == "Core Waimea"]

wai_otu <- sub_by_samples(samples = wai_meta$standardID,  otu_table = all_otus, OTU_name_column = "OTU_ID")

# Correlate ----------

wai_otu[1:5,1:5]

# for
for (empo in c("empo_1","empo_2","empo_3")) {
tm <- proc.time()


# for each OTU, count the number of sample types in which it occurs
sample_type_spread <- count_spread(
                                otu_table = wai_otu,
                                metadata_table = wai_meta,
                                group_column = empo
                              )
 proc.time() - tm 
# for each OTU, calculate the geographic range in which it occurs

tm <- proc.time()
geographic_range <- get_geographic_range(otu_table = all_otus,
                                         metadata_table = all_meta,
                                         coordinate_column = "latitude_deg")
proc.time() - tm 

# merge range and spread data
range_dat <- merge(sample_type_spread, geographic_range, by = "OTU_ID")

# save data
saveRDS(list(range_dat = range_dat,
             range_mod = range_mod,
             range_anova = range_anova,
             conf_int = conf_int), 
          paste0(empo,"_sample_range_correlation.rds")
        )
}


