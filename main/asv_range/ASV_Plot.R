## Local ASVs ------------------------------------------
library(gamlss)

# ASVs present in each habitat determined using the script ASV_Habitat.R 
# Ranges based on ranges from EM_Range.R
# Identify ASVs that are only present in Waimea, not larger EMP dataset

## Read In
# ASVs per variable (includes ASVs by habitat, trophic, site name)
var_asvs <- readRDS("../data/processed/ASVs_by_variable.rds")

# All ASV ranges
empo_3 <- readRDS("../data/processed/empo_ranges/empo_3_sample_range_correlation.rds")[[1]]

# All ASV occupancies and abundances
occ_abund <-fread("../data/processed/all_ASV_occupancy_abundance.csv")

# metadata
wai_meta <- fread("../data/processed/wai_meta.csv")

## Proportion Localized ASVs
# prop_endemic takes a list of asvs, subsets the range data, and calculates proportion present only in Waimea
# assumes columns "coord_min" and "coord_max"
get_local_asvs <-function(range_table, lat_range = c(21.59, 21.65)){
  new_table <- copy(range_table)
  new_table$local <- F
  # Pull out
  endemic_ASVs  <- new_table[coord_min >= min(lat_range) & coord_max <= max(lat_range), local := T]
  return(new_table)
}

# Identify all local and cosmopolitan ASVs in the data set
# ASVs can be looked up in this table to determine localness
all_local <- get_local_asvs(range_table = empo_3, lat_range = c(21.59, 21.65))


## Local ASVS --------------------------------------------------
# make data tables of ASV ranges/group counts for each habitat, trophic, site

# determine ranges for all variable groups
ranges_list <- lapply( names(var_asvs) , function(var_group) {
  local_by_group <- all_local[OTU_ID %in% var_asvs[[var_group]]]
  local_by_group$group <- var_group
  return(local_by_group)
})

# Convert list to data.table for plotting
ranges_dt <-rbindlist(ranges_list)

# Calculate proportion local vs. total
local_ratios <- lapply(unique(ranges_dt$group), function(var_group) {
                  ranges_dt[group %in% var_group & local, .N] /
                  ranges_dt[group %in% var_group, .N]
})

names(local_ratios) <- unique(ranges_dt$group)

print("Ratios of local ASVS by habitat, site, trophic")
print(local_ratios)

# Separate out by variable
var_list <- list( 
                habitat   = unique(wai_meta$habitat),
                site_name = unique(wai_meta[order(longitude_deg) ,site_name]),
                trophic    = unique(wai_meta$trophic),
                habitroph = unique(wai_meta[order(habitat,trophic), paste(habitat, trophic)])
                )

## Range Boxplots --------------------------------------------------

# get_comparisons takes a vector of factor and generates all possible pairwise comparisons 
# returns result list of integer vectors of length 2 that indicate pair positions

get_comparisons <- function(group_vec) {
  comp_mat <- combn(seq(1, length(unique(group_vec))), m = 2)
  comp_list <-
    lapply(seq_len(ncol(comp_mat)), function(i) comp_mat[, i])
  return(comp_list)
}

# make boxplots for each metadata variable (habitat, site, trophic)
lapply(names(var_list), function(var_name){
  
  # get ASV ranges for a given variable
  var_ranges_dt <- ranges_dt[group %in% var_list[[var_name]]]
  var_ranges_dt[ ,group:= factor(group, levels = var_list[[var_name]])]
  
  # set up comparisons based on number of unique values
  comps <- get_comparisons(var_ranges_dt$group)
  
  p <-
    ggplot(var_ranges_dt, aes(x = group, y = range_coord, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(width = .2) +
    theme_pubr() +
    ylab("ASV Latitudinal Range") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank())+
    labs(title = paste0("ASV ranges by ", var_name))
  
    # just do this part after the fact
    #+ stat_compare_means(comparisons = comps, method = "kruskal.test") # Add pairwise comparisons p-value
  p
  print(paste("saving",var_name, "plot"))
  ggsave(paste0("../outputs/figures/ASV_ranges_by_", var_name, ".png"), plot = p, width = 10, height = 5)
  
  
  #Pull out range and group count mean/sd for each habitat
  range_summary <- var_ranges_dt[, .(
    mean_range = mean(range_coord),
    sd_range = sd(range_coord),
    mean_types = mean(group_count),
    sd_types = sd(group_count)
  ),
  by = group]
  return(list(plot = p, summary = range_summary))
})


## Occupancy (presence in samples) vs. Abundance (total sequences) for ASVs ----------------------
plot(x = log10(occ_abund$occupancy), y = log10(occ_abund$abundance))

## Make a model (gamlss)
model_df <- data.frame(abund = log10(occ_abund$abundance), occ = log10(occ_abund$occupancy))
model_df <- model_df[model_df$abund >0 & model_df$occ > 0,]


gam1 <- gamlss(abund ~ occ, sigma.formula= ~ occ, family = "NO", data = model_df, weights=(occ^2) )
summary(gam1)
gam1rsq <- Rsq(gam1)

# plot the data and the model
plot(model_df$abund ~ model_df$occ, pch = 20, cex = 0.1, xlab="log10(occupancy)", ylab="log10(abundance)")
abline(a= (gam1$mu.coefficients[1]), b = (gam1$mu.coefficients[2]), col = "red", lwd=3)

