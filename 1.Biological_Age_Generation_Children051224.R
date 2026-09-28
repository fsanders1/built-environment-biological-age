# Creation of EpiAge Estimates

# Packages needed for generating EpiAge scores
## devtools::install_github("danbelsky/DunedinPACE")
## devtools::install_github("yiluyucheng/dnaMethyAge")

library(devtools)
library(missForest)
library(DunedinPACE)
library(aries)
library(dplyr)
library(meffil)
library(impute)
library(stringr)
library(matrixStats)
library(dnaMethyAge)
library(mice)
library(methylclock)
library(impute)
library(matrixStats)
library(meffonym)
library(foreign)
library(janitor)
library(arsenal)
library(openxlsx)



############################# Reading in DNAm data ###########################################################################################


setwd(".../DNAm_data_all_w24y/")
aries.dir="woc_release"

aries.samples.mod <- function(path) {
  samplesheet.filename <- file.path(path, "samplesheet2957", "b2957_samplesheet.csv")
  if (!file.exists(samplesheet.filename))
    stop("Invalid path for ARIES data: ", path)
  ret <- read.csv(samplesheet.filename, stringsAsFactors=F)
  rownames(ret) <- ret$Sample_Name
  ret
}

retrieve.gds.dims <- function(gds.filename) {
  stopifnot(file.exists(gds.filename))
  gds.file <- openfn.gds(gds.filename)
  on.exit(closefn.gds(gds.file))
  list(read.gdsn(index.gdsn(gds.file, "row.names")),
       read.gdsn(index.gdsn(gds.file, "col.names")))                  
}

aries.info.mod <- function(path) {
  samples <- aries.samples.mod(path)
  
  gds.filenames <- list.files(file.path(path, "betas"), pattern=".gds$", full.names=T)
  names(gds.filenames) <- sub(".gds$", "", basename(gds.filenames))
  sapply(names(gds.filenames), function(featureset) {
    gds.filename <- gds.filenames[[featureset]]
    ids <- retrieve.gds.dims(gds.filename)
    probe.names <- ids[[1]]
    sample.names <- ids[[2]]        
    samples <- samples[match(sample.names, samples$Sample_Name),]
    
    control.matrix <- read.table(file.path(path, "control_matrix", paste(featureset, "txt", sep=".")),
                                 row.names=1, header=T, sep="\t",check.names=F)
    control.matrix <- as.matrix(control.matrix)
    
    cell.count.filenames <- list.files(file.path(path, "derived", "cellcounts"), ".txt$", full.names=T)
    counts <- sapply(cell.count.filenames, read.table, sep="\t", row.names=1, header=T, simplify=F)
    names(counts) <- sub(".txt", "", sapply(names(counts), basename))
    counts <- lapply(counts, as.matrix)
    counts <- lapply(counts, function(counts) counts[match(sample.names, rownames(counts)),])
    names(counts) <- gsub(" ", "-", names(counts))
    
    list(samples=samples,
         probe.names=probe.names,
         control.matrix=control.matrix,
         cell.counts=counts)
  }, simplify=F)
}

aries.select.mod <- function(path, featureset, time.point=NULL, sample.names=NULL, replicates=F) {
  aries.info <- aries.info.mod(path=path)
  
  aries.select0 <- function(path, featureset, time.point, sample.names, replicates) {
    
    if (!featureset %in% names(aries.info))
      stop("'", featureset, "' is not a valid feature set name. ",
           "Try one of the following: ", paste(names(aries.info), collapse=", "))
    
    aries <- aries.info[[featureset]]
    
    if (!is.null(sample.names))
      replicates <- T
    
    if (is.null(sample.names)) 
      is.selected.sample <- rep(T, nrow(aries$samples))
    else {
      is.selected.sample <- aries$samples$Sample_Name %in% sample.names
      if (!all(sample.names %in% aries$samples$Sample_Name)) {
        if (all(!is.selected.sample))
          stop("None of the sample names appear in this dataset.")
        warning("Some of the sample names do not appear in this dataset.")
      }
    }
    
    if (!is.null(time.point)) {
      if (!all(time.point %in% aries$samples$time_point))
        stop("Time point(s) ", paste(setdiff(time.point, aries$samples$time_point), collapse="/"),
             " not available for feature set '", featureset, "'.")
      is.selected.sample <- is.selected.sample & aries$samples$time_point %in% time.point
    }
    
    if (!replicates)
      is.selected.sample <- is.selected.sample & !aries$samples$duplicate.rm
    
    cell.counts <- lapply(
      aries$cell.counts,
      function(counts) {
        counts[which(is.selected.sample),]
      })
    for (ref in names(cell.counts))
      if (all(is.na(cell.counts[[ref]])))
        cell.counts[[ref]] <- NULL
    
    list(path=path,
         featureset=featureset,
         samples=aries$samples[which(is.selected.sample),],
         probe.names=aries$probe.names,
         control.matrix=aries$control.matrix[which(is.selected.sample),],
         cell.counts=cell.counts)
  }
  
  if (missing(featureset)) {
    featuresets <- setdiff(aries.feature.sets(path), "common")
    nsamples <- sapply(featuresets, function(featureset) {
      tryCatch({
        selection <- aries.select0(
          path,
          featureset,
          time.point,
          sample.names,
          replicates)
        nrow(selection$samples)
      }, error=function(e) {
        0
      })                              
    })        
    if (length(featuresets) > 1 &&
        sort(nsamples,decreasing=T)[2] > 0)
      featureset <- "common"
    else
      featureset <- featuresets[which.max(nsamples)]
  }
  aries.select0(path, featureset, time.point, sample.names, replicates)
}

# Reading in a specific sample of methylation data.
## 1. Need to create samplesheet first with ID's fitting criteria.

samplesheet <- read.csv("...samplesheet.csv")

samplesheet <- subset(samplesheet, time_point=="15up")
samplesheet <- distinct(samplesheet,Sample_Name, .keep_all= TRUE)
samplesheet$IDC <- paste0(samplesheet$cidB2957, samplesheet$QLET) 

# add in siblngs column
data <- readRDS("....rds")
data <- complete(data, action =1)
data <- data[, c("IDC", "mz005l")]
data <- data[data$IDC %in% samplesheet$IDC, ]
samplesheet <- merge(samplesheet, data, by = "IDC", all.x = TRUE)

# Remove duplicates, twins etc
samplesheet_unique <- samplesheet %>% filter(QLET == "A")
message("Twins removed: ", nrow(samplesheet) - nrow(samplesheet_unique))
twins_removed <- setdiff(samplesheet$IDC, samplesheet_unique$IDC)

# remove duplicates
samplesheet_unique2 <- samplesheet_unique %>% filter(duplicate.rm == FALSE)
message("Duplicates removed: ", nrow(samplesheet_unique) - nrow(samplesheet_unique2))
duplicates_removed <- setdiff(samplesheet_unique$IDC, samplesheet_unique2$IDC)

# check sibling variable
message("Sibling status:\n", capture.output(print(table(samplesheet_unique2$mz005l))))

# remove siblings
samplesheet_unique3 <- samplesheet_unique2 %>% filter(mz005l == 2)
message("Siblings and NA removed: ", nrow(samplesheet_unique2) - nrow(samplesheet_unique3))
siblings_removed <- setdiff(samplesheet_unique2$IDC, samplesheet_unique3$IDC)

samplesheet <- samplesheet_unique3 # rename for rest of script (above was to check it worked)
rm(samplesheet_unique2, samplesheet_unique3, samplesheet_unique)

dataset <- readRDS("...rds") 
completed_list <- complete(dataset, "all")

############################## Merge imputed datasets with samplesheet ##############################

# completed_list: list of completed datasets from your mids object
merged_list <- lapply(completed_list, function(d) {
  
  # Only keep rows that exist in samplesheet
  samples <- samplesheet$IDC
  model_data <- d[d$IDC %in% samples, ]
  
  # Merge with Sample_Name from samplesheet
  IDs <- samplesheet[, c("IDC", "Sample_Name")]
  model_data <- merge(model_data, IDs, by = "IDC")
  
  # Keep distinct IDs
  model_data <- distinct(model_data, IDC, .keep_all = TRUE)
  
  return(model_data)
})

############################## Process methylation data ##############################
setwd("...")
horvath_cpgs <- read.csv("Horvath_2018_CpGs.csv", stringsAsFactors = FALSE)$Name
altum_cpgs   <- read.csv("Altum_CpGs.csv", stringsAsFactors = FALSE)$CpG
PACE_cpgs <- DunedinPACE::getRequiredProbes(backgroundList = TRUE)
PACE_cpgs <- unlist(PACE_cpgs)

# add in our smoking cpg site
combined_cpgs <- unique(c(horvath_cpgs, altum_cpgs, PACE_cpgs, "cg05575921"))

setwd("...")

aries <- aries.select.mod(aries.dir, sample.names = samplesheet$Sample_Name)
rm(dataset, completed_list, altum_cpgs, horvath_cpgs, PACE_cpgs)
gc()

aries$meth <- aries.methylation(aries, probe.names = combined_cpgs)
gc()
# Imputing missing cpgs using KMM
# Filter out probes with >30% missing
# ----------------------------------------------------------------------
aries.final <- as.data.frame(aries$meth)
nas.per.cpg <- rowSums(is.na(aries.final))
number.samples <- ncol(aries.final)
is.cpg.more.than.thirtyPercent.is.na <- (nas.per.cpg / number.samples) > 0.30
meth_filtered <- aries.final[!is.cpg.more.than.thirtyPercent.is.na, , drop = FALSE]

message("Filtered out probes with >30% missing values: ", nrow(aries.final) - nrow(meth_filtered))
message("Remaining probes: ", nrow(meth_filtered), ", samples: ", ncol(meth_filtered))
message("Missing values before imputation: ", sum(is.na(meth_filtered)))

# ----------------------------------------------------------------------
# Probe ordering 
# ----------------------------------------------------------------------
ann_450k <- meffil.get.features("450k")
ann_epic <- meffil.get.features("epic")
annotation_data <- rbind(ann_450k, ann_epic)
annotation_data <- annotation_data[!is.na(annotation_data$chr) & !is.na(annotation_data$pos), ]

# Drop duplicates, keeping the first occurrence
annotation_data <- annotation_data[!duplicated(annotation_data$name), ]

# Keep only probes in your data
annotation_data <- annotation_data[annotation_data$name %in% rownames(meth_filtered), ]

# Convert chrX/chrY to numeric for ordering
annotation_data$chr_num <- as.numeric(gsub("chr", "", annotation_data$chr))
annotation_data$chr_num[annotation_data$chr == "chrX"] <- 23
annotation_data$chr_num[annotation_data$chr == "chrY"] <- 24

# Order probes by chr then position
annotation_data <- annotation_data[order(annotation_data$chr_num, annotation_data$pos), ]
ordered_probes <- annotation_data$name
meth_filtered <- meth_filtered[ordered_probes, , drop = FALSE]

# ----------------------------------------------------------------------
# Remove samples with >30% missing
# ----------------------------------------------------------------------
sample_na_prop <- colMeans(is.na(meth_filtered))
high_missing_samples <- names(sample_na_prop[sample_na_prop > 0.30])

if (length(high_missing_samples) > 0) {
  message("Removing ", length(high_missing_samples), " samples with >30% missing values.")
  meth_filtered <- meth_filtered[, !(colnames(meth_filtered) %in% high_missing_samples), drop = FALSE]
}
message("Remaining probes: ", nrow(meth_filtered), ", samples: ", ncol(meth_filtered))

# ----------------------------------------------------------------------
# Chunked KNN imputation
# ----------------------------------------------------------------------
message("Imputing missing values using KNN in chunks...")

chunk_size <- 5000  # adjust based on your memory
n_probes <- nrow(meth_filtered)
chunks <- split(seq_len(n_probes), ceiling(seq_len(n_probes) / chunk_size))

imputed_chunks <- lapply(chunks, function(idx) {
  message("Imputing probes ", min(idx), " to ", max(idx))
  
  # Convert chunk to numeric matrix
  mat <- as.matrix(meth_filtered[idx, , drop = FALSE])
  storage.mode(mat) <- "double"  # ensures numeric type
  
  # Perform KNN imputation
  impute.knn(mat)$data
})

# Recombine chunks into single matrix
meth_imputed <- do.call(rbind, imputed_chunks)

# Restore row and column names
rownames(meth_imputed) <- rownames(meth_filtered)
colnames(meth_imputed) <- colnames(meth_filtered)

message("Missing after imputation: ", sum(is.na(meth_imputed)))
message("Value range after imputation: ", paste(range(meth_imputed), collapse = " - "))
# ----------------------------------------------------------------------
# Winsorize extreme values (>3*IQR per probe)
# ----------------------------------------------------------------------
q1 <- rowQuantiles(meth_imputed, probs = 0.25, na.rm = TRUE)
q3 <- rowQuantiles(meth_imputed, probs = 0.75, na.rm = TRUE)
iqr <- q3 - q1

lower <- q1 - 3 * iqr
upper <- q3 + 3 * iqr

meth_winsor <- meth_imputed
for (i in seq_len(nrow(meth_winsor))) {
  meth_winsor[i, ] <- pmin(pmax(meth_winsor[i, ], lower[i]), upper[i])
}

message("Winsorization complete.")
message("Value range after winsorization: ", paste(range(meth_winsor), collapse = " - "))

aries.final <- meth_winsor

# CpGs are in rows

######### Commented this section out as we have already ran the Python script,
######### but if you need to re-run it, then uncomment out this section ########

# BMIQ normalisation (this only needs to be done for DNAm data used for Altum age)
# So first, create a copy of DNAm data that we only use for horvath 2018 clock
aries.final.horvath <- aries.final

# Define the gold standard
#standard <- meffonym.horvath.standard()

# Ensure your dataframe is a matrix
#meth <- as.matrix(aries.final)

# Apply BMIQ normalization
#dnam_norm <- meffonym.bmiq.calibration(meth, standard)

#aries.final <- dnam_norm
#rm(dnam_norm, ann_450k, ann_epic, annotation_data, aries, chunks, imputed_chunks, meth, meth_imputed, meth_filtered, meth_winsor)
#gc()

# transpose so cpgs in columns for the python script to work
#aries.final <- t(aries.final)
# save this so we can use it to create Altum age in Python script
setwd("...")
#write.csv(aries.final, file = "cpgs_meth.csv", row.names = TRUE)


gc() # free memory

# transpose back
#aries.final <- t(aries.final)

# Run AltumAge python script and open output of that here
altum_ages <- read.csv("altumage_predictions.csv", header = TRUE, stringsAsFactors = FALSE)


############################## Generate EpiAge estimates ##############################

# Horvath 2018 clock
clock_name <- 'HorvathS2018'
Horvath_age <- methyAge(aries.final.horvath, clock = clock_name)
Horvath_age <- dplyr::rename(Horvath_age, Horvath_Age = mAge)
Horvath_age <- dplyr::rename(Horvath_age, Sample_Name = Sample)
altum_ages <- dplyr::rename(altum_ages, Sample_Name = X)

# DunedinPACE
clock_name = 'DunedinPACE'
PACE_age <- methyAge(aries.final.horvath, clock = clock_name)
PACE_age <- dplyr::rename(PACE_age, PACEAge = mAge)
PACE_age <- dplyr::rename(PACE_age, Sample_Name = Sample)

# Add the smokingDNAm column into our dataframe so we can merge
Horvath_age$cg05575921 <- aries.final.horvath["cg05575921", Horvath_age$Sample_Name]

# merge our AltumAge and Horvath 2018 estimates together
ages <- merge(Horvath_age, altum_ages, by = "Sample_Name")
ages <- merge(ages, PACE_age, by = "Sample_Name")

# rename our smoking variable
ages <- dplyr::rename(ages, Smoking_DNAm = cg05575921)

# Create epigenetic age acceleration measure
chronological_age <- samplesheet[, c("age", "Sample_Name")]
ages <- merge(ages, chronological_age, by = "Sample_Name")
ages <- ages %>% dplyr::rename(Chronological_Age = age)

horvath_model <- lm(Horvath_Age ~ Chronological_Age, data = ages)
ages$Horvath_EAA <- residuals(horvath_model)

#PACE_model <- lm(PACEAge ~ Chronological_Age, data = ages)
#ages$PACE_EAA <- residuals(PACE_model)

Altum_model <- lm(AltumAge ~ Chronological_Age, data = ages)
ages$Altum_EAA <- residuals(Altum_model)

############################## Merge biological age data with each imputed dataset ##############################

# 1. Add cell counts to samplesheet
cell_counts <- as.data.frame(aries$cell.counts)
cell_counts <- cell_counts[, (ncol(cell_counts)-5):ncol(cell_counts)]
colnames(cell_counts) <- gsub("^blood\\.gse35069.", "", colnames(cell_counts))
cell_counts$Sample_Name <- rownames(cell_counts)
samplesheet <- merge(samplesheet, cell_counts, by = "Sample_Name")
message("samplesheet rows: ", nrow(samplesheet))

# 2. Merge each dataset with biological ages AND samplesheet metadata
final_list <- lapply(merged_list, function(df) {
  # Merge Horvath/Altum ages
  df <- merge(df, ages, by = "Sample_Name", all.x = TRUE)
  
  # Merge samplesheet metadata (including cell counts)
  # Exclude Sample_Name to avoid duplicate column
  metadata <- samplesheet[, !names(samplesheet) %in% "Sample_Name"]
  df <- merge(df, metadata, by = "IDC", all.x = TRUE)
  
  return(df)
})

############################## Save result ##############################
setwd("...")
saveRDS(final_list, file = paste0("dataset_transport_with_biological_age_", Sys.Date(), ".rds"))



############ Brain Age ######################################################

PATH = "..."
setwd(PATH)
data_old <- read.spss('...sav', to.data.frame = TRUE, use.value.labels = FALSE)

# select data from old project we need (MRI variables and id to merge)
mri_variants <- names(data_old)[grep("\\.(1|2|3)$", names(data_old))]

# Combine with other non-MRI columns you want to select
other_cols <- c('cidB2957','qlet','kz021','ageyears.1','ageyears.2','ageyears.3', 'MRI.study1', 'MRI.study2', 'MRI.study3')

neuro_data <- data_old[, c(other_cols, mri_variants)]

neuro_data$MRI_study <- apply(
  neuro_data[, c("MRI.study1", "MRI.study2", "MRI.study3")],
  1,
  function(x) {
    first <- which(x == 1)[1]
    if (is.na(first)) NA else first
  }
)

# check it worked
sum(!is.na(neuro_data$MRI_study))

# Keep the id columns separate
id_cols <- c('cidB2957','qlet','kz021', 'MRI_study')

# Identify all MRI columns (everything not in id_cols)
mri_cols <- setdiff(names(neuro_data), id_cols)

# Create a function to merge .1/.2/.3 variants by first non-NA
collapse_mri <- function(df, cols) {
  # Get unique base names without .1/.2/.3
  base_names <- unique(gsub("\\.([1-3])$", "", cols))
  
  # Initialize result dataframe with ID columns
  result <- df[, id_cols, drop = FALSE]
  
  # Loop over base names and collapse
  for (base in base_names) {
    variants <- grep(paste0("^", base, "\\.(1|2|3)$"), names(df), value = TRUE)
    # Take the first non-NA value across the variants
    result[[base]] <- apply(df[, variants, drop = FALSE], 1, function(x) {
      x[!is.na(x)][1]  # first non-NA
    })
  }
  
  return(result)
}

# Apply the function
neuro_data <- collapse_mri(neuro_data, mri_cols)

# check it worked for one example feature
sum(!is.na(neuro_data$l_superiorparietal_thickavg))

# Rename sex and ages
neuro_data <- neuro_data %>% 
  dplyr::rename(Sex = kz021)

# Create single age measure (across all studies)
neuro_data$Age <- apply(
       neuro_data[, c("ageyears.1", "ageyears.2", "ageyears.3")],
       1,
       function(x) {
             x[!is.na(x)][1]   # take the first non-NA value
         }
   )

# Check worked (n = 891)
sum(!is.na(neuro_data$Age))

# Drop variables we don't need anymore
neuro_data$MRI.study1 <- NULL
neuro_data$MRI.study2 <- NULL
neuro_data$MRI.study3 <- NULL
neuro_data$ageyears.1 <- NULL
neuro_data$ageyears.2 <- NULL
neuro_data$ageyears.3 <- NULL

# recode sex so 0 and 1, not 1 and 2
neuro_data$Sex[neuro_data$Sex == 1] <- 0
neuro_data$Sex[neuro_data$Sex == 2] <- 1

# Create IDC column
# merge the cidB3067 and qlet to derive an ID column for each G1 offspring
neuro_data$SubjID <-str_c(neuro_data$cidB2957, '', neuro_data$qlet)

# Only keep people who took part in MRI study
neuro_data <- neuro_data %>%
  filter(!is.na(neuro_data$MRI_study))

# Impute missing neuroimaging data
## Specify which variables to temp remove
id_vars <- c("SubjID", "cidB2957", "qlet")

## Remove ID from data so it doesn't get imputed
id_data   <- neuro_data[ , id_vars, drop = FALSE]
to_impute <- neuro_data[ , !(names(neuro_data) %in% id_vars), drop = FALSE]

## Impute
imputed <- missForest(to_impute)

## Add ID back into data
neuro_data_imputed <- cbind(id_data, imputed$ximp)

# load/prepare headers for ENIGMA model input files
setwd("...")
headers <- read.csv("ENIGMA_headers.csv", header = FALSE)

## Cortical thickness                                   
CorticalMeasures_ALSPAC_ThickAvg  <- neuro_data_imputed[, c("SubjID", "l_bankssts_thickavg", "l_caudalanteriorcingulate_thickavg", "l_caudalmiddlefrontal_thickavg",    
                                                           "l_cuneus_thickavg", "l_entorhinal_thickavg", "l_fusiform_thickavg", "l_inferiorparietal_thickavg",       
                                                           "l_inferiortemporal_thickavg", "l_isthmuscingulate_thickavg", "l_lateraloccipital_thickavg", "l_lateralorbitofrontal_thickavg",   
                                                           "l_lingual_thickavg", "l_medialorbitofrontal_thickavg", "l_middletemporal_thickavg", "l_parahippocampal_thickavg",        
                                                           "l_paracentral_thickavg", "l_parsopercularis_thickavg", "l_parsorbitalis_thickavg", "l_parstriangularis_thickavg",       
                                                           "l_pericalcarine_thickavg", "l_postcentral_thickavg", "l_posteriorcingulate_thickavg", "l_precentral_thickavg",              
                                                           "l_precuneus_thickavg", "l_rostralanteriorcingulate_thick", "l_rostralmiddlefrontal_thickavg", "l_superiorfrontal_thickavg",     
                                                           "l_superiorparietal_thickavg", "l_superiortemporal_thickavg", "l_supramarginal_thickavg", "l_frontalpole_thickavg",            
                                                           "l_temporalpole_thickavg", "l_transversetemporal_thickavg", "l_insula_thickavg", "r_bankssts_thickavg",               
                                                           "r_caudalanteriorcingulate_thicka", "r_caudalmiddlefrontal_thickavg", "r_cuneus_thickavg", "r_entorhinal_thickavg",             
                                                           "r_fusiform_thickavg", "r_inferiorparietal_thickavg", "r_inferiortemporal_thickavg", "r_isthmuscingulate_thickavg",     
                                                           "r_lateraloccipital_thickavg", "r_lateralorbitofrontal_thickavg", "r_lingual_thickavg", "r_medialorbitofrontal_thickavg",  
                                                           "r_middletemporal_thickavg", "r_parahippocampal_thickavg", "r_paracentral_thickavg", "r_parsopercularis_thickavg",       
                                                           "r_parsorbitalis_thickavg", "r_parstriangularis_thickavg", "r_pericalcarine_thickavg", "r_postcentral_thickavg",            
                                                           "r_posteriorcingulate_thickavg", "r_precentral_thickavg", "r_precuneus_thickavg", "r_rostralanteriorcingulate_thick", 
                                                           "r_rostralmiddlefrontal_thickavg", "r_superiorfrontal_thickavg", "r_superiorparietal_thickavg", "r_superiortemporal_thickavg",       
                                                           "r_supramarginal_thickavg","r_frontalpole_thickavg", "r_temporalpole_thickavg", "r_transversetemporal_thickavg", 
                                                           "r_insula_thickavg", "lthickness", "rthickness", "lsurfarea",                         
                                                           "rsurfarea", "icv")]              
### Check for duplicated columns
duplicated(CorticalMeasures_ALSPAC_ThickAvg)                                            

## Cortical surface Area
CorticalMeasures_ALSPAC_SurfAvg <- neuro_data[, c("SubjID", "l_bankssts_surfavg", "l_caudalanteriorcingulate_surfav", "l_caudalmiddlefrontal_surfavg",   
                                                         "l_cuneus_surfavg", "l_entorhinal_surfavg", "l_fusiform_surfavg", "l_inferiorparietal_surfavg",      
                                                         "l_inferiortemporal_surfavg", "l_isthmuscingulate_surfavg", "l_lateraloccipital_surfavg", "l_lateralorbitofrontal_surfavg", 
                                                         "l_lingual_surfavg", "l_medialorbitofrontal_surfavg", "l_middletemporal_surfavg", "l_parahippocampal_surfavg",      
                                                         "l_paracentral_surfavg", "l_parsopercularis_surfavg", "l_parsorbitalis_surfavg", "l_parstriangularis_surfavg",      
                                                         "l_pericalcarine_surfavg", "l_postcentral_surfavg", "l_posteriorcingulate_surfavg", "l_precentral_surfavg",             
                                                         "l_precuneus_surfavg", "l_rostralanteriorcingulate_surfa", "l_rostralmiddlefrontal_surfavg", "l_superiorfrontal_surfavg",       
                                                         "l_superiorparietal_surfavg", "l_superiortemporal_surfavg", "l_supramarginal_surfavg", "l_frontalpole_surfavg",           
                                                         "l_temporalpole_surfavg", "l_transversetemporal_surfavg", "l_insula_surfavg", "r_bankssts_surfavg",              
                                                         "r_caudalanteriorcingulate_surfav", "r_caudalmiddlefrontal_surfavg", "r_cuneus_surfavg", "r_entorhinal_surfavg",            
                                                         "r_fusiform_surfavg", "r_inferiorparietal_surfavg", "r_inferiortemporal_surfavg", "r_isthmuscingulate_surfavg",      
                                                         "r_lateraloccipital_surfavg", "r_lateralorbitofrontal_surfavg", "r_lingual_surfavg", "r_medialorbitofrontal_surfavg",  
                                                         "r_middletemporal_surfavg", "r_parahippocampal_surfavg", "r_paracentral_surfavg", "r_parsopercularis_surfavg",       
                                                         "r_parsorbitalis_surfavg", "r_parstriangularis_surfavg", "r_pericalcarine_surfavg", "r_postcentral_surfavg",           
                                                         "r_posteriorcingulate_surfavg", "r_precentral_surfavg", "r_precuneus_surfavg", "r_rostralanteriorcingulate_surfa",
                                                         "r_rostralmiddlefrontal_surfavg", "r_superiorfrontal_surfavg", "r_superiorparietal_surfavg", "r_superiortemporal_surfavg",      
                                                         "r_supramarginal_surfavg", "r_frontalpole_surfavg", "r_temporalpole_surfavg", "r_transversetemporal_surfavg",    
                                                         "r_insula_surfavg", "lthickness", "rthickness", "lsurfarea",
                                                         "rsurfarea", "icv")]   
### Check for duplicated columns
duplicated(CorticalMeasures_ALSPAC_SurfAvg)

## Subcortical Volumes
SubcorticalMeasures_ALSPAC_VolAvg <- neuro_data[, c("SubjID", "llatvent", "rlatvent","lthal", "rthal", 
                                                           "lcaud", "rcaud", "lput", "rput", "lpal",    
                                                           "rpal", "lhippo", "rhippo", "lamyg", "ramyg", 
                                                           "laccumb", "raccumb", "icv")]
### Check for duplicated columns
duplicated(SubcorticalMeasures_ALSPAC_VolAvg)

# Compare with template headers and fix IDP column names as needed

## Cortical Thickness
headers_TC <- headers[2,]
headers_TC <- headers_TC %>% row_to_names(row_number = 1, remove_row = TRUE)
## Cortical Surface area
headers_SA <- headers[3,]
headers_SA <- headers_SA %>% row_to_names(row_number = 1, remove_row = TRUE)
## Cortical Surface area
headers_Vol <- headers[1,]
headers_Vol <- headers_Vol[,c(1:18)]
headers_Vol <- headers_Vol %>% row_to_names(row_number = 1, remove_row = TRUE)
rm(headers)

## Compare and fix Cortical thickness variables
summary(comparedf(headers_TC, CorticalMeasures_ALSPAC_ThickAvg))
### Capitalize l_(for Left) and r_(for Right)
names(CorticalMeasures_ALSPAC_ThickAvg) <- sub('l_', 'L_', names(CorticalMeasures_ALSPAC_ThickAvg))
names(CorticalMeasures_ALSPAC_ThickAvg) <- sub('r_', 'R_', names(CorticalMeasures_ALSPAC_ThickAvg))
### Uncapitalize "aL_" to "al_"
names(CorticalMeasures_ALSPAC_ThickAvg) <- sub('aL_', 'al_', names(CorticalMeasures_ALSPAC_ThickAvg))
# Rename individual cortical thickness variables to match those of the template
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "L_rostralanteriorcingulate_thick"] <- "L_rostralanteriorcingulate_thickavg"
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "R_caudalanteriorcingulate_thicka"] <- "R_caudalanteriorcingulate_thickavg"
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "R_rostralanteriorcingulate_thick"] <- "R_rostralanteriorcingulate_thickavg"
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "icv"] <- "ICV"
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "lthickness"] <- "LThickness"
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "rthickness"] <- "RThickness"
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "lsurfarea"] <- "LSurfArea"
names(CorticalMeasures_ALSPAC_ThickAvg)[names(CorticalMeasures_ALSPAC_ThickAvg) == "rsurfarea"] <- "RSurfArea"
# See revised column names and compare again with template
colnames(CorticalMeasures_ALSPAC_ThickAvg)
summary(comparedf(headers_TC, CorticalMeasures_ALSPAC_ThickAvg))
# variables names are now identical 
rm(headers_TC)

## Compare and fix Cortical surface area variables
summary(comparedf(headers_SA, CorticalMeasures_ALSPAC_SurfAvg))
### Capitalize l_(for Left) and r_(for Right)
names(CorticalMeasures_ALSPAC_SurfAvg) <- sub('l_', 'L_', names(CorticalMeasures_ALSPAC_SurfAvg))
names(CorticalMeasures_ALSPAC_SurfAvg) <- sub('r_', 'R_', names(CorticalMeasures_ALSPAC_SurfAvg))
### Uncapitalize "aL_" to "al_"
names(CorticalMeasures_ALSPAC_SurfAvg) <- sub('aL_', 'al_', names(CorticalMeasures_ALSPAC_SurfAvg))
### Rename individual cortical surface area variables to match those of the template
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "L_caudalanteriorcingulate_surfav"] <- "L_caudalanteriorcingulate_surfavg"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "L_rostralanteriorcingulate_surfa"] <- "L_rostralanteriorcingulate_surfavg"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "R_caudalanteriorcingulate_surfav"] <- "R_caudalanteriorcingulate_surfavg"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "R_rostralanteriorcingulate_surfa"] <- "R_rostralanteriorcingulate_surfavg"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "icv"] <- "ICV"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "lthickness"] <- "LThickness"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "rthickness"] <- "RThickness"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "lsurfarea"] <- "LSurfArea"
names(CorticalMeasures_ALSPAC_SurfAvg)[names(CorticalMeasures_ALSPAC_SurfAvg) == "rsurfarea"] <- "RSurfArea"
# See revised column names and compare again with template
colnames(CorticalMeasures_ALSPAC_SurfAvg)
summary(comparedf(headers_SA, CorticalMeasures_ALSPAC_SurfAvg))
# variables names are now identical
rm(headers_SA)

## Compare and fix subcortical volume variables
summary(comparedf(headers_Vol, SubcorticalMeasures_ALSPAC_VolAvg))
# Rename individuals columns to much those of the template
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "llatvent"] <- "LLatVent"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "rlatvent"] <- "RLatVent"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "lthal"] <- "Lthal"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "rthal"] <- "Rthal"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "lcaud"] <- "Lcaud"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "rcaud"] <- "Rcaud"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "lput"] <- "Lput"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "rput"] <- "Rput"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "lpal"] <- "Lpal"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "rpal"] <- "Rpal"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "lamyg"] <- "Lamyg"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "ramyg"] <- "Ramyg"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "lhippo"] <- "Lhippo"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "rhippo"] <- "Rhippo"
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "laccumb"] <- "Laccumb" 
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "raccumb"] <- "Raccumb" 
names(SubcorticalMeasures_ALSPAC_VolAvg)[names(SubcorticalMeasures_ALSPAC_VolAvg) == "icv"] <- "ICV"
# See column names compare with tamplate,
colnames(SubcorticalMeasures_ALSPAC_VolAvg)
summary(comparedf(headers_Vol, SubcorticalMeasures_ALSPAC_VolAvg))
# variables names are now identical
rm(headers_Vol)

Covariates <- neuro_data[,c("SubjID", "Age", "Sex")]

write.csv(CorticalMeasures_ALSPAC_SurfAvg,'CorticalMeasuresENIGMA_SurfAvg.csv', row.names=FALSE)
write.csv(CorticalMeasures_ALSPAC_ThickAvg,'CorticalMeasuresENIGMA_ThickAvg.csv', row.names=FALSE)
write.csv(SubcorticalMeasures_ALSPAC_VolAvg,'SubcorticalMeasuresENIGMA_VolAvg.csv', row.names=FALSE)
write.csv(Covariates,'Covariates.csv', row.names=FALSE)

load.lib <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      #require( i , character.only = TRUE )
      library(i)
    }
  }
}

#  Then try/install packages...
load.lib( c("ppcor" , "lsmeans" , "multcomp","data.table","plyr","ModelMetrics",
            "caret","gridExtra","Hmisc","pastecs","psych","ggplot2") )

# source functions
cat("Prep: sourcing functions\n")

source("get.means.R")  
source("prepare.files.R")

# derive mean volumnes / thickness and surfarea

cat("Step 2: Obtaining measures of brain age\n")
cat("deriving mean values for thickness, surface and volume\n") 

Thick=get.means("CorticalMeasuresENIGMA_ThickAvg.csv") 
Thick$ICV=NULL  

Surf=get.means("CorticalMeasuresENIGMA_SurfAvg.csv") 
Surf$ICV=NULL 

Vol=get.means("SubcorticalMeasuresENIGMA_VolAvg.csv") 

# merged all together
TS=merge(Thick,Surf,by="row.names")

TSV=merge(TS,Vol,by.x="Row.names",by.y="row.names")

# read in covariates
Covs <- read.csv("Covariates.csv"); #Read in the covariates file

# Check that all of the required columns are present
mcols=c("SubjID","Sex","Age")
colind=match(mcols,names(Covs))
if(length(which(is.na(colind))) > 0){
  stop('At least one of the required columns in your Covariates.csv file is missing. Make sure that the column names are spelled exactly as listed\n
       It is possible that the problem column(s) is: ', mcols[which(is.na(colind))])
}

# Check for duplicated SubjIDs that may cause issues with merging data sets.
if(anyDuplicated(Covs[,c("SubjID")]) != 0) { stop('You have duplicate SubjIDs in your Covariates.csv file.\nMake sure there are no repeat SubjIDs.') }

#combine the files into one dataframe
data = merge(Covs, TSV, by.x="SubjID", by.y="Row.names")

# Only keep those with MRI data
data <- data[complete.cases(data), ]

cat("creating csv files for brainAge estimation\n")
df=prepare.files(data,names(Covs))

males=df$males
females=df$females
rm(df)

invisible(readline(prompt="You should see 2 csv files in your working directory:\n
-females_raw.csv\n
-males_raw.csv.\n

These should be uploaded on the PHOTON-AI platform as described in README file.
Link to ENIGMA model (PHOTON-AI): https://photon-ai.com/enigma_brainage 
Once you have dowloaded the output files, re-name these to 
'males_raw_out.csv' and 'females_raw_out' for males and females and store them in your working directory, respectively"))

# ---- CentileBrain prep ----

# Load data
data_CB <- neuro_data_imputed

# load column headers for CentileBrain model input files
setwd("...")
headers <- read.csv("CentileBrain_headers.csv", header = TRUE)
headers <- headers %>% dplyr::rename(SITE = `o..SITE`)

# Drop variables not required for CentileBrain (list can be adjusted depending on your headers)
data_CB <- data_CB %>% dplyr::select(-c("lthickness", "rthickness", "lsurfarea", "rsurfarea")) # add/remove based on template

# Rename covariates
data_CB <- data_CB %>% 
  dplyr::rename(age = Age,
                sex = Sex,
                SubjectID = SubjID)

# Apply CentileBrain-specific renaming
# compare IDP variables with template headers and fix column names
summary(comparedf(headers, data_CB))
## Capitalize l_(for Left) and r_(for Right)
names(data_CB) <- sub('l_', 'L_', names(data_CB))
names(data_CB) <- sub('r_', 'R_', names(data_CB))
## Uncapitalize "aL_" to "al_"
names(data_CB) <- sub('aL_', 'al_', names(data_CB))
## Rename individual IDPs to match those of the template
names(data_CB)[names(data_CB) == "L_rostralanteriorcingulate_thick"] <- "L_rostralanteriorcingulate_thickavg"
names(data_CB)[names(data_CB)  == "R_caudalanteriorcingulate_thicka"] <- "R_caudalanteriorcingulate_thickavg"
names(data_CB)[names(data_CB)  == "R_rostralanteriorcingulate_thick"] <- "R_rostralanteriorcingulate_thickavg"
names(data_CB)[names(data_CB)  == "L_entorhinal_thickavg"] <- "L_entorhil_thickavg"
names(data_CB)[names(data_CB)  == "L_supramarginal_thickavg"] <- "L_supramargil_thickavg"
names(data_CB)[names(data_CB)  == "R_entorhinal_thickavg"] <- "R_entorhil_thickavg"
names(data_CB)[names(data_CB)  == "R_supramarginal_thickavg"] <- "R_supramargil_thickavg"
names(data_CB)[names(data_CB)  == "L_caudalanteriorcingulate_surfav"] <- "L_caudalanteriorcingulate_surfavg"
names(data_CB)[names(data_CB)  == "L_rostralanteriorcingulate_surfa"] <- "L_rostralanteriorcingulate_surfavg"
names(data_CB)[names(data_CB)  == "R_caudalanteriorcingulate_surfav"] <- "R_caudalanteriorcingulate_surfavg"
names(data_CB)[names(data_CB) == "R_rostralanteriorcingulate_surfa"] <- "R_rostralanteriorcingulate_surfavg"
names(data_CB)[names(data_CB)  == "L_entorhinal_surfavg"] <- "L_entorhil_surfavg"
names(data_CB)[names(data_CB) == "L_supramarginal_surfavg"] <- "L_supramargil_surfavg"
names(data_CB)[names(data_CB) == "R_entorhinal_surfavg"] <- "R_entorhil_surfavg"
names(data_CB)[names(data_CB)  == "R_supramarginal_surfavg"] <- "R_supramargil_surfavg"
names(data_CB)[names(data_CB)  == "lthal"] <- "Lthal"
names(data_CB)[names(data_CB)  == "rthal"] <- "Rthal"
names(data_CB)[names(data_CB)  == "lcaud"] <- "Lcaud"
names(data_CB)[names(data_CB)  == "rcaud"] <- "Rcaud"
names(data_CB)[names(data_CB)  == "lput"] <- "Lput"
names(data_CB)[names(data_CB)  == "rput"] <- "Rput"
names(data_CB)[names(data_CB)  == "lpal"] <- "Lpal"
names(data_CB)[names(data_CB)  == "rpal"] <- "Rpal"
names(data_CB)[names(data_CB) == "lamyg"] <- "Lamyg"
names(data_CB)[names(data_CB)  == "ramyg"] <- "Ramyg"
names(data_CB)[names(data_CB)  == "lhippo"] <- "Lhippo"
names(data_CB)[names(data_CB)  == "rhippo"] <- "Rhippo"
names(data_CB)[names(data_CB)  == "laccumb"] <- "Laccumb" 
names(data_CB)[names(data_CB)  == "raccumb"] <- "Raccumb" 

## compare again with headers
summary(comparedf(headers, data_CB))
# Add required non-imaging metadata
data_CB$SITE <- "ALSPAC_PE"
data_CB$ScannerType <- "3T_General_Electric_HDx"
data_CB$FreeSurfer_Version <- "6.0.0"

# Reorder columns (make sure SITE etc. are in the correct place)
data_CB <- data_CB[, names(headers)]

# Split by sex
CB_male_raw <- data_CB %>% filter(sex == 0)
CB_female_raw <- data_CB %>% filter(sex == 1)

# Export to Excel (CentileBrain needs xlsx not csv)
write.xlsx(CB_male_raw, 'CB_males_raw.xlsx', row.names = FALSE)
write.xlsx(CB_female_raw, 'CB_females_raw.xlsx', row.names = FALSE)

# Input these files into Centile brain age model (https://centilebrain.org/#/brainAge_global)

# Load your outputs here
CB_females_output <- read.csv("...Centile_Outputs/output_file_2025-10-09-12-19-01_MR_predicted_age_female.csv")
CB_females_output_adjusted <- read.csv("...Centile_Outputs/output_file_2025-10-09-12-19-01_Adjusted_MR_predicted_age_female.csv")

CB_males_output <- read.csv("...Centile_Outputs/output_file_2025-10-09-12-22-13_MR_predicted_age_male.csv")
CB_males_output_adjusted <- read.csv("...Centile_Outputs/output_file_2025-10-09-12-22-13_Adjusted_MR_predicted_age_male.csv")

# Rename variables
CB_females_output <- dplyr::rename(CB_females_output, Centile_Pred = x)
CB_females_output_adjusted <- dplyr::rename(CB_females_output_adjusted, Centile_Pred_Adj = x)

CB_males_output <- dplyr::rename(CB_males_output, Centile_Pred = x)
CB_males_output_adjusted <- dplyr::rename(CB_males_output_adjusted, Centile_Pred_Adj = x)

# Merge age predictions into our data
CB_female_raw$Centile_Pred <- CB_females_output$Centile_Pred
CB_female_raw$Centile_Pred_Adj <- CB_females_output_adjusted$Centile_Pred_Adj

CB_male_raw$Centile_Pred <- CB_males_output$Centile_Pred
CB_male_raw$Centile_Pred_Adj <- CB_males_output_adjusted$Centile_Pred_Adj

# Now let's merge females and males together into one dataframe
CB_ages <- rbind(CB_female_raw, CB_male_raw)

# And now let's merge this into our main dataframe with all phenotypes etc.
CB_ages$SubjectID <- str_trim(CB_ages$SubjectID)
completed_list <- lapply(completed_list, function(df) {
  df$IDC <- str_trim(as.character(df$IDC))  
  df
})
dataset_list <- lapply(completed_list, function(x) {
  merge(x, CB_ages, by.x = "IDC", by.y = "SubjectID", all.x = TRUE)
})

# Merge MRI_Study back in
neuro_subset <- neuro_data[, c("SubjID", "MRI_study", "Age")]
neuro_subset <- neuro_data %>% dplyr::rename(Age_at_MRI = Age)
neuro_subset <- neuro_subset[!is.na(neuro_subset$MRI_study), ]
neuro_subset$SubjID <- str_trim(neuro_subset$SubjID)
dataset_list <- lapply(dataset_list, function(df) {
  merge(df, neuro_subset, by.x = "IDC", by.y = "SubjID", all.x = TRUE)
  })

# only keep those with mri data
dataset_list <- lapply(dataset_list, function(df) {
  df[!is.na(df$Centile_Pred), ]
})


# save 
setwd("...RF")
saveRDS(df_list_with_BAA, file = paste0("dataset_transport_with_biological_age_centile", Sys.Date(), ".rds"))

