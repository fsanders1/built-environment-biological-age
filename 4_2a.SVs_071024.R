# SVA

library(meffil) #for SVA
library(dplyr)
library(mice)

# paths
PHENO_path="..."

impute.matrix <- function(x, margin=1, fun=function(x) mean(x, na.rm=T)) {
  if (margin == 2) x <- t(x)
  
  idx <- which(is.na(x) | !is.finite(x), arr.ind=T)
  if (length(idx) > 0) {
    na.idx <- unique(idx[,"row"])
    v <- apply(x[na.idx,,drop=F],1,fun) ## v = summary for each row
    v[which(is.na(v))] <- fun(v)      ## if v[i] is NA, v[i] = fun(v)
    x[idx] <- v[match(idx[,"row"],na.idx)] ##
    stopifnot(all(!is.na(x)))
  }
  
  if (margin == 2) x <- t(x)
  x
}


random.seed <- set.seed(83)


# Read in data
dataset.transport <- readRDS(paste0(PHENO_path, "..."))


############################## SVs ######################################

# Read in DNAm data (we can't read in all DNAm at once as file too large)
library(devtools)
library(aries)

setwd("...") 
aries.dir=".../woc_release"

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

samplesheet <- read.csv("...samplesheet.csv", sep=",")
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

# remove duplicates
samplesheet_unique2 <- samplesheet_unique %>% filter(duplicate.rm == FALSE)
message("Duplicates removed: ", nrow(samplesheet_unique) - nrow(samplesheet_unique2))

# check sibling variable
message("Sibling status:\n", capture.output(print(table(samplesheet_unique2$mz005l))))

# remove siblings
samplesheet_unique3 <- samplesheet_unique2 %>% filter(mz005l == 2)
message("Siblings and NA removed: ", nrow(samplesheet_unique2) - nrow(samplesheet_unique3))

samplesheet <- samplesheet_unique3 # rename for rest of script (above was to check it worked)
rm(samplesheet_unique2, samplesheet_unique3, samplesheet_unique)

samplesheet_450 <- subset(samplesheet, chip=="450k")
samplesheet_epic <- subset(samplesheet, chip=="epic")

dataset <- readRDS("....rds") 
completed_list <- complete(dataset, "all")

############################## Merge imputed datasets with samplesheet ##############################

# Merge for 450k
merged_list_450 <- lapply(dataset.transport, function(d) {
  samples <- samplesheet_450$IDC
  model_data <- d[d$IDC %in% samples, ]
  IDs <- samplesheet_450[, c("IDC", "Sample_Name")]
  model_data <- merge(model_data, IDs, by="IDC")
  model_data <- distinct(model_data, IDC, .keep_all=TRUE)
  return(model_data)
})

# Merge for EPIC
merged_list_epic <- lapply(dataset.transport, function(d) {
  samples <- samplesheet_epic$IDC
  model_data <- d[d$IDC %in% samples, ]
  IDs <- samplesheet_epic[, c("IDC", "Sample_Name")]
  model_data <- merge(model_data, IDs, by="IDC")
  model_data <- distinct(model_data, IDC, .keep_all=TRUE)
  return(model_data)
})

# For 450k
aries_450 <- aries.select.mod(aries.dir, featureset="450", sample.names = samplesheet_450$Sample_Name)
aries_450$meth <- aries.methylation(aries_450)
beta_450 <- aries_450[["meth"]]

# For EPIC
aries_epic <- aries.select.mod(aries.dir, featureset="epic", sample.names = samplesheet_epic$Sample_Name)
aries_epic$meth <- aries.methylation(aries_epic)
beta_epic <- aries_epic[["meth"]]

# For 450k
model_data_450 <- merged_list_450$`1`[merged_list_450$`1`$Sample_Name.x %in% colnames(beta_450), ]
variable_450 <- model_data_450$Horvath_Age
covariates_450 <- model_data_450[, c("Bcell","CD4T","CD8T","Gran","Mono","NK","age",
                                     "contextual_age41mum","child_sex","smoking_child")]

# Similarly for EPIC
model_data_epic <- merged_list_epic$`1`[merged_list_epic$`1`$Sample_Name.x %in% colnames(beta_epic), ]
variable_epic <- model_data_epic$Horvath_Age
covariates_epic <- model_data_epic[, c("Bcell","CD4T","CD8T","Gran","Mono","NK","age",
                                       "contextual_age41mum","child_sex","smoking_child")]


# 450k

featureset="450k"
features <- meffil.get.features(featureset)

stopifnot(length(rownames(beta_450)) > 0 && all(rownames(beta_450) %in% features$name))
stopifnot(ncol(beta_450) == length(variable_450))
stopifnot(is.null(covariates_450) || is.data.frame(covariates_450) && nrow(covariates_450) == ncol(beta_450))

original.variable <- variable_450
original.covariates <- covariates_450


sample.idx <- which(!is.na(variable_450))

cat("Removing", ncol(beta_450) - length(sample.idx), "missing case(s).")

beta_450 <- beta_450[,sample.idx]
variable_450 <- variable_450[sample.idx]

covariates_450 <- covariates_450[sample.idx,,drop=F]

surrogates.ret <- NULL 
beta.sva <- beta_450

autosomal.sites <- meffil.get.autosomal.sites(featureset)
autosomal.sites <- intersect(autosomal.sites, rownames(beta.sva))

most.variable <- length(autosomal.sites)

beta.sva <- beta.sva[autosomal.sites,]

var.idx <- order(rowVars(beta.sva, na.rm=T), decreasing=T)[1:most.variable]
rm(aries, completed_list, beta, beta450, data, dataset, dataset.transport, features, merged_list, original.covariates)
gc()
memory.limit(size = 80000) 

beta.sva <- impute.matrix(beta.sva[var.idx,,drop=F])
gc()

cov.frame <- model.frame(~., data.frame(covariates_450, stringsAsFactors=F), na.action=na.pass)
mod0 <- model.matrix(~., cov.frame)
mod <- cbind(mod0, variable_450)


set.seed(random.seed)
sva.ret <- sva(beta.sva, mod=mod, mod0=mod0, n.sv=10)

#Check SVs aren't associated with Trait (you need to remove those SVs that associate with outcome)
SVs<-sva.ret$sv
SV_check<-apply(SVs,2,function(x) summary(lm(x~variable_450))$coef[2,]) 
SV_check 

row.names(SVs)=samplesheet_450[sample.idx,c("Sample_Name")]

# here, only include those SVs that were not significant in SV_check
samplesheet_450=merge(samplesheet_450,SVs[,c(1,2,3,4,5,6,7,8,9,10)], by.x="Sample_Name", by.y="row.names", all.x=T)

saveRDS(samplesheet_450, file = "samplesheet_450.rds")

# EPIC

featureset="epic"
features <- meffil.get.features(featureset)

stopifnot(length(rownames(beta_epic)) > 0 && all(rownames(beta_epic) %in% features$name))
stopifnot(ncol(beta_epic) == length(variable_epic))
stopifnot(is.null(covariates_epic) || is.data.frame(covariates_epic) && nrow(covariates_epic) == ncol(beta_epic))

original.variable <- variable_epic
original.covariates <- covariates_epic


sample.idx <- which(!is.na(variable_epic))

cat("Removing", ncol(beta_epic) - length(sample.idx), "missing case(s).")

beta_epic <- beta_epic[,sample.idx]
variable_epic <- variable_epic[sample.idx]

covariates_epic <- covariates_epic[sample.idx,,drop=F]

surrogates.ret <- NULL 
beta.sva <- beta_epic

autosomal.sites <- meffil.get.autosomal.sites(featureset)
autosomal.sites <- intersect(autosomal.sites, rownames(beta.sva))

most.variable <- length(autosomal.sites)

beta.sva <- beta.sva[autosomal.sites,]

var.idx <- order(rowVars(beta.sva, na.rm=T), decreasing=T)[1:most.variable]
rm(aries, completed_list, beta, beta_450, beta_epic, data, dataset, dataset.transport, features, merged_list, original.covariates)
gc()
memory.limit(size = 80000) 

beta.sva <- impute.matrix(beta.sva[var.idx,,drop=F])
gc()

cov.frame <- model.frame(~., data.frame(covariates_epic, stringsAsFactors=F), na.action=na.pass)
mod0 <- model.matrix(~., cov.frame)
mod <- cbind(mod0, variable_epic)


set.seed(random.seed)
sva.ret <- sva(beta.sva, mod=mod, mod0=mod0, n.sv=10)

#Check SVs aren't associated with Trait (you would need to remove those SVs that associate with outcome)
SVs<-sva.ret$sv
SV_check<-apply(SVs,2,function(x) summary(lm(x~variable_epic))$coef[2,]) 
SV_check 

row.names(SVs)=samplesheet_epic[sample.idx,c("Sample_Name")]

# here, only include those SVs that were not significant in SV_check
samplesheet_epic=merge(samplesheet_epic,SVs[,c(1,2,3,4,5,6,7,8,9,10)], by.x="Sample_Name", by.y="row.names", all.x=T)

saveRDS(samplesheet_epic, file = "samplesheet_epic.rds")

# Merge them together
samplesheet <- rbind(samplesheet_450, samplesheet_epic)

# Save this final dataset
setwd("...")
saveRDS(samplesheet, file = paste0('SVs_data', Sys.Date(), '.rds'))

