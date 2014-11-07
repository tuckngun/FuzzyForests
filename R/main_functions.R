## WGCNA
# Sample data
load(file = "test.data.small.RData")
load(file = "test.data.large.RData")
load(file = "test.data.vlarge.RData")

# To determine the power threshold for WCGNA:
# Data has to be in a data frame.
# Each row is a sample and each column is a feature. 
# Row names have to be the sample names.
# Need to silence output from pickSoftThreshold

### This is the master function and will automatically choose between the two sub functions based number of features
library(WGCNA)
options(stringsAsFactors = FALSE)

power.threshold = function(data){
  allowWGCNAThreads()
  enableWGCNAThreads()
  powers = c(c(1:10)
             , seq(from = 12
                   , to = 20
                   , by = 2
             ))
  assign("powers", powers, envir = .GlobalEnv)
  num.features = ncol(data)
  if (num.features <= 20000){
    power.threshold.small(data)
    print("Thresholding on full data set")
    return(sft.results)    
  }
  if (num.features > 20000){
    power.threshold.large(data)
    print("Thresholding on a subset of data")
    return(sft.results)
  }
}

#### Sub function
# Power threshold function for small data sets
power.threshold.small = function(data){
  sft = pickSoftThreshold(data
                          , powerVector = powers
  )
  sft.results = cbind(powers
                      , sft$fitIndices$SFT.R.sq
  )
  assign("sft.results", sft.results, envir = .GlobalEnv)
}

#### Sub function
# Power threshold function for large data sets
# Best if the num.chunks variable could adapt to num.features in a smart way. 
# 20000 features at a time is the most I can get into an 8GB machine without crashing.
power.threshold.large = function(data){
  if (num.features %in% 20001:300000){
    num.chunks = floor(num.features/20000)
    sample.size = floor(num.features/num.chunks)
    if (sample.size > 20000){
      # sample.size = 20000
      sample.size = 5000
    }
  }
  else{
    num.chunks = 15
    sample.size = 20000
  }
  for (i in 1:num.chunks){
    sub = sample(data, sample.size)
    # num.features = ncol(test.data.large)
    # sub = sample(test.data.large, sample.size)
    sft = pickSoftThreshold(sub
                            , powerVector = powers
    )
    if (i == 1) {
      sft.results = cbind(powers, sft$fitIndices$SFT.R.sq)
    }
    else {
      sft.results = cbind(sft.results, sft$fitIndices$SFT.R.sq)
    }
  }
  mean.sft = rep(NA, nrow(sft.results))
    for (i in 1:nrow(sft.results)){
    mean.sft[i] = mean(sft.results[i, -1])
  }
  sft.results = cbind(powers
                      , mean.sft
  )
  assign("sft.results", sft.results, envir = .GlobalEnv)
}
---------------------------------------------------------------------
# Need to automate picking of power threshold.
# Plot SFT results to pick best power value.
cex1 = 0.9
plot(thresholding.large[, 1]
     , thresholding.large[, 2]
     , xlab="Soft Threshold (power)"
     , ylab="Scale Free Topology Model Fit,signed R^2"
     , type="n"
     , main = paste("Scale independence")
)
text(thresholding.large[, 1]
     , thresholding.large[, 2]
     , labels=powers
     , cex=cex1
     , col="red"
)
abline(h=0.90
       , col="red"
)
-----------------------------------------------------------------------
# Functions to make WCGNA modules
library(WGCNA)
options(stringsAsFactors = FALSE)  
  
make.modules = function(type # required "default", "d", "loose", "l"
                        , data # required
                        , power # required
                        , maxBlockSize # required
                        , networkType # required "signed", "unsigned"
                        , minModuleSize = 100
                        , detectCutHeight # required for custom WGCNA
                        , mergeCutHeight # required for custom WGCNA
                        , minKMEtoStay # required for custom WGCNA
                        , minCoreKME # required for custom WGCNA
                        ){
  allowWGCNAThreads()
  enableWGCNAThreads()
  if (type %in% c("default", "d")){
    detectCutHeight = 0.995
    mergeCutHeight = 0.25
    minKMEtoStay = 0.3
    minCoreKME = 0.5
  }
  if (type %in% c("loose", "l")){
    detectCutHeight = 0.999
    mergeCutHeight = 0.15
    minKMEtoStay = 0.15
    minCoreKME = 0.3
  }
  bwnet = blockwiseModules(data
                         , power = power
                         , maxBlockSize = maxBlockSize
                         , TOMType = "unsigned"
                         , minModuleSize = minModuleSize
                         , reassignThreshold = 0                            
                         , numericLabels = TRUE
                         , saveTOMs = FALSE
                         , networkType = networkType
                         , verbose = 3
                         , detectCutHeight = detectCutHeight
                         , mergeCutHeight = mergeCutHeight
                         , minKMEtoStay = minKMEtoStay
                         , minCoreKME = minCoreKME
                         )
  module.info.out = sprintf("%s_%s_modules_info.RData"
                            , deparse(substitute(data))
                            , type
                            )
  save(bwnet
       , file = module.info.out)
  print(sprintf("WGCNA results saved to '%s'", module.info.out))
  module.membership = rbind(data, bwnet$colors)
  t.module.membership = as.data.frame(t(module.membership))
  module.membership = t.module.membership
  module.membership.out = sprintf("%s_%s_module_membership.RData"
                                  , deparse(substitute(data))
                                  , type
                                  )
  save(module.membership
       , file = module.membership.out)
  assign("module.membership", module.membership, envir = .GlobalEnv)
  print(sprintf("Module membership saved to '%s'", module.membership.out))
  return(table(bwnet$colors))
}
-------------------------------------------------------------------
## Iterative random forest function and merging of modules to get surviving features
library(randomForest)
library(foreach)
library(doParallel)

# Trait/phenotype should be in a separate object from the data.
# Sample traits
load("factors.RData")

iterative.RF = function(data = module.membership
                        , type # Required "classification" "regression" The code currently only fully works for classification forests
                        , trait # Required
                        , ntree = 50000
                        , num.processors # Required
                        , mtry.factor = 5
                        , drop.fraction = 0.25
                        , stop.fraction = 0.05
                        ){
  print("This step can take upwards of 24 hours if your data set contains large modules (>100,000 members)")
  module.list = names(table(data[, ncol(data)]))
  dir.create(sprintf("top_%sper", stop.fraction*100))
  for (i in 1:length(module.list)){
    module = data[data[, ncol(data)]==module.list[i], ]
    module = as.data.frame(t(module[, -ncol(module)])) # Initial module state
    num.features = ncol(module)
    mtry = sqrt(num.features) * mtry.factor
    sampled.times = ntree * (mtry/num.features)
    if (sampled.times < 10){
      ntree = 10/(mtry/num.features)
      print(sprintf("The number of trees is too low and has been set to %s.", ntree))
    }
    cutoff = ceiling(num.features * stop.fraction)
    while (num.features > cutoff){
      cl = makeCluster(num.processors)
      registerDoParallel(cl)
      # module = module[complete.cases(module), ] # this and the following line are causing a lot of problems. don't think they are really needed
      # module = module[!duplicated(module), ]
      rf = foreach(ntree = rep(ntree/num.processors, num.processors)
                   , .combine = combine
                   , .packages = 'randomForest'
      ) %dopar%
        randomForest(module
                     , trait
                     , type = type
                     , ntree = ntree
                     , mtry = mtry
                     , importance = TRUE
        )
      stopCluster(cl)
      var.importance = importance(rf, type = 1)
      quantile = quantile(var.importance, 0.25)
      trimmed.varlist = as.data.frame(var.importance[var.importance > quantile(var.importance, 0.25), ])
      regions = row.names(trimmed.varlist)
      module = data[regions, ]
      module = as.data.frame(t(module[, -ncol(module)]))
      num.features = ncol(module)
      mtry = sqrt(num.features) * mtry.factor
    }
    outfile = sprintf("top_%sper/module_%s.RData", stop.fraction*100, module.list[i])
    save(module
         , file = outfile)
    progress = sprintf("Done with module %s. Saved as '%s'.", module.list[i], outfile)
    print(progress)
  }
  # Merge modules:
  file.list = list.files(sprintf("top_%sper", stop.fraction*100))
  survivors = NULL
  for (i in 1:length(file.list)) {
    file.to.load = sprintf("top_%sper/%s", stop.fraction*100, file.list[i])
    if (grep("module", file.to.load) == TRUE){
      load(file.to.load)  
      survivors = c(survivors, module)
    }
  }
  survivors = as.data.frame(survivors)
  rownames(survivors) = row.names(module)
  save(survivors
       , file = sprintf("top_%sper_survivors.RData", stop.fraction*100)
  )
  assign("survivors", survivors, envir = .GlobalEnv)
  print(sprintf("Recursive feature selection is complete. Modules have been merged and the surviving features have been saved in 'top_%sper_survivors.RData'", stop.fraction*100))
}
----------------------------------------------
# Iterative RF on survivors
library(randomForest)
library(foreach)
library(doParallel)

survivors.RF = function(data = survivors
                        , type # Required "classification" "regression" Only classification works right now
                        , trait # Required
                        , ntree = 50000
                        , num.processors # Required
                        , mtry.factor = 5
                        , drop.fraction = 0.25
                        , stop.fraction = 0.05
){
  print("This step can take upwards of 12 hours if the number of survivors is over 100,000")
  # Initial module state
  num.features = ncol(data)
  mtry = sqrt(num.features) * mtry.factor
  sampled.times = ntree * (mtry/num.features)
  if (sampled.times < 10){
    ntree = 10/(mtry/num.features)
    print(sprintf("The number of trees is too low and has been set to %s.", ntree))
  }
  cutoff = ceiling(num.features * stop.fraction)
  print(sprintf("This will stop when the number of features <= %s", cutoff))
  while (num.features > cutoff){
    cl = makeCluster(num.processors)
    registerDoParallel(cl)
    rf = foreach(ntree = rep(ntree/num.processors, num.processors)
                 , .combine = combine
                 , .packages = 'randomForest'
                 ) %dopar%
      randomForest(data
                   , trait
                   , type = type
                   , ntree = ntree
                   , mtry = mtry
                   , importance = TRUE
                   )
    stopCluster(cl)
    var.importance = importance(rf, type = 1)
    quantile = quantile(var.importance, 0.25)
    trimmed.varlist = as.data.frame(var.importance[var.importance > quantile(var.importance, 0.25), ])
    regions = row.names(trimmed.varlist)
    data = as.data.frame(t(data))
    data = data[regions, ]    
    data = as.data.frame(t(data))
    num.features = ncol(data)
    print(sprintf("Number of features remaining: %s", num.features))
    mtry = sqrt(num.features) * mtry.factor
  }
  assign("top.survivors", data, envir = .GlobalEnv)
}
----------------------------------------------  
# Looping RF function (to get variable importance and OOB of survivors)
library(randomForest)

post.iterative.importance = function(data = survivors
                                     , type # Required  "classification" "regression" Only classification works right now
                                     , trait # Required
                                     , num.runs = 25
                                     , ntree = 50000
                                     , mtry.factor = 5
                                     , do.trace = ntree/2
                                     ){
  regions = names(data)
  num.features = length(regions)
  mtry = sqrt(num.features) * mtry.factor
  sampled.times = ntree * (mtry/num.features)
  if (sampled.times < 10){
    ntree = 10/(mtry/num.features)
    print(sprintf("The number of trees is too low and has been set to %s.", ntree))
  }
  dummy = rep(1, num.features)
  imp.results = as.data.frame(cbind(regions
                                    , dummy
                                    )
                              , stringsAsFactors = FALSE
                              )
  num.cat = length(levels(trait))
  cat.list = rep(NA, num.cat)
  for (i in 1:num.cat){
   category = paste("cat.", levels(trait)[i], ".error", sep = "")
   assign(category, rep(NA, num.runs))
   cat.list[i] = category
  }
  for (j in 1:num.runs){
    print(sprintf("Starting run #%s.", j))
    rf = randomForest(data
                      , trait
                      , type = type
                      , ntree = ntree
                      , mtry = mtry
                      , importance = TRUE
                      , do.trace = do.trace
                      )
    imp = as.data.frame(importance(rf
                                   , type = 1
                                   )
                        , stringsAsFactors = FALSE
                        )
    imp.results = merge(imp.results
                        , imp
                        , by.x = "regions"
                        , by.y = "row.names"
                        )
    for (k in 1:num.cat){
      current.cat = cat.list[k]
      entries = get(cat.list[k])
      entries[j] = rf$confusion[num.cat^2+k]
      assign(current.cat, entries)      
    }
  }
  # save.image()
  imp.results = imp.results[, -2]
  names(imp.results) = c("region"
                         , seq(1:num.runs)
                         )
  mean.imp = apply(imp.results[, 2:(num.runs+1)], 1, mean)
  imp.results = cbind(imp.results
                      , mean.imp
                      )
  assign("imp.results", imp.results, envir = .GlobalEnv)
  write.table(imp.results[, c(1, (num.runs+2))]
              , file = "varImp_post_iterativeRF.txt"
              , row.names = FALSE
              , sep="\t"
              , quote = FALSE
              )
  OOB.table = rep(NA, num.runs)
  for (l in 1:num.cat){
    entries = get(cat.list[l])
    OOB.table = cbind(OOB.table, entries)
  }
  OOB.table = as.data.frame(OOB.table)
  names(OOB.table) = c(NA, cat.list)
  OOB.table = OOB.table[, -1]
  overall.OOB = apply(OOB.table, 1, mean)
  OOB.table = cbind(OOB.table
                    , overall.OOB)
  OOB.table = rbind(OOB.table
                      , apply(OOB.table, 2, mean)
                      , apply(OOB.table, 2, sd)
  )
  row.names(OOB.table) = c(seq(1:num.runs)
                           , "mean"
                           , "sd"
                           )
  assign("OOB.table", OOB.table, envir = .GlobalEnv) 
  return(OOB.table)
}