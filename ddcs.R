## DDCS -- Desiderata for Drug Classification Systems
# By Fabricio Kury and Olivier Bodenreider
# Contact: fabricio.kury at nih.gov
# Start: 2016-09-02 02:38
# CCSD v2.R start: 2016-09-09 4:20
# 

library(plyr)
library(XML)
library(RJSONIO)

setwd('C:/Users/kuryfs/Documents/NLM/Projects/Medicare/ddcs/')
resource_dir <- 'Resources/16-09-13/results-20160913/'
output_dir <- 'Output/Kury/16-09-13/'
curtime <- function() format(Sys.time(), "%Y_%m_%d %H_%M")
exec_start_time <- curtime()
message('Execution started at ', exec_start_time, '.')

classes_input <- c('ATC4', 'EPC', 'MoA', 'PE', 'Chem', 'VAClass')

computeMaxDCSDepth <- function(dcs) {
  message('Computing the maximum depth of ', dcs, '.')
  i <<- 0
  countMaxDepth <- function(tree, classes, cur_count = 0) {
    i <<- i+1
    # The variable "classes" is just for debugging. It keeps track of the line of hierarchy up until the current node.
    if(i%%1000 == 0)
      message(i, ' classes processed.')
    if(!exists('rxclassTree', where = tree))
      cur_count
    else
      max(unlist(lapply(tree$rxclassTree, countMaxDepth, classes = c(classes,
        tree$rxclassTree[[1]]$rxclassMinConceptItem[[1]], tree$rxclassTree[[1]]$rxclassMinConceptItem[[2]]),
        cur_count = cur_count+1)))
  }
  
  # A bit of munging around how to write "VA Class" and "ATC"
  if(tolower(dcs) == 'vac' || tolower(dcs) == 'vaclass')
    dcs <- 'va'
  if(tolower(substr(dcs, 1, 3)) == 'atc')
    dcs <- 'ATC1-4'
  
  classes <- xmlToDataFrame(getNodeSet(xmlParse(RCurl::getURL(paste0(
    'https://rxnav.nlm.nih.gov/REST/rxclass/allClasses?classTypes=', dcs))),
    "//rxclassdata//rxclassMinConceptList//rxclassMinConcept//classId"))

  max(unlist(lapply(unlist(classes), function(e) { countMaxDepth(fromJSON(RCurl::getURL(paste0(
      "https://rxnav.nlm.nih.gov/REST/rxclass/classTree.json?classId=", as.character(e)))), as.character(e)) })))
}

computeMinDCSDepth <- function(dcs) {
  message('Computing the minimum depth of ', dcs, '.')
  cur_min_depth <- 1000
  i <<- 0
  countMaxDepth <- function(tree, cur_count = 0) {
    i <<- i+1
    if(i%%1000 == 0)
      message(i, ' classes processed.')
    if(!exists('rxclassTree', where = tree))
      cur_count
    else
      min(unlist(lapply(tree$rxclassTree, countMaxDepth, cur_count = cur_count+1)))
  }
  
  # A bit of munging around how to write "VA Class" and "ATC"
  if(tolower(dcs) == 'vac' || tolower(dcs) == 'vaclass')
    dcs <- 'va'
  if(tolower(substr(dcs, 1, 3)) == 'atc')
    dcs <- 'ATC1-4'
  
  classes <- xmlToDataFrame(getNodeSet(xmlParse(RCurl::getURL(paste0(
    'https://rxnav.nlm.nih.gov/REST/rxclass/allClasses?classTypes=', dcs))),
    "//rxclassdata//rxclassMinConceptList//rxclassMinConcept//classId"))

  min(unlist(lapply(unlist(classes), function(e) { countMaxDepth(fromJSON(RCurl::getURL(paste0(
    "https://rxnav.nlm.nih.gov/REST/rxclass/classTree.json?classId=", e)))) })))
}


desiderataA <- function(classes = classes_input) {
  # Coverage
  readRxNpToClassSummary <- function(class) {
    d <- read.csv(paste0(resource_dir, 'res-RxNprod_to_', class, '-summary-201606.txt'), header=F, sep='|')
    colnames(d) <- c('complete_mapping', 'ts_summary', 'min_drug',
      'map_type', 'nb_class', 'rxn_drug_id', 'rxn_drug_label')
    retval <- data.frame(RxNp_cvrg = sum(d$complete_mapping=='Y')/nrow(d))
    rownames(retval) <- class
    retval
  }
  
  readNDCToClassSummary <- function(class) {
    d <- read.csv(paste0(resource_dir, 'res-NDC_to_', class, '-summary-201606.txt'), header=F, sep='|')
    colnames(d) <- c('ndc_id', 'pyear', 'pmonth', 'total_freq', 'complete_mapping', 'ts_summary', 'min_drug',
      'map_type', 'nb_class')
    retval <- data.frame(ndc_cvrg = length(unique(d$ndc_id[d$complete_mapping=='Y']))/length(unique(d$ndc_id)),
      claims_cvrg = sum(d$total_freq[d$complete_mapping=='Y'])/sum(d$total_freq))
    rownames(retval) <- class
    retval
  }
  
  retval1 <- cbind(do.call(rbind, lapply(classes, readRxNpToClassSummary)),
    do.call(rbind, lapply(classes, readNDCToClassSummary)))
  
  #### THE PART BELOW is completely independent of the part above
  master_ndc <- read.csv('Resources/16-05-18/16-05-18 MASTER_NDC_INFO.csv')
  readNDCtoClass <- function(class) {
    ndc_to_class <- read.csv(paste0(resource_dir, 'res-NDC_to_', class, '-201606.txt'), header=F, sep='|')
    colnames(ndc_to_class) <- c('ndc_id', 'pyear', 'pmonth', 'atc_id', 'atc_name', 'total_freq', 'complete_mapping',
      'min_drug', 'topical_or_systemic', 'clinically_significant', 'status', 'rxn_drug_id', 'rxn_drug_label',
      'rxn_ing_id', 'rxn_ing_label')
    ndc_to_class <- ndc_to_class[ndc_to_class$complete_mapping == 'Y',]
    ndc_to_class
  }
  
  coverageofClass <- function(class)
    length(unique(readNDCtoClass(class)$ndc_id))/length(unique(master_ndc$NDC))
  
  ndc_coverage <- lapply(classes_input, coverageofClass)
  retval2 <- ndc_coverage
}


desiderataB <- function(classes = classes_input) {
  # Granularity
  readClxProd <- function(class) {
    d <- read.csv(paste0(resource_dir, 'res-ClxProd-', class, '-201606.txt'), header=F, sep='|')
    colnames(d) <- c('class_id', 'class_name', 'sum_freq_ok')
    retval <- data.frame(classes = nrow(d), non_empty = sum(d$sum_freq_ok>0),
      non_empty_p = sum(d$sum_freq_ok>0)/nrow(d))
    rownames(retval) <- class
    retval
  }
  do.call(rbind, lapply(classes, readClxProd))
}


desiderataC <- function(classes = classes_input, outfile = 'ddcs - C') {
  # Unambiguity
  
  dcs_ingredients <- list()
  readRxNprodToClass <- function(class) {
    M <- read.csv(paste0(resource_dir, 'res-RxNprod_to_', class, '-201606.txt'), header=F, sep='|')
    colnames(M) <- c('rxn_drug_id', 'rxn_drug_label', 'class_id', 'class_name', 'complete_mapping',
      'min_drug', 'topical_or_systemic', 'clinically_significant', 'step', 'rxn_ing_id', 'rxn_ing_label')
    M <- M[M$min_drug==0 & M$clinically_significant==1,]
    assign('dcs_ingredients', c(dcs_ingredients, length(unique(M$rxn_drug_id))), envir = parent.env(environment()))
    ddply(M, ~ rxn_drug_id, function(e) length(unique(e$class_id)))$V1
  }
  
  message('Processing map files...')
  out_list <- lapply(classes, readRxNprodToClass)
  message('Complete.')
  names(out_list) <- classes
  
  groupByNumberOfClasses <- function(i, number_of_classes)
    return(round(sum(out_list[[i]]==number_of_classes)/dcs_ingredients[[i]], 4))

  max_n_classes <- max(unlist(out_list))
  for(i in 1:length(out_list)) {
    out_list[[i]] <- lapply((max_n_classes+1):1, groupByNumberOfClasses, i = i)
    names(out_list[[i]]) <- c(paste0(max_n_classes+1, ' or more'), paste0(max_n_classes:2, ' classes'), '1 class')
  }
  names(out_list) <- classes

  retval <- do.call(cbind, out_list)
  browser()
}

#
## Computation starts here.
desiderataA()
desiderataB()
desiderataC()
max_class_depths <- data.frame(DCS=classes_input, `Max depth`=unlist(lapply(classes_input, computeMaxDCSDepth)))
min_class_depths <- data.frame(DCS=classes_input, `Min depth`=unlist(lapply(classes_input, computeMinDCSDepth)))
message('Script execution completed at ', curtime(), '.')
