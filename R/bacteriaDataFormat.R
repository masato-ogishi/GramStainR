#' Bacteria data formatting for machine learning
#' 
#' @param bacteriaDetectionResult The result by \code{bacteriaDetection}.
#' @param spl A ratio of training data. Remaining data will be designated as testing data. Disable splitting by providing 0, NA, or NULL.
#' @param seed A random seed.
#' @param withReference Whether the \code{bacteriaDetectionResult} contains bacteria labels for referencing purpose. If set TRUE, letters from original image file names are used as bacteria labels.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom data.table rbindlist
#' @importFrom caret createDataPartition
#' @export
#' @rdname bacteriaDataFormat
#' @name bacteriaDataFormat
bacteriaDataFormat <- function(bacteriaDetectionResult, spl=0.7, seed=12345, withReference=T){
  set.seed(seed)
  if(withReference==F){  ## Splitting must be done with the bacteria labels of provided dataset
    formatBacteriaData.Default <- function(res){
      bactData <- data.table::rbindlist(lapply(res, function(d){d$"PixelDF"})) %>%
        dplyr::select("Source", setdiff(colnames(.), c("Source", "ObjID")))
      return(bactData)
    }
    return(formatBacteriaData.Default(bacteriaDetectionResult))
  }else{
    bactImageLabels <- as.factor(gsub("[[:digit:]]+$", "", names(bacteriaDetectionResult)))
    formatBacteriaData.Default <- function(res){
      bactData <- data.table::rbindlist(lapply(res, function(d){d$"PixelDF"})) %>%
        dplyr::mutate("Bacteria"=factor(gsub("[[:digit:]]+$", "", Source), levels=levels(bactImageLabels))) %>%
        dplyr::select("Bacteria", "Source", setdiff(colnames(.), c("Bacteria", "Source", "ObjID")))
      return(bactData)
    }
    if(any(c(spl==0, is.na(spl), is.null(spl)))){
      return(formatBacteriaData.Default(bacteriaDetectionResult))
    }else{
      trainIDs <- as.numeric(caret::createDataPartition(bactImageLabels, p=spl, list=F))
      testIDs <- setdiff(1:length(bactImageLabels), trainIDs)
      res <- list("TrainDF"=formatBacteriaData.Default(bacteriaDetectionResult[trainIDs]),
                  "TestDF"=formatBacteriaData.Default(bacteriaDetectionResult[testIDs]))
      return(res)
    }
  }
}
