#' Classification of bacteria based on their morphological features.
#' 
#' Note: only binary classification (two types of bacteria) is currently supported.
#'
#' @param trainDF A training dataframe.
#' @param testDF A testing dataframe.
#' @param validDF (Optional) A dataframe for external validation.
#' @param predictType At which level prediction is summarized. Can be either "Bacteria" or "Image". 
#' @param seed A random seed.
#' @param max_mem_size A memory limit for H2O Java machine learning engine.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr if_else
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom stringr str_detect
#' @importFrom caret downSample
#' @importFrom pROC roc
#' @importFrom pROC coords
#' @import h2o
#' @export
#' @rdname bacteriaClassification
#' @name bacteriaClassification
bacteriaClassification <- function(trainDF, testDF, validDF=NULL,
                                   predictType=c("Bacteria","Image"),
                                   seed=12345, 
                                   max_mem_size="6G"){
  set.seed(seed)
  lev <- levels(trainDF$"Bacteria")
  if(is.null(lev)){
    trainDF$"Bacteria" <- as.factor(trainDF$"Bacteria")
    lev <- levels(trainDF$"Bacteria")
  }
  
  ## Down-sampling to avoid class imbalance
  trainDF <- caret::downSample(dplyr::select(trainDF, -Bacteria), trainDF$"Bacteria", yname="Bacteria") %>%
    dplyr::select(Bacteria, setdiff(colnames(.), "Bacteria"))
  
  ## Train bacteria classifiers
  h2oServer <- h2o::h2o.init(ip="localhost", port=seed, max_mem_size=max_mem_size, nthreads=6)
  df_train <- h2o::as.h2o(trainDF)
  df_test <- h2o::as.h2o(testDF)
  aml <- h2o::h2o.automl(
    x=setdiff(colnames(trainDF), c("Bacteria", "Source")),  ## remove "Source" metadata to avoid leakage
    y="Bacteria",
    training_frame=df_train,
    validation_frame=df_test,
    max_models=10,
    stopping_metric="mean_per_class_error", stopping_tolerance=0.01, stopping_rounds=1,
    seed=seed
  )
  df_lb <- as.data.frame(aml@leaderboard)
  View(df_lb)
  best_model <- h2o::h2o.getModel(df_lb[["model_id"]][[1]])
  try(h2o::h2o.saveModel(object=best_model, path=getwd(), force=T), silent=T)
  
  ## Evaluate the best classifier
  predEval <- function(evalH2ODF, H2OMod, predictType=c("Bacteria","Image")){
    dir.create(outputHeader, showWarnings=F, recursive=T)
    evalDF <- as.data.frame(evalH2ODF)
    predDF <- as.data.frame(predict(H2OMod, evalH2ODF))
    colnames(predDF)[1] <- "PredictedBacteria"
    predDF <- data.frame("Bacteria"=evalDF$"Bacteria", "Source"=evalDF$"Source", predDF)
    thr <- pROC::coords(pROC::roc(evalDF$"Bacteria", predDF[[lev[1]]]), "best", ret="threshold")
    predDF$"PredictedBacteria" <- dplyr::if_else(predDF[[lev[1]]]>=thr, lev[1], lev[2])
    if(predictType=="Bacteria"){
      return(predDF)
    }else if(predictType=="Image"){
      predDF <- predDF %>% 
        tidyr::gather(ProbBacteria, Probability, -Bacteria, -Source, -PredictedBacteria) %>%
        dplyr::group_by(Bacteria, Source, ProbBacteria) %>%
        dplyr::summarize(Probability=mean(Probability)) %>%
        dplyr::ungroup() %>%
        tidyr::spread(ProbBacteria, Probability)
      return(predDF)
    }
  }
  predDFList <- list()
  predDFList[["LeaderBoard"]] <- df_lb
  predDFList[["Train"]] <- predEval(df_train, best_model, predictType=predictType)
  predDFList[["Test"]] <- predEval(df_test, best_model, predictType=predictType)
  if(!is.null(validDF)){
    df_valid <- h2o::as.h2o(validDF)
    predDFList[["Valid"]] <- predEval(df_valid, best_model, predictType=predictType)
  }
  gc();gc()
  return(predDFList)
}
