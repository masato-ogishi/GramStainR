#' Classification of bacteria based on their morphological features.
#'
#' @param trainDF A training dataframe.
#' @param testDF A testing dataframe.
#' @param validDF (Optional) A dataframe for external validation.
#' @param seed A random seed.
#' @param outputDir A directory for saving outputs.
#' @param max_mem_size A memory limit for H2O Java machine learning engine.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr if_else
#' @importFrom dplyr bind_cols
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
                                   seed=12345,
                                   outputDir="./Results/",
                                   max_mem_size="6G"){
  set.seed(seed)
  dir.create(outputDir, showWarnings=F, recursive=T)

  ## Define bacteria species to be classified
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
    x=setdiff(colnames(trainDF), c("Bacteria", "Source", "ObjID")),  ## remove object metadata to avoid leakage
    y="Bacteria",
    training_frame=df_train,
    validation_frame=df_test,
    max_models=100,
    stopping_metric="AUTO",
    stopping_tolerance=0.05,
    stopping_rounds=5,
    seed=seed
  )
  df_lb <- as.data.frame(aml@leaderboard)
  View(df_lb)
  best_model <- h2o::h2o.getModel(df_lb[["model_id"]][[1]])
  try(h2o::h2o.saveModel(object=best_model, path=outputDir, force=T), silent=T)

  ## Evaluate the best classifier
  predEval <- function(evalH2ODF, H2OMod){
    evalDF <- as.data.frame(evalH2ODF)
    predDF <- as.data.frame(predict(H2OMod, evalH2ODF))
    colnames(predDF)[1] <- "PredictedBacteria"
    predDF <- dplyr::bind_cols(dplyr::select(evalDF, c("Bacteria", "Source", "ObjID")), predDF)
    predDF_Bacteria <- predDF
    predDF_Image <- predDF %>%
        tidyr::gather(ProbBacteria, Probability, -Bacteria, -Source, -ObjID, -PredictedBacteria) %>%
        dplyr::group_by(Bacteria, Source, ProbBacteria) %>%
        dplyr::summarize(Probability=mean(Probability)) %>%  ## dplyr::summarize(Probability=quantile(Probability, 0.9, names=F)) %>%
        dplyr::ungroup() %>%
        tidyr::spread(ProbBacteria, Probability)
    return(list("Bacteria"=predDF_Bacteria, "Image"=predDF_Image))
  }
  predDFList <- list()
  predDFList[["LeaderBoard"]] <- df_lb
  predDFList[["Train"]] <- predEval(df_train, best_model)
  predDFList[["Test"]] <- predEval(df_test, best_model)
  if(!is.null(validDF)){
    df_valid <- h2o::as.h2o(validDF)
    predDFList[["Valid"]] <- predEval(df_valid, best_model)
  }
  rm(list=setdiff(ls(), "predDFList"))
  gc();gc()
  return(predDFList)
}
