#' Classification of bacteria based on their morphological features.
#' 
#' @param trainDF A training dataframe.
#' @param testDF A testing dataframe.
#' @param validDF (Optional) A dataframe for external validation.
#' @param withReference Whether the \code{bacteriaDetectionResult} contains bacteria labels for referencing purpose. If set TRUE, letters from original image file names are used as bacteria labels. Currently, only a binary classification is implemented.
#' @param predictType At which level the predictive performance is evaluated. Can be either "Bacteria" or "Image".
#' @param seed A random seed.
#' @param max_mem_size A memory limit for H2O Java machine learning engine.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom tidyr gather
#' @importFrom readr write_csv
#' @importFrom stringr str_detect
#' @importFrom caret confusionMatrix 
#' @importFrom classifierplots classifierplots_folder
#' @import h2o
#' @export
#' @rdname bacteriaClassification
#' @name bacteriaClassification
bacteriaClassification <- function(trainDF, testDF, validDF=NULL, 
                                   withReference=T, predictType=c("Bacteria","Image"), 
                                   seed=12345, max_mem_size="30G"){
  ## Train a bacteria image classifier
  h2oServer <- h2o::h2o.init(ip="localhost", port=seed, max_mem_size=max_mem_size, nthreads=6)
  df_train <- h2o::as.h2o(dplyr::select(trainDF, -Source))
  df_test <- h2o::as.h2o(testDF)
  aml <- h2o::h2o.automl(
    y="Bacteria",
    training_frame=df_train,
    validation_frame=df_test,
    max_models=10,
    stopping_metric="mean_per_class_error", stopping_tolerance=0.01, stopping_rounds=1,
    seed=seed
  )
  df_lb <- as.data.frame(aml@leaderboard)
  readr::write_csv(df_lb, "./Results/LeaderBoard.csv")
  View(df_lb)
  best_model <- h2o::h2o.getModel(df_lb[["model_id"]][[1]])
  
  ## Evaluate the best classifier
  predEval <- function(evalH2ODF, H2OMod, withReference=T, predictType=c("Bacteria","Image"), outputHeader=""){
    dir.create(outputHeader, showWarnings=F, recursive=T)
    evalDF <- as.data.frame(evalH2ODF)
    predDF <- as.data.frame(predict(H2OMod, evalH2ODF))
    colnames(predDF)[1] <- "PredictedBacteria"
    lev <- levels(predDF$"PredictedBacteria")
    if(withReference==F){
      predDF <- data.frame("Source"=evalDF$"Source", predDF)
      if(predictType=="Bacteria"){
        return(predDF)
      }else if(predictType=="Image"){
        predDF <- predDF %>%
          tidyr::gather(Prob.PredBacteria, Probability, -Source, -PredictedBacteria) %>%
          dplyr::group_by(Source, Prob.PredBacteria) %>%
          dplyr::summarise(Probability=mean(Probability)) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(Source) %>%
          dplyr::summarise(PredictedBacteria=lev[which.max(Probability)], Probability=Probability[1]) %>%
          dplyr::mutate(PredictedBacteria=factor(PredictedBacteria, levels=lev))
        return(predDF)
      }
    }else{
      predDF <- data.frame("Bacteria"=evalDF$"Bacteria", "Source"=evalDF$"Source", predDF)
      if(predictType=="Bacteria"){
        sink(paste0(outputHeader, "/ConfusionMatrix.txt"))
        print(caret::confusionMatrix(predDF$"PredictedBacteria", predDF$"Bacteria"))
        sink()
        classifierplots::classifierplots_folder(test.y=2-as.numeric(predDF$"Bacteria"), pred.prob=predDF[[lev[1]]], folder=outputHeader)
        return(predDF)
      }else if(predictType=="Image"){
        predDF <- predDF %>%
          tidyr::gather(Prob.PredBacteria, Probability, -Bacteria, -Source, -PredictedBacteria) %>%
          dplyr::group_by(Bacteria, Source, Prob.PredBacteria) %>%
          dplyr::summarise(Probability=mean(Probability)) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(Source) %>%
          dplyr::summarise(PredictedBacteria=lev[which.max(Probability)], Probability=Probability[1]) %>%
          dplyr::mutate(PredictedBacteria=factor(PredictedBacteria, levels=lev))
        sink(paste0(outputHeader, "/ConfusionMatrix.txt"))
        print(caret::confusionMatrix(predDF$"PredictedBacteria", predDF$"Bacteria"))
        sink()
        classifierplots::classifierplots_folder(test.y=2-as.numeric(predDF$"Bacteria"), pred.prob=predDF$"Probability", folder=outputHeader)
        return(predDF)
      }
    }
  }
  
  predDFList <- list()
  if(is.null(validDF)){
    predDFList[["Train"]] <- predEval(df_train, best_model, withReference=withReference, predictType=predictType, outputHeader="./Results/Training/")
    predDFList[["Test"]] <- predEval(df_test, best_model, withReference=withReference, predictType=predictType, outputHeader="./Results/HoldOutValidation/")
  }else{
    df_valid <- h2o::as.h2o(validDF)
    predDFList[["Train"]] <- predEval(df_train, best_model, withReference=withReference, predictType=predictType, outputHeader="./Results/Training/")
    predDFList[["Test"]] <- predEval(df_test, best_model, withReference=withReference, predictType=predictType, outputHeader="./Results/HoldOutValidation/")
    predDFList[["Valid"]] <- predEval(df_valid, best_model, withReference=withReference, predictType=predictType, outputHeader="./Results/ExternalValidation/")
  }
  gc();gc()
  return(predDFList)
}
