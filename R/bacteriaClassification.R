#' Classification of bacteria based on their morphological features using H2O AutoML instances
#'
#' @param trainDF A dataframe for model training.
#' @param testDF A dataframe for hold-out validation.
#' @param evalDF A dataframe for evaluation.
#' @param targetBacteria The name of the bacteria of interest.
#' @param destDir A directory for saving outputs.
#' @param LeaderBoardName The file name of the leader board to be saved. The file format must be CSV. If you want to skip saving, set \code{NULL}.
#' @param H2OModelName The file name of the best H2O model to be saved. No file extension is required. If you want to skip saving, set \code{NULL}.
#' @param PredictionHeader The file name header of the prediction results to be saved.
#' @param seed A random seed.
#' @param max_mem_size The upper limit of memory for H2O Java machine learning engine.
#' @param nthreads The number of threads for H2O Java machine learning engine.
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
#' @importFrom readr write_csv
#' @importFrom caret downSample
#' @importFrom pROC multiclass.roc
#' @importFrom pROC coords
#' @importFrom plotROC geom_roc
#' @importFrom plotROC style_roc
#' @importFrom plotROC calc_auc
#' @importFrom plotUtility theme_Publication
#' @importFrom plotUtility savePDF
#' @importFrom ggsci pal_d3
#' @import h2o
#' @import ggplot2
#' @export
#' @rdname bacteriaClassification
#' @name bacteriaClassification
bacteriaClassification_AutoML <- function(
  trainDF, testDF,
  destDir="./Results/", LeaderBoardName="LeaderBoard.csv", H2OModelName="BestH2OModel",
  seed=12345, max_mem_size="6G", nthreads=6
){
  set.seed(seed)
  dir.create(destDir, showWarnings=F, recursive=T)
  h2o::h2o.init(ip="localhost", port=seed, max_mem_size=max_mem_size, nthreads=nthreads)

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
  if(!is.null(LeaderBoardName)) readr::write_csv(df_lb, file.path(destDir, LeaderBoardName))

  # Save the best classifier
  best_model_name <- df_lb[["model_id"]][[1]]
  best_model <- h2o::h2o.getModel(best_model_name)
  if(!is.null(H2OModelName)){
    h2o::h2o.saveModel(object=best_model, path=destDir, force=T)
    invisible(file.copy(
      from=file.path(destDir, best_model_name),
      to=file.path(destDir, H2OModelName),
      overwrite=T
    ))
    invisible(file.remove(file.path(destDir, best_model_name)))
  }

  return(df_lb)
}

#' @export
#' @rdname bacteriaClassification
#' @name bacteriaClassification
bacteriaClassification_Prediction <- function(
  evalDF, destDir="./Results/", H2OModelName="BestH2OModel", seed=12345, max_mem_size="6G", nthreads=6
){
  set.seed(seed)
  dir.create(destDir, showWarnings=F, recursive=T)
  h2o::h2o.init(ip="localhost", port=seed, max_mem_size=max_mem_size, nthreads=nthreads)

  evalH2ODF <- h2o::as.h2o(evalDF)
  H2OMod <- h2o::h2o.loadModel(file.path(destDir, H2OModelName))
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

#' @export
#' @rdname bacteriaClassification
#' @name bacteriaClassification
bacteriaClassification_Evaluation <- function(
  evalDF, targetBacteria="Paeruginosa", destDir="./Results/", H2OModelName="BestH2OModel", PredictionHeader="Predictions_",
  seed=12345, max_mem_size="6G", nthreads=6
){
  set.seed(seed)
  dir.create(destDir, showWarnings=F, recursive=T)
  h2o::h2o.init(ip="localhost", port=seed, max_mem_size=max_mem_size, nthreads=nthreads)

  predDFList <- bacteriaClassification_Prediction(evalDF, destDir, H2OModelName, seed, max_mem_size, nthreads)
  saveRDS(predDFList, file.path(destDir, paste0(PredictionHeader, "Results.rds")))

  cm <- caret::confusionMatrix(predDFList$"Bacteria"$"PredictedBacteria", predDFList$"Bacteria"$"Bacteria")
  print(cm)
  sink(file.path(destDir, paste0(PredictionHeader, "BacteriaLevel_ConfusionMatrix_Seed", seed, ".txt")))
  print(cm)
  sink()

  lev <- levels(evalDF$"Bacteria")
  colPal <- ggsci::pal_d3()(length(lev))
  colPal[grep(targetBacteria, lev, value=F)] <- "black"
  thrList <- suppressWarnings(lapply(pROC::multiclass.roc(predDFList$"Image"$"Bacteria", predDFList$"Image"[[targetBacteria]])$rocs, function(r){pROC::coords(r, "best", ret="threshold")})) ## Two or more threshold could sometimes be returned... The first one is the lowest, meaning maximum sensitivity for P.aeruginosa.
  probPlot <- ggplot(data=predDFList$"Image")
  for(i in 1:length(thrList)){
    probPlot <- probPlot + geom_hline(yintercept=thrList[[i]][1], color=setdiff(colPal, "black")[i], size=1)
  }
  probPlot <- probPlot +
    geom_jitter(aes_string(x="Bacteria", y=targetBacteria, color="Bacteria"), width=0.2) +
    xlab(NULL) + scale_color_manual(values=colPal, guide=F) +
    plotUtility::theme_Publication()
  plotUtility::savePDF(probPlot, outputFileName=file.path(destDir, paste0(PredictionHeader, "ImageLevel_JitterPlot_Seed", seed, ".pdf")))

  if(length(lev)==2){
    rocPlot <- predDFList$"Image" %>%
      dplyr::mutate(Bacteria=as.numeric(Bacteria)-1) %>%
      ggplot(aes_string(d="Bacteria", m=targetBacteria)) +
      plotROC::geom_roc(n.cuts=0) +
      plotROC::style_roc() +
      plotUtility::theme_Publication()
    rocPlot <- rocPlot + annotate("text", x=.75, y=.25, size=8, label=paste("AUC =", round(plotROC::calc_auc(rocPlot)$"AUC", 2)))
    plotUtility::savePDF(rocPlot, outputFileName=file.path(destDir, paste0(PredictionHeader, "ImageLevel_ROCPlot_Seed", seed, ".pdf")))
  }

  return(predDFList)
}
