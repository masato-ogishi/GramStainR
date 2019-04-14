GramStainR: Automatic classification of bacterial species using Gram-stained images
===============================================

The 'GramStainR' package provides some easy-to-use functions for analyzing digital Gram stain images.

Installation
------------------------

-   Install the latest version from [GitHub](https://github.com/masato-ogishi/GramStainR) as follows:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/plotUtility")  ## internally used in GramStainR
devtools::install_github("masato-ogishi/GramStainR")
```

-   You might be prompted to install some packages before installling GramStainR. Follow the message(s).

Loading
------------------

``` r
library(tidyverse)
library(GramStainR)
```

Usage
-----------------------------------
``` r
# Working Environment
library(tidyverse)
library(GramStainR)
set.seed(12345)

# Detect bacteria
bact_ref <- bacteriaDetection(list.files(path="./Images/ReferenceImages/", pattern=".JPG", full.names=T))
saveRDS(bact_ref, "./Results/BacteriaDetection_Reference.rds")

bact_dibas <- bacteriaDetection(list.files(path="./Images/DIBaSImages/", pattern=".jpg", full.names=T))
saveRDS(bact_dibas, "./Results/BacteriaDetection_DIBaS.rds")

bact_clin <- bacteriaDetection(list.files(path="./Images/ClinicalImages/", pattern=".JPG", full.names=T))
saveRDS(bact_clin, "./Results/BacteriaDetection_Clinical.rds")

# Format bacteria datasets
bact_ref_PixelDF <- bacteriaDataFormat(bact_ref, spl=NULL, withReference=T)
saveRDS(bact_ref_PixelDF, "./Results/BacteriaData_Reference_All.rds")
bact_dibas_PixelDF <- bacteriaDataFormat(bact_dibas, spl=NULL, withReference=T)
saveRDS(bact_dibas_PixelDF, "./Results/BacteriaData_DIBaS_All.rds")
bact_clin_PixelDF <- bacteriaDataFormat(bact_clin, spl=NULL, withReference=T)
saveRDS(bact_clin_PixelDF, "./Results/BacteriaData_Clinical_All.rds")

bacteriaDataFormat_Split <- function(bactData, bactDataLabel, seedSet=1:5){
  for(s in seedSet){
    bact_PixelDFList <- bacteriaDataFormat(bactData, spl=0.7, seed=s, withReference=T)
    saveRDS(bact_PixelDFList, paste0("./Results/BacteriaData_", bactDataLabel, "_Seed", s, ".rds"))
  }
}
bacteriaDataFormat_Split(bact_ref, "Reference", seedSet=1:5)
bacteriaDataFormat_Split(bact_dibas, "DIBaS", seedSet=1:5)

# Train classifiers using the H2O AutoML instances
bacteriaClassification_AutoML_Wrapper <- function(bactDataLabel_TrainTest, seedSet=1:5){
  singleWrapper <- function(bactDataLabel_TrainTest, s=12345){
    bact_PixelDFList <- readRDS(paste0("./Results/BacteriaData_", bactDataLabel_TrainTest, "_Seed", s, ".rds"))
    leaderboardDF <- bacteriaClassification_AutoML(
      trainDF=bact_PixelDFList$"TrainDF",
      testDF=bact_PixelDFList$"TestDF",
      destDir="./Results/",
      LeaderBoardName=NULL,
      H2OModelName=paste0("BestH2OModel_", bactDataLabel_TrainTest, "_Seed", s),
      seed=s, max_mem_size="6G", nthreads=3
    ) %>% dplyr::mutate(DataLabel=bactDataLabel_TrainTest, Seed=s)
    return(leaderboardDF)
  }
  leaderboardDF <- data.table::rbindlist(lapply(seedSet, function(s){singleWrapper(bactDataLabel_TrainTest, s)}))
  readr::write_csv(leaderboardDF, paste0("./Results/LeaderBoard_", bactDataLabel_TrainTest, ".csv"))
  return(leaderboardDF)
}
bacteriaClassification_AutoML_Wrapper("Reference", seedSet=1:5)
bacteriaClassification_AutoML_Wrapper("DIBaS", seedSet=1:5)

# Apply classifiers to different datasets
bacteriaClassification_Evaluation_Wrapper <- function(bactDataLabel_model, bactDataLabel_eval, testQ=T, seedSet=1:5){
  testString <- ifelse(testQ, "_Test", "")
  singleWrapper <- function(bactDataLabel_model, bactDataLabel_eval, testQ=T, s=12345){
    if(testQ==T){
      bact_PixelDF <- readRDS(paste0("./Results/BacteriaData_", bactDataLabel_eval, "_Seed", s, ".rds"))[["TestDF"]]
    }else{
      bact_PixelDF <- readRDS(paste0("./Results/BacteriaData_", bactDataLabel_eval, "_All.rds"))
    }
    predDFList <- bacteriaClassification_Evaluation(
      evalDF=bact_PixelDF,
      targetBacteria="Paeruginosa",
      destDir="./Results/",
      H2OModelName=paste0("BestH2OModel_", bactDataLabel_model, "_Seed", s),
      PredictionHeader=paste0("Predictions_By_", bactDataLabel_model, "_", bactDataLabel_eval, testString, "_"),
      seed=s, max_mem_size="6G", nthreads=3
    )
    predDFList$"Bacteria" <- predDFList$"Bacteria" %>%
      dplyr::mutate(ModelLabel=bactDataLabel_model, DataLabel=paste0(bactDataLabel_eval, testString), Seed=s)
    predDFList$"Image" <- predDFList$"Image" %>%
      dplyr::mutate(ModelLabel=bactDataLabel_model, DataLabel=paste0(bactDataLabel_eval, testString), Seed=s)
    return(predDFList)
  }
  predDFList <- lapply(seedSet, function(s){singleWrapper(bactDataLabel_model, bactDataLabel_eval, testQ, s)})
  predDFList <- list(
    "Bacteria"=data.table::rbindlist(lapply(1:length(seedSet), function(i){predDFList[[i]][["Bacteria"]]})),
    "Image"=data.table::rbindlist(lapply(1:length(seedSet), function(i){predDFList[[i]][["Image"]]}))
  )
  readr::write_csv(predDFList$"Bacteria", paste0("./Results/Predictions_By_", bactDataLabel_model, "_", bactDataLabel_eval, testString, "_Bacteria.csv"))
  readr::write_csv(predDFList$"Image", paste0("./Results/Predictions_By_", bactDataLabel_model, "_", bactDataLabel_eval, testString, "_Image.csv"))
  return(predDFList)
}

bacteriaClassification_Evaluation_Wrapper("Reference", "Reference", testQ=T, seedSet=1:5); h2o::h2o.shutdown(F); gc(); gc()
bacteriaClassification_Evaluation_Wrapper("DIBaS", "DIBaS", testQ=T, seedSet=1:5); h2o::h2o.shutdown(F); gc(); gc()

bacteriaClassification_Evaluation_Wrapper("Reference", "Clinical", testQ=F, seedSet=1:5); h2o::h2o.shutdown(F); gc(); gc()
bacteriaClassification_Evaluation_Wrapper("DIBaS", "Clinical", testQ=F, seedSet=1:5); h2o::h2o.shutdown(F); gc(); gc()

# Performance
library(fbroc)
ROCPlot <- function(df, title=NULL){
  df <- df %>%
    dplyr::filter(Bacteria %in% c("Ecoli", "Paeruginosa")) %>%
    dplyr::mutate(
      Bacteria=factor(Bacteria, levels=c("Ecoli", "Paeruginosa"), labels=c("E.coli", "P.aeruginosa")),
      PredictedBacteria=factor(PredictedBacteria, levels=c("Ecoli", "Paeruginosa"), labels=c("E.coli", "P.aeruginosa"))
    )
  x <- df$Paeruginosa
  y <- dplyr::if_else(df$Bacteria=="P.aeruginosa", T, F)
  result.boot <- boot.roc(x, y, n.boot=100)
  plot(result.boot, show.metric="auc") +
    ggtitle(title) +
    plotUtility::theme_Publication(16)
}
violinPlot <- function(df, title=NULL){
  df <- df %>%
    dplyr::filter(Bacteria %in% c("Ecoli", "Paeruginosa")) %>%
    dplyr::mutate(
      Bacteria=factor(Bacteria, levels=c("Ecoli", "Paeruginosa"), labels=c("E.coli", "P.aeruginosa")),
      PredictedBacteria=factor(PredictedBacteria, levels=c("Ecoli", "Paeruginosa"), labels=c("E.coli", "P.aeruginosa"))
    )
  df$Paeruginosa <- scales::rescale(df$Paeruginosa)
  ggpubr::ggviolin(df, x="Bacteria", y="Paeruginosa", fill="Bacteria", 
                   add="boxplot", add.params=list(fill="white")) +
    ggpubr::stat_compare_means(aes(label=..p.signif..), 
                               comparisons=list(c("E.coli", "P.aeruginosa")),
                               label.y=1.2) +
    xlab(NULL) + ylab("Predicted probability") + 
    scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1.2)) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(3, "Dark2")[2:1]) +
    ggtitle(title) +
    plotUtility::theme_Publication(16) + theme(legend.position="none")
}

df_test_bact <- readr::read_csv("./Results/Predictions_By_Reference_Reference_Test_Bacteria.csv")
p1 <- ROCPlot(df_test_bact, "Bacteria-level")
p2 <- violinPlot(df_test_bact, "Bacteria-level")

df_test_img <- readr::read_csv("./Results/Predictions_By_Reference_Reference_Test_Image.csv")
p3 <- ROCPlot(df_test_img, "Image-level")
p4 <- violinPlot(df_test_img, "Image-level")

ggpubr::ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, 
                  labels=c("A","B","C","D"), font.label=list(size=18))
plotUtility::savePDF(
  o="./Results/Predictions_By_Reference_Reference_Test.pdf",
  w=10, h=10
)

df_clin_img <- readr::read_csv("./Results/Predictions_By_Reference_Clinical_Image.csv")
p5 <- ROCPlot(df_clin_img)
p6 <- violinPlot(df_clin_img)
ggpubr::ggarrange(p5, p6, ncol=2, nrow=1, 
                  labels=c("A","B"), font.label=list(size=18))
plotUtility::savePDF(
  o="./Results/Predictions_By_Reference_Clinical.pdf",
  w=10, h=5
)

df_dibas_test_img <- readr::read_csv("./Results/Predictions_By_DIBaS_DIBaS_Test_Image.csv")
p7 <- ROCPlot(df_dibas_test_img)
p8 <- violinPlot(df_dibas_test_img)
ggpubr::ggarrange(p7, p8, ncol=2, nrow=1, 
                  labels=c("A","B"), font.label=list(size=18))
plotUtility::savePDF(
  o="./Results/Predictions_By_DIBaS_DIBaS_Test.pdf",
  w=10, h=5
)

```
