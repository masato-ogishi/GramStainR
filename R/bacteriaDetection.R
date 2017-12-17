#' Automatic bacteria detection.
#' 
#' @param originalImageFileName The input image file name. Should be provided as a combination of letters and numbers, since letters are interpreted as image source labels and numbers as image IDs. 
#' @param originalImageFileNameList A list of image file names.
#' @importFrom dplyr %>%
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom stringr str_split
#' @importFrom stringr fixed
#' @importFrom stringr str_replace_all
#' @importFrom EBImage readImage
#' @importFrom EBImage Image
#' @importFrom EBImage imageData
#' @importFrom EBImage getFrames
#' @importFrom EBImage combine
#' @importFrom EBImage otsu
#' @importFrom EBImage dilate
#' @importFrom EBImage makeBrush
#' @importFrom EBImage fillHull
#' @importFrom EBImage rotate
#' @importFrom EBImage flop
#' @importFrom EBImage as.Image
#' @importFrom EBImage display
#' @importFrom EBImage bwlabel
#' @importFrom EBImage computeFeatures.moment
#' @importFrom EBImage computeFeatures.shape
#' @importFrom OpenImageR edge_detection
#' @importFrom OpenImageR rgb_2gray
#' @importFrom imager as.cimg
#' @importFrom imager threshold
#' @importFrom imager clean
#' @importFrom data.table last
#' @importFrom tibble as_tibble
#' @importFrom abind adrop
#' @importFrom psych describe
#' @importFrom pbapply pblapply
#' @importFrom plotUtility savePDF
#' @export
#' @rdname bacteriaDetection
#' @name bacteriaDetection
bacteriaDetection_Single <- function(originalImageFileName){
  # Parameters
  windowSize=100
  ws=round(windowSize/2)
  brushSize=5
  thresholdAdjustRatio=1.2
  minEccentricity=0.6
  minAreaSize=300
  maxAreaSize=1000
  minOutlierFeature=0.1
  maxOutlierFeature=0.5
  #robustRange <- function(x, trim=0.2, n.sd=2){
  #  x.log.dist <- psych::describe(log10(x), trim=0.2)
  #  x.log.range <- x.log.dist$trimmed + x.log.dist$mad/0.675 * n.sd * c(-1, 1)  ## mean +/- 1sd => 67%, 2sd -> 95%
  #  x.range <- 10^x.log.range
  #  return(x.range)
  #}
  
  # Load the image data [RGB]
  ext <- data.table::last(unlist(stringr::str_split(basename(originalImageFileName), stringr::fixed("."))))
  header <- stringr::str_replace_all(originalImageFileName, paste0(".", ext, "$"), "")
  originalImageFileName <- paste0(header, ".", tolower(ext))
  img.Image <- EBImage::readImage(originalImageFileName)[,,1:3]
  
  # Thresholding
  img.thr <- EBImage::Image(EBImage::combine(lapply(EBImage::getFrames(img.Image), function(fr){
    img.neg <- 1-fr
    thr <- EBImage::otsu(img.neg)*thresholdAdjustRatio
    img.thr <- img.neg > thr
    img.thr <- img.thr %>% EBImage::dilate(EBImage::makeBrush(brushSize, shape="disc")) %>% EBImage::fillHull()
    img.thr <- 1- img.thr*img.neg
    return(img.thr)
  })), colormode='Color')
  
  # Detect edges
  img.array <- img.thr %>% EBImage::rotate(90) %>% EBImage::flop() %>% as.array()
  edg <- OpenImageR::edge_detection(img.array, method="Scharr", conv_mode="same")
  edg[1,,1:3] <- 0; edg[dim(img.array)[1],,1:3] <- 0; edg[,1,1:3] <- 0; edg[,dim(img.array)[2],1:3] <- 0
  edg <- OpenImageR::rgb_2gray(edg)
  
  # Create an object mask based on edges
  msk <- imager::as.cimg(edg) %>% imager::threshold() %>% imager::clean(1) %>% 
    EBImage::as.Image() %>% EBImage::dilate(EBImage::makeBrush(brushSize, shape="disc")) %>% EBImage::fillHull() %>% 
    EBImage::rotate(90) %>% EBImage::flop()
  if(length(dim(msk))>=4){
    dm <- 1:length(dim(msk))
    dm <- rev(dm[-1:-3])
    for(d in dm){
      msk <- abind::adrop(msk, d)
    }
  }
  img.r <- dim(img.Image) %>% (function(X){X[2]/X[1]})
  img.msk <- EBImage::Image(EBImage::combine(lapply(EBImage::getFrames(img.Image), function(fr){fr*msk})), colormode='Color')
  EBImage::display(img.msk, method="raster")
  suppressMessages(plotUtility::savePDF(outputFileName=paste0(header, "_Masked.pdf"), w=8, h=8*img.r))
  
  # Label segmented objects
  img.bw <- EBImage::bwlabel(msk)
  #img.seg <- EBImage::paintObjects(img.bw, img.Image, col=c('#007eff', NA))
  #EBImage::display(img.seg, method="raster")
  #suppressMessages(plotUtility::savePDF(outputFileName=paste0(header, "_Segmented.pdf"), w=8, h=8*img.r))
  
  # Filter objects based on their morphological features
  featureDF <- dplyr::bind_cols(
    as.data.frame(EBImage::computeFeatures.moment(img.bw)),
    as.data.frame(EBImage::computeFeatures.shape(img.bw))
  ) %>% 
    dplyr::mutate(ObjID=1:nrow(.)) %>%
    dplyr::mutate(m.theta=m.theta*180/pi)   ## convert radian to degree
  featureDF[["OutlierFeature"]] <- (featureDF$s.perimeter/featureDF$s.area)*(featureDF$s.radius.max/featureDF$s.radius.min)
  featureDF <- featureDF %>%
    dplyr::filter(m.majoraxis<=windowSize) %>%  
    dplyr::filter(minEccentricity<=m.eccentricity) %>% 
    dplyr::filter(minAreaSize<=s.area & s.area<=maxAreaSize) %>%   
    dplyr::filter(minOutlierFeature<=OutlierFeature & OutlierFeature<=maxOutlierFeature) %>%
    dplyr::filter(1<=round(m.cx-ws) & round(m.cx+ws)<=dim(img.bw)[1] & 1<=round(m.cy-ws) & round(m.cy+ws)<=dim(img.bw)[2]) %>%
    dplyr::select(ObjID, setdiff(colnames(.), "ObjID"))
  if(nrow(featureDF)==0) return(NULL)
  
  # Retrieve single bacteria images
  singleObjectImageList <- as.list(as.numeric(featureDF$"ObjID"))
  for(r in 1:nrow(featureDF)){
    i <- singleObjectImageList[[r]]
    msk.sobj <- msk*(img.bw==i)
    ft <- dplyr::filter(featureDF, ObjID==i)
    msk.sobj <- msk.sobj[round(ft[["m.cx"]]+(-ws:ws)), round(ft[["m.cy"]]+(-ws:ws))]
    msk.sobj <- EBImage::rotate(msk.sobj, -ft[["m.theta"]], bg.col="black")
    dm <- dim(msk.sobj)[1]
    rm.dm <- dm-windowSize
    msk.sobj <- msk.sobj[seq(floor(rm.dm/2)+1, floor(dm-rm.dm/2)), seq(floor(rm.dm/2)+1, floor(dm-rm.dm/2))]
    singleObjectImageList[[r]] <- msk.sobj
  }
  
  # Generate a combined view
  singleObjectImageCombined <- EBImage::combine(singleObjectImageList)
  EBImage::display(singleObjectImageCombined, all=T, method="raster")
  suppressMessages(plotUtility::savePDF(outputFileName=paste0(header, "_SingleObjectImageCombined.pdf"), w=8, h=8))
  
  # Convert to pixel matrices
  singleBacteriaImageToPixelVector <- function(singleBacteriaImage){
    if(EBImage::colorMode(singleBacteriaImage)=="Color") singleBacteriaImage <- singleBacteriaImage[,,1:3] 
    img <- EBImage::imageData(singleBacteriaImage)
    img <- as.vector(img)
    return(img)
  }
  dm <- dim(singleObjectImageList[[1]])
  pixelMatrix <- matrix(nrow=length(singleObjectImageList), ncol=dm*dm)
  for(i in 1:length(singleObjectImageList)){
    pixelMatrix[i, ] <- singleBacteriaImageToPixelVector(singleObjectImageList[[i]])
  }
  pixelDF <- data.frame("ObjID"=featureDF$"ObjID", as.data.frame(pixelMatrix))
  
  # Output
  imageSource <- sub(pattern="(.*)\\..*$", replacement="\\1", basename(originalImageFileName))
  featureDF <- tibble::as_tibble(data.frame("Source"=imageSource, featureDF))
  pixelDF <- tibble::as_tibble(data.frame("Source"=imageSource, pixelDF))
  res <- list("FeatureDF"=featureDF, "PixelDF"=pixelDF, "SingleObjectImageSet"=singleObjectImageList)
  rm(list=setdiff(ls(), "res"))
  gc();gc()
  return(res)
}

#' @export
#' @rdname bacteriaDetection
#' @name bacteriaDetection
bacteriaDetection <- function(originalImageFileNameList){
  imageSourceList <- lapply(
    originalImageFileNameList, 
    function(originalImageFileName){sub(pattern="(.*)\\..*$", replacement="\\1", basename(originalImageFileName))}
  )
  res <- pbapply::pblapply(originalImageFileNameList, bacteriaDetection_Single)
  names(res) <- imageSourceList
  return(res)
}
