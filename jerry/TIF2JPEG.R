#This script converts all tifs in a parent folder and all nested folders to jpegs and saves them to a specified location
#It also automatically detects if a tif has already been converted and skips to the next one if so
#09/12/23 now works with nested folders!
#Written by TH 07/03/23


library('jpeg')
library('tiff')
library('tidyverse')

#change the path here to the folder that contains ALL subfolders with images you want to convert
files <- list.files(path='/Volumes/All_staff/home/jerry/data/slice_imaging/',pattern='.tif',recursive = TRUE,full.names = TRUE)

#change the path here to the folder where you want to save you jpegs
converted <- list.files(path='/Volumes/All_staff/home/jerry/data/slice_imaging/jpegs/',pattern='.jpg',recursive = TRUE, full.names = FALSE)

if (length(converted)>0) {
  for (n in 1:length(converted)) {
    if (grepl('archive/',converted[n],fixed=TRUE) == TRUE){
      converted[n] <- gsub('archive/','',converted[n])
    } else {next}
  }
else{
  stop('No new images to be converted')
}
}

#change the path here to the folder where you want to save you jpegs
setwd('/Volumes/All_staff/home/jerry/data/slice_imaging/jpegs/')

for (file in files) {
  metdat <- str_split_1(file,pattern='/')
  ori_imgname <- metdat[11] #change the number here to the n-th item that contains your .tif file name
  imgid <- gsub('.tif','',ori_imgname)
  newname <- paste0(imgid,'.jpg')
  if (newname %in% converted){
    next} else {
      img <- readTIFF(file, native=TRUE)
      writeJPEG(img, target = newname, quality = 1)
  }
}


  