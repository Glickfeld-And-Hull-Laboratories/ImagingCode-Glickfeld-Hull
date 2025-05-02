#This script converts all tifs in a parent folder and all nested folders to svgs and saves them to a specified location
#It also automatically detects if a tif has already been converted and skips to the next one if so
#Written by TH 07/04/23

library(jpeg)
library(tiff)
library(tidyverse)
library(magick)

#change the path here to the folder that contains ALL subfolders with images you want to convert
files <- list.files(path='/Volumes/All_staff/home/jerry/data/slice_imaging/',pattern='.tif',recursive = TRUE,full.names = TRUE)

#change the path here to the folder where you want to save you jpegs
converted <- list.files(path='/Volumes/All_staff/home/jerry/data/slice_imaging/svgs/',pattern='.svg',full.names = FALSE)

#change the path here to the folder where you want to save you jpegs
setwd('/Volumes/All_staff/home/jerry/data/slice_imaging/svgs/')

for (file in files) {
  metdat <- str_split_1(file,pattern='/')
  ori_imgname <- metdat[11] #change the number here to the n-th item that contains your .tif file name
  newname <- gsub('.tif','.svg',ori_imgname)
  if (newname %in% converted){
    next} else {
      img <- image_read(file)
      my_svg <- image_convert(img, format="svg")
      image_write(my_svg, newname)
    }
}

