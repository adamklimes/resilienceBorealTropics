## climate_download
library(terra)

## data_download
getPrec <- function(selTime, selVar, extArea, folderEnd = NULL){
  urlChelsa <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/",
  selVar, "/CHELSA_", selVar, "_", selTime, "_V.2.1.tif")
  out <- crop(rast(urlChelsa), extArea)
  writeRaster(out, paste0("data/clim/", selVar, folderEnd, "/dat", selTime, ".tif"))
  print(selTime)
  NULL
}

# RU
selMonths <- paste0("0", 6:8)
selTime <- c(paste(selMonths, rep(1980:1999, each = length(selMonths)), sep = "_"),
  paste(selMonths, rep(2012:2016, each = length(selMonths)), sep = "_"))

sapply(selTime, getPrec, selVar = "pr", extArea = c(68,80,45,74))
sapply(selTime, getPrec, selVar = "tas", extArea = c(68,80,45,74))

#SA
selMonths <- c("01", "02", "12")
selTime <- c(paste(selMonths, rep(1980:1999, each = length(selMonths)), sep = "_"),
  paste(selMonths, rep(2012:2017, each = length(selMonths)), sep = "_")[-c(1,2,18)],
  "01_2000", "02_2000")

sapply(selTime, getPrec, selVar = "pr", extArea = c(-62,-50,-29,0))
sapply(selTime, getPrec, selVar = "tas", extArea = c(-62,-50,-29,0))

#SAsupp
selMonths <- paste0("0", 6:8)
selTime <- c(paste(selMonths, rep(1995:1999, each = length(selMonths)), sep = "_"),
  paste(selMonths, rep(2012:2016, each = length(selMonths)), sep = "_"))

sapply(selTime, getPrec, selVar = "pr", extArea = c(-62,-50,-29,0), folderEnd = "Supp")
sapply(selTime, getPrec, selVar = "tas", extArea = c(-62,-50,-29,0), folderEnd = "Supp")

#_
