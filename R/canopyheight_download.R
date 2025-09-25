## canopyheight
library(terra)
# https://langnico.github.io/globalcanopyheight/
# Lang, N., Jetz, W., Schindler, K., & Wegner, J. D. (2023). A high-resolution canopy height model of the Earth. Nature Ecology & Evolution, 1-12.

## data download
getCHM <- function(coors){
  fileN <- paste0("ETH_GlobalCanopyHeight_10m_2020_", coors, "_Map.tif")
  urlETH <- paste0("https://libdrive.ethz.ch/index.php/s/cO8or7iOe5dT2Rt/download?path=%2F3deg_cogs&files=", fileN)
  download.file(urlETH, paste0(getwd(),"/data/CHM/", fileN), mode = "wb")
  print(coors)
  NULL
}
fillZ <- function(x, n = 2){
  z <- pmax(n - nchar(x), 0)
  paste0(vapply(z, function(x) paste(rep("0", x), collapse = ""), FUN.VALUE = "00"), x)
}

cRU <- expand.grid(fillZ(15:24*3), fillZ(22:26*3, 3))
coorsRU <- paste0("N", cRU$Var1, "E", cRU$Var2)
cSA <- expand.grid(fillZ(1:10*3), fillZ(17:21*3, 3))
coorsSA <- paste0("S", cSA$Var1, "W", cSA$Var2)

# lapply(coorsRU, getCHM)
# lapply(coorsSA, getCHM)

## data preparation
prepCHM <- function(coors, refSet, name){
  chmVrt <- vrt(paste0("data/CHM/ETH_GlobalCanopyHeight_10m_2020_", coors, "_Map.tif"))
  chmVrtA <- aggregate(chmVrt, fact = 70, na.rm = TRUE)
  chmVrtAR <- resample(chmVrtA, refSet)
  names(chmVrtAR) <- "chm"
  writeRaster(chmVrtAR, paste0("data/CHM/", name, ".tif"))
}
eviRU <- rast("data/modis/RU/eviSeason.tiff")
eviSA <- rast("data/modis/SA/eviSeason.tiff")

# prepCHM(coorsRU, eviRU, "chmRU")
# prepCHM(coorsSA, eviSA, "chmSA")

#_
