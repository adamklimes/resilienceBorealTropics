## data preparation
library(raster)
library(terra)

## data_preparation
prepVegFiles <- function(fold, seasonDates, expSeasonDates, extArea, thres, nameAdd = NULL, doLC = TRUE){
  foldVI <- paste0(fold, "VI_16Days_1Km_v61/")
  eviRawList <- lapply(paste0(foldVI, "EVI/", dir(paste0(foldVI, "EVI/"))), rast)
  QARawList <- lapply(paste0(foldVI, "QA_qual/", dir(paste0(foldVI, "QA_qual/"))), rast)
  eviExt <- vapply(eviRawList, function(x) as.character(round(ext(x), 3)), FUN.VALUE = "aa")
  eviRaw <- rast(eviRawList[eviExt == names(sort(table(eviExt), decreasing = TRUE)[1])])
  QAnamesFull <- vapply(QARawList, sources, FUN.VALUE = "aa")
  QAnames <- sub("QA_qual", "EVI", substr(QAnamesFull, nchar(QAnamesFull) - 27,
    nchar(QAnamesFull) - 4))
  QARaw <- rast(QARawList[QAnames %in% names(eviRaw)])

  eviC <- crop(eviRaw, extent(extArea))
  QAC <- crop(QARaw, extent(extArea))
  excludeQA <- function(dat, QA) dat + log10((QA < 0.5) * 2 - 1)
  eviQA <- rast(mapply(excludeQA, as.list(eviC), as.list(QAC)))
  writeRaster(eviQA, paste0(fold, "eviQA.tiff"))
  # eviQA <- rast(paste0(fold, "eviQA.tiff"))
  year <- as.numeric(regmatches(names(eviQA), regexpr(paste0("[[:digit:]]{4}"), names(eviQA))))
  day <- regmatches(names(eviQA), regexpr(paste0("[[:digit:]]{3}$"), names(eviQA)))

  eviSeason <- eviQA[[day %in% seasonDates]]
  eviSeasonM <- tapp(eviSeason, (year[day %in% seasonDates]) > 2010,
    function(x) mean(x, na.rm = TRUE) + log10((sum(is.na(x)) < thres[1]) * 2 - 1))
  names(eviSeasonM) <- paste0("eviM", nameAdd, c("Old", "New"))
  writeRaster(eviSeasonM, paste0(fold, "eviSeason", nameAdd, ".tiff"))

  eviYear <- eviQA[[day %in% expSeasonDates]]
  eviYearSD <- tapp(eviYear, (year[day %in% expSeasonDates]) > 2010, function(x)
    sd(x, na.rm = TRUE) + log10((sum(is.na(x)) < thres[2]) * 2 - 1))
  names(eviYearSD) <- paste0("eviSD", nameAdd, c("Old", "New"))
  writeRaster(eviYearSD, paste0(fold, "eviSD", nameAdd, ".tiff"))

  if (doLC){
    foldLC <- paste0(fold, "LandCover_Type_Yearly_500m_v6/")
    lcRaw <- rast(paste0(foldLC, "LC2/", dir(paste0(foldLC, "LC2/"))))
    lcA <- aggregate(lcRaw, fact = 2, fun = "modal", ties = "random", na.rm = TRUE)
    lcC <- crop(lcA, extent(extArea))
    yearLC <- as.numeric(regmatches(names(lcC), regexpr(paste0("[[:digit:]]{4}"), names(lcC))))
    lc <- lcC[[yearLC %in% c(2002, 2019)]]
    names(lc) <- c("lcOld", "lcNew")
    writeRaster(lc, paste0(fold, nameAdd, "lc.tiff"))
  }
}

# RU
prepVegFiles(
  fold = "data/modis/RU/",
  seasonDates = c("161","177","193","209","225","241"), # Julian days for June-August
  expSeasonDates = c("097","113","129","145","161","177","193","209","225","241","257","273","289"), # Julian days for April-October
  extArea = c(68,80,45,74),
  thres = c(10, 36) # number of missing values for a pixel to be excluded from EVI mean and var
)

# SA
prepVegFiles(
  fold = "data/modis/SA/",
  seasonDates = c("337","353","001","017","033","049"), # Julian days for December-February
  expSeasonDates = c("289","305","321","337","353","001","017","033","049","065","081","097","113"), # Julian days for October-April
  extArea = c(-62,-50,-29,0),
  thres = c(28, 54) # number of missing values for a pixel to be excluded from EVI mean and var
)
# SA Supp
prepVegFiles(
  fold = "data/modis/SA/",
  seasonDates = c("161","177","193","209","225","241"), # Julian days for June-August
  expSeasonDates = c("097","113","129","145","161","177","193","209","225","241","257","273","289"), # Julian days for April-October
  extArea = c(-62,-50,-29,0),
  thres = c(28, 54), # number of missing values for a pixel to be excluded from EVI mean and var
  nameAdd = "Supp", doLC = FALSE
)

# clim
precRU <- rast(grep("dat06|dat07|dat08", paste0("data/clim/pr/", dir("data/clim/pr/")), value = TRUE)) # in kg*m-2*-month/100
tempRU <- rast(grep("dat06|dat07|dat08", paste0("data/clim/tas/", dir("data/clim/tas/")), value = TRUE)) # in K/10
precSA <- rast(grep("dat01|dat02|dat12", paste0("data/clim/pr/", dir("data/clim/pr/")), value = TRUE))
tempSA <- rast(grep("dat01|dat02|dat12", paste0("data/clim/tas/", dir("data/clim/tas/")), value = TRUE))
# plot(temp[[1]]/10 - 273.15) # in °C
eviRU <- rast("data/modis/RU/eviSeason.tiff")
eviSA <- rast("data/modis/SA/eviSeason.tiff")

selRU <- list(Full = "_19(8|9)[0-9]_", Old = "_199[5-9]_", New = "_201[2-6]_")
selSA <- list(Full = "_19(8|9)[0-9]_", Old = "_199[6-9]_|_12_1995_|_0(1|2)_2000_", New = "_201[3-6]_|_12_2012_|_0(1|2)_2017_")

prepClim <- function(dat, selDates, resFile, fn = mean, newName = NULL){
  n <- substitute(dat)
  out <- lapply(selDates, function(x) fn(dat[[grep(x, names(dat))]]))
  out <- resample(rast(out), resFile)
  names(out) <- paste0(substr(n, 1, nchar(n) - 2), names(out))
  if (!is.null(newName)) names(out) <- newName
  writeRaster(out, paste0("data/clim/", substitute(dat), ".tiff"))
  out
}

precRUL <- prepClim(precRU, selRU, eviRU)
tempRUL <- prepClim(tempRU, selRU, eviRU)
precSAL <- prepClim(precSA, selSA, eviSA)
tempSAL <- prepClim(tempSA, selSA, eviSA)

# SA Supp
precSuppSA <- rast(grep("dat06|dat07|dat08", paste0("data/clim/prSupp/", dir("data/clim/prSupp/")), value = TRUE)) # in kg*m-2*-month/100
tempSuppSA <- rast(grep("dat06|dat07|dat08", paste0("data/clim/tasSupp/", dir("data/clim/tasSupp/")), value = TRUE)) # in K/10

precSuppSAL <- prepClim(precSuppSA, selRU[2:3], eviSA)
tempSuppSAL <- prepClim(tempSuppSA, selRU[2:3], eviSA)

#_
