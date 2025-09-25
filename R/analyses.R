## analyses
library(mixglm)
library(raster)
library(terra)

# data_loading
climFiles <- c("prec", "temp", "precSupp", "tempSupp")
vegFiles <- c("eviSeason", "eviSD", "lc", "eviSeasonSupp", "eviSDSupp")
filesRU <- rast(c(paste0("data/clim/", climFiles[1:2], "RU.tiff"), paste0("data/modis/RU/", vegFiles[1:3], ".tiff"), "data/chm/chmRU.tif"))
filesSA <- rast(c(paste0("data/clim/", climFiles, "SA.tiff"), paste0("data/modis/SA/", vegFiles, ".tiff"), "data/chm/chmSA.tif"))

# data_preparation
st <- function(x, y = x) (x - mean(y)) / sd(y)
stData <- function(files, R = 4000, fixedID = NULL){
  vals <- terra::values(files, dataframe = TRUE)
  selID <- which(rowSums(is.na(vals)) == 0 &
    (!vals$lcOld %in% c(0,11:13,15)) & (!vals$lcNew %in% c(0,11:13,15)))
  set.seed(10)
  sampID <- sample(selID, R)
  if (!is.null(fixedID)) sampID <- fixedID
  datSel <- vals[sampID, ]
  datSelSt <- data.frame(apply(datSel, 2, st))
  list(sampID = sampID, datSel = datSel, datSelSt = datSelSt)
}
getOrigID <- function(files, coor){
  (rowFromY(files, coor$ycoor)-1) * (dim(files)[2]) + colFromX(files, coor$xcoor) - 1
}

## loading previously calculated model outputs
load(file = "data/analyses/modOuts.RData")
  # fixedID to recover original sampID
datRU <- stData(filesRU, fixedID = getOrigID(filesRU, modOuts[[1]]$coor))
datSA <- stData(filesSA, fixedID = getOrigID(filesSA, modOuts[[20]]$coor))

## analyses_functions
modsPars <- data.frame(
  transect = c(rep(c("RU", "SA"), each = 8), rep("SA", 4)),
  respType = c(rep(c("eviM", "eviSD"), times = 2, each = 4), rep(c("eviMSupp", "eviSDSupp"), each = 2)),
  respCat = c(rep(c("Old", "New"), times = 4, each = 2), rep(c("Old", "New"), 2)),
  predsType = c(rep(c("Full", "Old", "Full", "New"), 4), rep(c("SuppOld", "SuppNew"), 2))
)
runMod <- function(pars){
  formPred <- formula(paste0("~ prec", pars["predsType"], " + temp", pars["predsType"]))
  form <- formula(paste0(pars["respType"], pars["respCat"], " ", paste(as.character(formPred), collapse = " ")))
  inputData <- paste0("dat", pars["transect"])
  numStates <- 6
  setInit = list(intercept_stateVal = c(-1, rep(0.5, numStates - 1)),
                 intercept_statePrec = rep(2, numStates),
                 intercept_stateProb = c(0, rep(0.01, numStates - 1)),
                 precX_stateVal = rep(0.01, numStates),
                 precX_statePrec = rep(0.01, numStates),
                 precX_stateProb = rep(0.01, numStates),
                 tempX_stateVal = rep(0.01, numStates),
                 tempX_statePrec = rep(0.01, numStates),
                 tempX_stateProb = rep(0.01, numStates))
  setPriors = list(statePrec = list(int = "dnorm(0, 0.1)", pred = "dnorm(0, 0.1)"))
  names(setInit) <- paste0(rep(c("intercept", "prec", "temp"), each = 3),
    rep(c("", pars["predsType"], pars["predsType"]), each = 3), c("_stateVal", "_statePrec", "_stateProb"))
  set.seed(11)
  mod <- mixglm(
    stateValModels = form,
    stateProbModels = formPred,
    statePrecModels = formPred,
    inputData = get(inputData)$datSelSt,
    numStates = numStates,
    mcmcChains = 2,
    setInit = setInit,
    setPriors = setPriors
  )
  save(mod, file = paste0("data/analyses/mod", paste(pars, collapse = ""), ".RData"))
  mod
}
# model_run
# done: 1:20
# runMod(modsPars[1, ])
# apply(modsPars, 1, runMod)

# figure outputs
assignColor <- function(dat, col){
  if (all(is.na(dat))) return(NA)
  sq <- seq(min(dat, na.rm = TRUE), max(dat, na.rm = TRUE), length.out = length(col)+1)
  col[findInterval(dat, sq, all.inside = TRUE)]
}
invSt <- function(x, y) x * sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)
calcDistState <- function(ID, pred, dat, selVar, cf = 1){
  states <- pred$tipStable[[ID]]
  selStates <- states[states$catSt == 1, ]
  selStates$dist <- (invSt(selStates$resp, dat$datSel[, selVar]) - dat$datSel[, selVar][ID]) * cf
  tips <- selStates[selStates$state == 0, ]
  if (nrow(tips) == 0) tips <- setNames(data.frame(NA,NA,NA,NA,NA), names(tips))
  distToTip <- tips$dist[abs(tips$dist) == min(abs(tips$dist))]
  if (is.na(distToTip)) distToTip <- Inf
  distToStates <- selStates$dist[selStates$state == 1]
  selDists <- distToStates[abs(distToStates) < abs(distToTip) | sign(distToStates) != sign(distToTip)]
  selDists[abs(selDists) == min(abs(selDists))]
}
calcOutputs <- function(mod, transect, thres = 0.1){
  pred <- predict(mod, threshold = thres)
  sampID <- get(paste0("dat", transect))$sampID
  files <- get(paste0("files", transect))
  coor <- data.frame(xcoor = xFromCol(files)[(sampID %% dim(files)[2]) + 1], ycoor = yFromRow(files)[1 + trunc(sampID/dim(files)[2])])
  colsPotE <- assignColor(pred$obsDat$potentEn, rev(heat.colors(130)[1:100])) # higher values are red
  colsDistT <- assignColor(pred$obsDat$distToTip, heat.colors(130)[1:100]) # lower values are red
  distToState <- vapply(1:ncol(pred$probCurves), calcDistState, pred, get(paste0("dat", transect)), names(mod$data), FUN.VALUE = 1.0)
  colsDistS <- assignColor(abs(distToState), rev(heat.colors(130)[1:100])) # higher values are red
  list(coor = coor, colsPotE = colsPotE, colsDistT = colsDistT, pred = pred, modPars = c(transect, names(mod$data)), distToState = distToState, colsDistS = colsDistS)
}
calcModOut <- function(i, modsPars){
  load(file = paste0("data/analyses/mod", paste(modsPars[i, ], collapse = ""), ".RData"))
  print(i)
  calcOutputs(mod, modsPars$transect[i], 0.1)
}

# modOuts <- lapply(1:nrow(modsPars), calcModOut, modsPars)
# save(modOuts, file = "data/analyses/modOuts.RData")
# load(file = "data/analyses/modOuts.RData")

# i <- 1
# load(file = paste0("data/analyses/mod", paste(modsPars[i, ], collapse = ""), ".RData"))
# modOut <- calcOutputs(mod, modsPars$transect[i], 0.1)

# change ~ potE
calcEffPotE <- function(modOut, ylab = NULL, cf = 1, pos, axis1 = -2:2/5, sepPlots = FALSE, xlab = TRUE){
  transect <- modOut$modPars[1]
  respType <- sub("Old|New", "", modOut$modPars[2])
  charNum <- function(x, ndec = 2, pval = FALSE){
    if (x == 0 & pval) return("< 0.001")
    x <- as.character(x)
    if (!grepl("\\.", x)) x <- paste0(x, ".")
    n <- nchar(x)
    paste0(x, paste(rep(0, ndec + 2 - n), collapse = ""))
  }
  potE <- modOut$pred$obsDat$potentEn
  potEst <- st(potE)
  dat <- get(paste0("dat", transect))
  respNames <- paste0(respType, c("Old", "New"))
  change <- (dat$datSel[, respNames[2]] - dat$datSel[, respNames[1]]) * cf
  if (is.null(ylab)) ylab <- paste("Change in", respType)
  if (!sepPlots) par(mfrow = c(1, 2), mai = c(1,0.8,0.1,0.1)) else par(mai = par("mai")[c(1,4,3,2)])
  plot(change ~ potEst, xlab = "", ylab = ylab,
       axes = FALSE, cex = 0.5)
  if (xlab) axis(1, labels = "Potential energy", at = mean(par("usr")[1:2]), tick = FALSE, line = 2)
  axis(1, labels = 0:2/2, at = st(0:2/2, potE))
  axis(2, las = 2)
  #if (modOut$modPars[1] == "RU") box(col = trcol[1], lwd = 3) else box(col = trcol[2], lwd = 3)
  box()
  m <- lm(change ~ potEst)
  xx <- seq(min(potEst), max(potEst), length.out = 100)
  pr <- predict(m, newdata = list(potEst = xx), se.fit = TRUE)
  lines(xx, pr$fit, lwd = 2, col = "red")
  lines(xx, pr$fit + pr$se.fit, lwd = 2, lty = 2, col = "red")
  lines(xx, pr$fit - pr$se.fit, lwd = 2, lty = 2, col = "red")
  expVar <- charNum(round(summary(m)$adj.r.squared, 2))
  pval <- charNum(round(summary(m)$coef[2, 4], 3), 3, pval = TRUE)
  text(pos[1], pos[2],
    bquote("Adj-"*R^2*": "*.(expVar)), adj = c(0,0), col = "red", cex = 0.8)
  text(pos[1], pos[2] - pos[3],
       bquote("P-valu"*e["slope"]*": "*.(pval)), adj = c(0,0), col = "red", cex = 0.8)

  distSS <- modOut$distToState * cf
  distSSst <- st(distSS)
  if (!sepPlots) par(mai = c(1,0.1,0.1,0.8)) else par(mai = par("mai")[c(1,4,3,2)])
  plot(change ~ distSSst, xlab = "",
    ylab = "", axes = FALSE, cex = 0.5)
  if (xlab) axis(1, labels = "Dist. to stable state", at = mean(par("usr")[1:2]), tick = FALSE, line = 2)
  axis(1, labels = axis1, at = st(axis1, distSS))
  #if (modOut$modPars[1] == "RU") box(col = trcol[1], lwd = 3) else box(col = trcol[2], lwd = 3)
  box()
  m2 <- lm(change ~ distSSst)
  xx <- seq(min(distSSst), max(distSSst), length.out = 100)
  pr <- predict(m2, newdata = list(distSSst = xx), se.fit = TRUE)
  lines(xx, pr$fit, lwd = 2, col = "red")
  lines(xx, pr$fit + pr$se.fit, lwd = 2, lty = 2, col = "red")
  lines(xx, pr$fit - pr$se.fit, lwd = 2, lty = 2, col = "red")
  expVar <- charNum(round(summary(m2)$adj.r.squared, 2))
  pval <- charNum(round(summary(m2)$coef[2, 4], 3), 3, pval = TRUE)
  text(pos[4], pos[5],
       bquote("Adj-"*R^2*": "*.(expVar)), adj = c(0,0), col = "red", cex = 0.8)
  text(pos[4], pos[5] - pos[6],
       bquote("P-valu"*e["slope"]*": "*.(pval)), adj = c(0,0), col = "red", cex = 0.8)
  list(m, m2)
}

## figures
modsParsS <- apply(modsPars, 1, paste, collapse = "")
addCurve <- function(modOut, id, pos, addCols = TRUE){
  cols <- if(addCols) c("black", "grey40") else c("grey50", "grey70")
  selVar <- modOut$modPars[2]
  stXcoor <- function(x, vals = modOut$pred$sampledResp, tr = function(x) x*5-2.5) tr((x - min(vals)) / max(vals - min(vals)))
  lines(stXcoor(modOut$pred$sampledResp) + pos[1], -modOut$pred$probCurves[[id]]*4 + pos[2] + 2, col = cols[2], lwd = 3, xpd = NA)
  arrows(modOut$coor$xcoor[id], modOut$coor$ycoor[id], pos[1]-2.5, pos[2], col = "grey40", length = 0, lwd = 2, lty = 2, xpd = NA)
  points(stXcoor(modOut$pred$obsDat$respVal[id]) + pos[1], -(1-modOut$pred$obsDat$potentEn[id])*4 + pos[2] + 2, pch = 16, col = cols[1], cex = 1.3, xpd = NA)
}
colTipDirection <- function(tab, dat, cols = c("#000000","#FF0000", "#0000FF"), BR = FALSE) {
  tipResp <- tab$resp[tab$state < 0.5 & tab$catSt == 1]
  tipResp <- tipResp[!is.na(tipResp)]
  if (length(tipResp) < 1) NA else if (all(tipResp < dat))
    cols[2] else if (all(tipResp > dat)) cols[3] else if (BR) NA else cols[1]
}
boxCalc <- function(files){
  boxExt <- ext(files)
  boxExtPlus <- boxExt + c((boxExt[2] - boxExt[1]) * 0.02, (boxExt[2] - boxExt[1]) * 0.02, (boxExt[2] - boxExt[1]) * 0.02, (boxExt[2] - boxExt[1]) * 0.02)
  list(boxExt = boxExt, boxExtPlus = boxExtPlus)
}
plotMapBack <- function(files, borders, enlarge = FALSE, ann = TRUE, main = NULL, extVec = c(0,15,-5,5)){
  be <- boxCalc(files)
  addEnlarge <- if (enlarge) extVec else rep(0, 4)
  ylab <- if (ann) "Latitude" else ""
  plot(0,0, type = "n", xlim = be$boxExtPlus[1:2] + addEnlarge[1:2], ylim = be$boxExtPlus[3:4] + addEnlarge[3:4], asp = 1, xlab = "", ylab = ylab, axes = FALSE)
  auxUsr <- par("usr")
  auxUsr[2] <- be$boxExt[1] - auxUsr[1] + be$boxExt[2]
  axis(3, label = bquote(bold(.(main))), at = mean(auxUsr[1:2]), lwd.ticks = 0, line = -0.5)
  rect(auxUsr[1], auxUsr[3], auxUsr[2], auxUsr[4], col = "#81B5FF")
  bordersC <- crop(borders, auxUsr)
  polys(bordersC, col = "white", border = NA)
  axis(2, labels = c("", ""), at = c(-100,100))
  axis(1, labels = c("", ""), at = c(-100,auxUsr[2]), lwd.ticks = 0)
  axis(3, labels = c("", ""), at = c(-100,auxUsr[2]), lwd.ticks = 0)
  invisible(bordersC)
}
addBiomes <- function(biomes, boxExtPlus, biomeCols){
  biomesC <- crop(biomes, boxExtPlus)
  polys(biomesC, col = biomeCols[seq_along(biomesC$BIOME)])
}
plotTr <- function(modOut, files, borders, col, enlarge = FALSE, curves = NULL,
  dat = NULL, extVec = c(0,15,-5,5), main = NULL, ann = TRUE, biomes = NULL,
  biomeCols = c("#28A15A", "#72C935", "#CEDF88", "#88DFCE", "#FFF6D6", "#81B5FF"),
  onlyTipPoints = FALSE, addXlab = FALSE, colTip = NA){
  if (all(is.na(colTip))) colTip <- unlist(Map(colTipDirection, modOut$pred$tipStable, dat$datSelSt[, modOut$modPars[2]]))
  be <- boxCalc(files)
  bordersC <- plotMapBack(files, borders, enlarge, ann, main, extVec)
  auxUsr <- par("usr")
  auxUsr[2] <- be$boxExt[1] - auxUsr[1] + be$boxExt[2]
  if (!is.null(biomes)) addBiomes(biomes, be$boxExtPlus, biomeCols)
  if (!onlyTipPoints) points(modOut$coor, col = "grey", pch = 15, cex = 0.2)
  points(modOut$coor, col = colTip, pch = 15, cex = 0.4)
  polys(bordersC)
  rect(be$boxExtPlus[1], be$boxExtPlus[3], be$boxExtPlus[2], be$boxExtPlus[4], border = col, lwd = 2)
  if (addXlab) ann <- TRUE
  if (ann) axis(1, labels = "Longitude", at = mean(be$boxExt[1:2]), tick = FALSE, line = 2)
  if (!is.null(curves)){
    ySeq <- seq(mean(c(be$boxExt[4], auxUsr[4])), mean(c(be$boxExt[3], auxUsr[3])), length.out = 5)
    coors <- data.frame(t(data.frame(xcoor = c(rep(auxUsr[2] + auxUsr[2] - be$boxExt[2], 5), rep(auxUsr[2] + (auxUsr[2] - be$boxExt[2])*2.5, 4)),
                                     ycoor = c(ySeq, ySeq[-1] - diff(ySeq)/2))))
    Map(addCurve, list(modOut), curves$ID, coors, curves$addCols)
  }
  invisible(colTip)
}
plotSc <- function(dat, modOut, prec = FALSE, eviSD = FALSE, xlab = "", ylab = "", ax2 = FALSE, xticks = "", supp = ""){
  colTip <- unlist(Map(colTipDirection, modOut$pred$tipStable, dat$datSelSt[, modOut$modPars[2]]))
  xvar <- if (prec) dat$datSel$precNew/100 else dat$datSel$tempNew/10 - 273.15
  yvar <- if (eviSD) dat$datSel$eviSDNew/10000 else dat$datSel$eviMNew/10000
  if (supp != ""){
  xvar <- if (prec) dat$datSel[[paste0("precSupp", supp)]]/100 else dat$datSel[[paste0("tempSupp", supp)]]/10 - 273.15
  yvar <- if (eviSD) dat$datSel$eviSDSuppNew/10000 else dat$datSel$eviMSuppNew/10000
  }
  plot(xvar, yvar, col = "grey", axes = FALSE, cex = 0.5, xlab = xlab, ylab = ylab)
  points(xvar, yvar, col = colTip, pch = 1, cex = 0.5)
  axis(1, labels = xticks, at = xticks)
  if (ax2) axis(2, las = 2)
  if (substitute(dat) == "datRU") box(col = trcol[1], lwd = 3) else box(col = trcol[2], lwd = 3)
}
plotMap <- function(coors, cols, files, borders, borderCol, biomes = NULL, biomeCols = NULL){
  bordersC <- plotMapBack(files, borders, ann = FALSE)
  be <- boxCalc(files)
  auxUsr <- par("usr")
  auxUsr[2] <- be$boxExt[1] - auxUsr[1] + be$boxExt[2]
  if (!is.null(biomes)) addBiomes(biomes, be$boxExtPlus, biomeCols)
  points(coors, col = cols, pch = 15, cex = 0.4)
  polys(bordersC)
  rect(be$boxExtPlus[1], be$boxExtPlus[3], be$boxExtPlus[2], be$boxExtPlus[4], border = borderCol, lwd = 2)
}

# biomeEco <- vect("data/wwf_biomes/6kcchn7e3u_official_teow/official/wwf_terr_ecos.shp")
# biomes <- aggregate(biomeEco, by = "BIOME")
# writeVector(biomes, "data/wwf_biomes/biomes.gpkg")
biomes <- vect("data/wwf_biomes/biomes.gpkg")
selBiomesRU <- crop(biomes[biomes$BIOME %in% c(1,2,4,6,7,8,9,11,13,98)], c(67.9, 80.1, 44.9, 74.1))
selBiomesSA <- biomes[biomes$BIOME %in% c(1,2,4,6,7,9,11,13,98)]
cols <- c("#28A15A", "#72C935", "#CEDF88", "#88DFCE", "#FFF6D6", "#81B5FF")
countriesRU <- geodata::gadm(c("RUS", "KAZ", "CHN", "KGZ", "UZB", "TJK", "TKM", "AFG", "MNG", "PAK"), path = "data/gadm/", level = 0, resolution = 2)
countriesSA <- geodata::gadm(c("BRA", "ARG", "PRY", "BOL", "SUR", "GUY", "VEN", "URY", "GUF"), path = "data/gadm/", level = 0, resolution = 2)
trcol <- c("orange3", "darkgreen")

# Figure 2 - boreal transect
# png("figures/Fig2_boreal.png", height = 4800*2, width = 480*16, res = 72*16)
par(mai = c(0.8,0.55,0.5,0.01))
layout(matrix(c(1,3,1,3,1,4,2,4,2,5,2,5,6,7,6,7), nrow = 2))
plotTr(modOuts[[which(modsParsS == "RUeviMNewNew")]], filesRU, countriesRU,
       trcol[1], enlarge = TRUE, dat = datRU, main = "Growing season EVI", extVec = c(0,9,-5,5),
       curves = list(ID = c(591,551,670,1335,18,809,9,14,150),
                     addCols = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)),
       biomes = selBiomesRU, onlyTipPoints = TRUE)
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
par(mai = c(0.8,0.35,0.5,0.21))
plotTr(modOuts[[which(modsParsS == "RUeviSDNewNew")]], filesRU, countriesRU,
       trcol[1], enlarge = TRUE, dat = datRU, main = "EVI variability", extVec = c(0,9,-5,5),
       ann = FALSE, addXlab = TRUE, curves = list(ID = c(591,1542,22,173,87,809,48,14,278),
                     addCols = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)),
       biomes = selBiomesRU, onlyTipPoints = TRUE)
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
par(mai = c(0.8,0.55,0.5,0.01))
plotTr(modOuts[[which(modsParsS == "RUeviMNewNew")]], filesRU, countriesRU,
       trcol[1], dat = datRU, main = "Dist. to tipping point",
       colTip = modOuts[[which(modsParsS == "RUeviMNewNew")]]$colsDistT)
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotTr(modOuts[[which(modsParsS == "RUeviMNewNew")]], filesRU, countriesRU,
       trcol[1], dat = datRU, main = "Dist. to stable state", ann = FALSE, addXlab = TRUE,
       colTip = modOuts[[which(modsParsS == "RUeviMNewNew")]]$colsDistS)
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
plotTr(modOuts[[which(modsParsS == "RUeviMNewNew")]], filesRU, countriesRU,
       trcol[1], dat = datRU, main = "Potential energy", ann = FALSE, addXlab = TRUE,
       colTip = modOuts[[which(modsParsS == "RUeviMNewNew")]]$colsPotE)
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
par(mai = c(0,0,0,0))
plot(c(-180, 195), c(-650,163), type = "n", axes = FALSE, ann = FALSE)
rect(-180,-90,195,90,col = "#81B5FF")
maps::map("world", fill=TRUE, col="white", add = TRUE, border = NA)
rect(ext(filesRU)[1], ext(filesRU)[3], ext(filesRU)[2], ext(filesRU)[4], col = trcol[1], border = NA)
#rect(ext(filesSA)[1], ext(filesSA)[3], ext(filesSA)[2], ext(filesSA)[4], col = trcol[2])
legend(-180,-150, pch = 15, col = c("#FF0000", "#0000FF"), legend = c("Higher state variable", "Lower state variable"), title = "Alternative stable states", bty = "n", cex = 0.8, title.adj = 0)
bOrd <- c(4,2,1,3,5,6)
legend(-180,-300, fill = cols[bOrd], legend = c("Temp. broadleaf and mixed forests", "Boreal forests/taiga", "Temp. grass., savannas, shrub.", "Tundra", "Deserts, xeric shrublands", "Water")[bOrd], bty = "n", title = "Biomes", cex = 0.8, title.adj = 0)

plot(0:1, 0:1, type = "n", axes = FALSE, ann = FALSE)
xSq <- seq(0.1, 0.9, length.out = 101)
rect(head(xSq, -1),0.47,tail(xSq, -1),0.53, col = heat.colors(130)[1:100], border = NA)
rect(0.1,0.47,0.9,0.53, lwd = 1)
text(0.5,0.58,"Resilience")
text(c(0.2,0.8),0.42,c("Low", "High"))

par(new = TRUE, fig = c(0,1,0,1), mai = c(0,0,0.2,0))
plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "", cex.main = 1.5)
text(0.34, 1.055, "Stable states", cex = 1.5, font = 2, xpd = NA)
text(0.02,1.02,"A", cex = 2)
text(0.39,1.02,"B", cex = 2)
par(new = TRUE, fig = c(0,1,0,0.5))
plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "", cex.main = 1.5)
text(0.40, 1.06, "Ecological resilience", cex = 1.5, font = 2, xpd = NA)
text(0.02,0.995,"C", cex = 2)
text(0.29,0.995,"D", cex = 2)
text(0.555,0.995,"E", cex = 2)

# Figure 3 - tropical transect
# png("figures/Fig3_tropical.png", height = 4800*2, width = 480*16, res = 72*16)
par(mai = c(0.8,0.55,0.5,0.01))
layout(matrix(c(1,3,1,3,1,4,2,4,2,5,2,5,6,7,6,7), nrow = 2))
plotTr(modOuts[[which(modsParsS == "SAeviMNewNew")]], filesSA, countriesSA,
       trcol[2], enlarge = TRUE, dat = datSA, main = "Growing season EVI", extVec = c(0,9,-5,5),
       curves = list(ID = c(210,600,166,216,2938,164,34,388,361),
                     addCols = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)),
       biomes = selBiomesSA, onlyTipPoints = TRUE)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
axis(2, labels = paste0(-3:0*10, "°"), at = -3:0*10, las = 2)
par(mai = c(0.8,0.35,0.5,0.21))
plotTr(modOuts[[which(modsParsS == "SAeviSDNewNew")]], filesSA, countriesSA,
       trcol[2], enlarge = TRUE, dat = datSA, main = "EVI variability", extVec = c(0,9,-5,5),
       ann = FALSE, addXlab = TRUE, curves = list(ID = c(210,600,166,90,880,164,34,388,2117),
                     addCols = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)),
       biomes = selBiomesSA, onlyTipPoints = TRUE)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
par(mai = c(0.8,0.55,0.5,0.01))
plotTr(modOuts[[which(modsParsS == "SAeviMNewNew")]], filesSA, countriesSA,
       trcol[2], dat = datRU, main = "Dist. to tipping point",
       colTip = modOuts[[which(modsParsS == "SAeviMNewNew")]]$colsDistT)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
axis(2, labels = paste0(-3:0*10, "°"), at = -3:0*10, las = 2)
plotTr(modOuts[[which(modsParsS == "SAeviMNewNew")]], filesSA, countriesSA,
       trcol[2], dat = datRU, main = "Dist. to stable state", ann = FALSE, addXlab = TRUE,
       colTip = modOuts[[which(modsParsS == "SAeviMNewNew")]]$colsDistS)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
plotTr(modOuts[[which(modsParsS == "SAeviMNewNew")]], filesSA, countriesSA,
       trcol[2], dat = datRU, main = "Potential energy", ann = FALSE, addXlab = TRUE,
       colTip = modOuts[[which(modsParsS == "SAeviMNewNew")]]$colsPotE)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
par(mai = c(0,0,0,0))
plot(c(-180, 195), c(-650,163), type = "n", axes = FALSE, ann = FALSE)
rect(-180,-90,195,90,col = "#81B5FF")
maps::map("world", fill=TRUE, col="white", add = TRUE, border = NA)
#rect(ext(filesRU)[1], ext(filesRU)[3], ext(filesRU)[2], ext(filesRU)[4], col = trcol[1], border = NA)
rect(ext(filesSA)[1], ext(filesSA)[3], ext(filesSA)[2], ext(filesSA)[4], col = trcol[2], border = NA)
legend(-180,-150, pch = 15, col = c("#FF0000", "#0000FF"), legend = c("Higher state variable", "Lower state variable"), title = "Alternative stable states", bty = "n", cex = 0.8, title.adj = 0)
legend(-180,-300, fill = cols[1:4], legend = c("(Sub)tropical moist broadl. forests", "(Sub)tropical dry broadleaf forests", "(Sub)tropical grasslands, savannas", "Flooded grasslands, savannas"), bty = "n", title = "Biomes", cex = 0.8, title.adj = 0)
plot(0:1, 0:1, type = "n", axes = FALSE, ann = FALSE)
xSq <- seq(0.1, 0.9, length.out = 101)
rect(head(xSq, -1),0.47,tail(xSq, -1),0.53, col = heat.colors(130)[1:100], border = NA)
rect(0.1,0.47,0.9,0.53, lwd = 1)
text(0.5,0.58,"Resilience")
text(c(0.2,0.8),0.42,c("Low", "High"))

par(new = TRUE, fig = c(0,1,0,1), mai = c(0,0,0.2,0))
plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "", cex.main = 1.5)
text(0.34, 1.055, "Stable states", cex = 1.5, font = 2, xpd = NA)
text(0.02,1.02,"A", cex = 2)
text(0.39,1.02,"B", cex = 2)
par(new = TRUE, fig = c(0,1,0,0.5))
plot(0:1,0:1, type = "n", axes = FALSE, xlab = "", ylab = "", cex.main = 1.5)
text(0.40, 1.06, "Ecological resilience", cex = 1.5, font = 2, xpd = NA)
text(0.02,0.995,"C", cex = 2)
text(0.29,0.995,"D", cex = 2)
text(0.555,0.995,"E", cex = 2)

# Figure 4 - vegetation height
yLabPos <- function(){
  p <- par("usr")
  p[3] - (p[4]-p[3])/12
}
boxplotCHM <- function(modOut, dat, tit, yaxis = TRUE, ylim = range(dat$datSel$chm[!is.na(cTip)])){
  cTip <- factor(unlist(Map(colTipDirection, modOut$pred$tipStable, dat$datSelSt[, modOut$modPars[2]], BR = TRUE)), levels = c("#FF0000", "#0000FF"))
  ylab <- if (yaxis) "Vegetation height (m)" else ""
  boxplot(dat$datSel$chm ~ cTip, col = c("red", "blue"), xlab = "", ylab = ylab, names = rep("", 2), axes = FALSE, ylim = ylim)
  box()
  axis(1, labels = c("", ""), at = 1:2)
  if (yaxis) axis(2, las = 2)
  axis(3, at = mean(par("usr")[1:2]), label = tit, tick = FALSE, font = 2, cex.axis = 1.2)
  points(c(1, 2), rep(yLabPos(), 2), bg = "grey", pch = 21, col = c("red", "blue"), xpd = TRUE, cex = 2, lwd = 2)
  b <- dat$datSel$chm[cTip == "#0000FF"]
  r <- dat$datSel$chm[cTip == "#FF0000"]
  p <- par("usr")
  if (t.test(b, r)$p.value < 0.001) text(1.5, p[3] + (p[4]-p[3]) * 0.95, "***", font = 2, cex = 2)
  invisible(ylim)
}

# png("figures/Fig4_chm.png", height = 4800, width = 4800, res = 720)
par(mfrow = c(2,2), mai = c(0.5,0.7,0.8,0.1))
boxplotCHM(modOuts[[which(modsParsS == "RUeviMNewNew")]], datRU, tit = "Growing season EVI")
boxplotCHM(modOuts[[which(modsParsS == "RUeviSDNewNew")]], datRU, tit = "EVI variability")
legend(1, -6, pch = 21, col = c("red", "blue"), pt.bg = "grey", pt.cex = 2, pt.lwd = 2, xpd = NA, legend = c("Higher state variable", "Lower state variable"), title = "Alternative stable states")
par(mai = c(0.3,0.7,1.0,0.1))
boxplotCHM(modOuts[[which(modsParsS == "SAeviMNewNew")]], datSA, tit = "Growing season EVI")
boxplotCHM(modOuts[[which(modsParsS == "SAeviSDNewNew")]], datSA, tit = "EVI variability")
par(new = TRUE, fig = c(0,1,0,1), mai = c(0,0,0.2,0))
plot(0:1,0:1, type = "n", axes = FALSE, main = "Boreal transect", xlab = "", ylab = "", cex.main = 1.5)
text(0.02,0.97,"A", cex = 2)
text(0.56,0.97,"B", cex = 2)
par(new = TRUE, fig = c(0,1,0,0.47))
plot(0:1,0:1, type = "n", axes = FALSE, main = "Tropical transect", xlab = "", ylab = "", cex.main = 1.5)
text(0.02,0.89,"C", cex = 2)
text(0.56,0.89,"D", cex = 2)

# Info
altArea <- setNames(lapply(modOuts, function(x) !is.na(x$colsDistT)), modsParsS)
sum(altArea$RUeviMNewNew | altArea$RUeviSDNewNew) / length(altArea[[1]])
sum(altArea$SAeviMNewNew | altArea$SAeviSDNewNew) / length(altArea[[1]])
sum(altArea$SAeviMSuppNewSuppNew | altArea$SAeviSDSuppNewSuppNew) / length(altArea[[1]])
sum(altArea$SAeviMNewNew | altArea$SAeviSDNewNew | altArea$SAeviMSuppNewSuppNew | altArea$SAeviSDSuppNewSuppNew) / length(altArea[[1]])

getChm <- function(modName, dat){
  mod <- modOuts[[which(modsParsS == modName)]]
  colTip <- unlist(Map(colTipDirection, mod$pred$tipStable, dat$datSelSt[, mod$modPars[2]]))
  tapply(dat$datSel$chm, colTip, mean)
  #cbind(dat$datSel$chm, colTip)
}
getChm("RUeviMNewNew", datRU)
getChm("SAeviMNewNew", datSA)
getChm("SAeviMSuppNewSuppNew", datSA)

out <- setNames(lapply(modOuts, function(x) round(sum(!is.na(x$colsDistT))/length(x$colsDistT), 3)), modsParsS)
out$RUeviMNewNew
out$RUeviSDNewNew
out$SAeviMNewNew
out$SAeviSDNewNew
out$SAeviMSuppNewSuppNew
out$SAeviSDSuppNewSuppNew
# change in precipitation
round(sum((datSA$datSel$precNew - datSA$datSel$precOld) > 0) / length(datSA$datSel$precNew), 3)

  # area with alternative stable state
out <- setNames(lapply(modOuts, function(x) round(sum(!is.na(x$colsDistT))/length(x$colsDistT), 3)), modsParsS)
range(out[1:8]) # RU
range(out[9:16]) # SA
  # change in distance to stable state
calcChange <- function(old, new, ttest = FALSE, ...) {
  out <- c(old = mean(old), new = mean(new), change = mean(new) / mean(old) - 1)
  if (ttest) out <- list(t.test(old, new, ...), out)
  out
}
evalChangeDistToState <- function(modOut1, modOut2, dat){
  dist1 <- modOut1$distToState * 0.0001
  dist2 <- modOut2$distToState * 0.0001
  list(t.test(abs(dist1), abs(dist2), paired = TRUE), calcChange(abs(dist1), abs(dist2)))
}
evalChangeDistToState(modOuts[[which(modsParsS == "RUeviMOldOld")]], modOuts[[which(modsParsS == "RUeviMNewNew")]], datRU)
evalChangeDistToState(modOuts[[which(modsParsS == "RUeviSDOldOld")]], modOuts[[which(modsParsS == "RUeviSDNewNew")]], datRU)
evalChangeDistToState(modOuts[[which(modsParsS == "SAeviMOldOld")]], modOuts[[which(modsParsS == "SAeviMNewNew")]], datSA)
evalChangeDistToState(modOuts[[which(modsParsS == "SAeviSDOldOld")]], modOuts[[which(modsParsS == "SAeviSDNewNew")]], datSA)
  # change in ecosystem state variables
calcChange(datRU$datSel$eviMOld/10000, datRU$datSel$eviMNew/10000, ttest = TRUE, paired = TRUE)
calcChange(datRU$datSel$eviSDOld/10000, datRU$datSel$eviSDNew/10000, ttest = TRUE, paired = TRUE)
calcChange(datSA$datSel$eviMOld/10000, datSA$datSel$eviMNew/10000, ttest = TRUE, paired = TRUE)
calcChange(datSA$datSel$eviSDOld/10000, datSA$datSel$eviSDNew/10000, ttest = TRUE, paired = TRUE)
  # change in temperature
round(sum((datRU$datSel$tempNew - datRU$datSel$tempOld) / 10 > 0) / length(datRU$datSel$tempNew), 3)
round(mean((datSA$datSel$tempNew - datSA$datSel$tempOld) / 10), 2)

# supplementary figures
# Figure S1
# png("figures/FigS1_boreal.png", height = 4800*2, width = 480*18, res = 72*14)
par(mfrow = c(2,4), mai = c(0.3,0.3,0.5,0.02))
plotTr(modOuts[[which(modsParsS == "RUeviMNewNew")]], filesRU, countriesRU,
       trcol[1], main = "2017-2022, short-term climate", ann = FALSE, dat = datRU)
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotTr(modOuts[[which(modsParsS == "RUeviMNewFull")]], filesRU, countriesRU,
       trcol[1], main = "2017-2022, long-term climate", ann = FALSE, dat = datRU)
plotTr(modOuts[[which(modsParsS == "RUeviMOldOld")]], filesRU, countriesRU,
       trcol[1], main = "2000-2005, short-term climate", ann = FALSE, dat = datRU)
plotTr(modOuts[[which(modsParsS == "RUeviMOldFull")]], filesRU, countriesRU,
       trcol[1], main = "2000-2005, long-term climate", ann = FALSE, dat = datRU)
plotTr(modOuts[[which(modsParsS == "RUeviSDNewNew")]], filesRU, countriesRU,
       trcol[1], main = "2017-2022, short-term climate", ann = FALSE, dat = datRU)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotTr(modOuts[[which(modsParsS == "RUeviSDNewFull")]], filesRU, countriesRU,
       trcol[1], main = "2017-2022, long-term climate", ann = FALSE, dat = datRU)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
plotTr(modOuts[[which(modsParsS == "RUeviSDOldOld")]], filesRU, countriesRU,
       trcol[1], main = "2000-2005, short-term climate", ann = FALSE, dat = datRU)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
plotTr(modOuts[[which(modsParsS == "RUeviSDOldFull")]], filesRU, countriesRU,
       trcol[1], main = "2000-2005, long-term climate", ann = FALSE, dat = datRU)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
par(new = TRUE, fig = c(0,1,0,1), mai = c(0,0,0.2,0))
plot(0,0, type = "n", axes = FALSE, main = "Growing season EVI", xlab = "", ylab = "", cex.main = 1.5)
par(new = TRUE, fig = c(0,1,0,0.5))
plot(0,0, type = "n", axes = FALSE, main = "EVI variability", xlab = "", ylab = "", cex.main = 1.5)

# Figure S2
# png("figures/FigS2_predictions.png", height = 4800, width = 4800, res = 720)
par(mfrow = c(2,4), mai = c(0.3,0.1,0.4,0.55))
modPotE <- calcEffPotE(modOuts[[which(modsParsS == "RUeviMOldOld")]], ylab = "Change in Growing season EVI", cf = 0.0001, pos = c(-1, -0.2, 0.03, -9.6, -0.2, 0.03), axis1 = c(-0.25,0,0.25), sepPlots = TRUE, xlab = FALSE)
modPotE <- calcEffPotE(modOuts[[which(modsParsS == "SAeviMOldOld")]], ylab = "Change in Growing season EVI", cf = 0.0001, pos = c(-1, -0.25, 0.035, -2.9, -0.25, 0.035), sepPlots = TRUE, xlab = FALSE)
par(mai = c(0.55,0.1,0.15,0.55))
modPotE <- calcEffPotE(modOuts[[which(modsParsS == "RUeviSDOldOld")]], ylab = "Change in EVI variability", cf = 0.0001, pos = c(-1.1, -0.05, 0.009, -8.3, 0.058, 0.009), axis1 = -2:2/10, sepPlots = TRUE)
modPotE <- calcEffPotE(modOuts[[which(modsParsS == "SAeviSDOldOld")]], ylab = "Change in EVI variability", cf = 0.0001, pos = c(-1, -0.13, 0.025, -9.0, -0.13, 0.025), axis1 = -2:2/10, sepPlots = TRUE)
par(new = TRUE, fig = c(0,1,0,1), mai = c(0,0,0,0), mfrow = c(1,1))
plot(0,0, type = "n", axes = FALSE, main = "", xlab = "", ylab = "")
text(-0.55, 0.99, "Boreal zone", cex = 1, font = 2)
text(0.55, 0.99, "Tropical zone", cex = 1, font = 2)

# Figure S3
# png("figures/FigS3_tropical.png", height = 4800*2, width = 480*18, res = 72*14)
par(mfrow = c(2,4), mai = c(0.3,0.3,0.5,0.02))
plotTr(modOuts[[which(modsParsS == "SAeviMNewNew")]], filesSA, countriesSA,
       trcol[2], main = "2017-2022, short-term climate", ann = FALSE, dat = datSA)
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotTr(modOuts[[which(modsParsS == "SAeviMNewFull")]], filesSA, countriesSA,
       trcol[2], main = "2017-2022, long-term climate", ann = FALSE, dat = datSA)
plotTr(modOuts[[which(modsParsS == "SAeviMOldOld")]], filesSA, countriesSA,
       trcol[2], main = "2000-2005, short-term climate", ann = FALSE, dat = datSA)
plotTr(modOuts[[which(modsParsS == "SAeviMOldFull")]], filesSA, countriesSA,
       trcol[2], main = "2000-2005, long-term climate", ann = FALSE, dat = datSA)
plotTr(modOuts[[which(modsParsS == "SAeviSDNewNew")]], filesSA, countriesSA,
       trcol[2], main = "2017-2022, short-term climate", ann = FALSE, dat = datSA)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotTr(modOuts[[which(modsParsS == "SAeviSDNewFull")]], filesSA, countriesSA,
       trcol[2], main = "2017-2022, long-term climate", ann = FALSE, dat = datSA)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
plotTr(modOuts[[which(modsParsS == "SAeviSDOldOld")]], filesSA, countriesSA,
       trcol[2], main = "2000-2005, short-term climate", ann = FALSE, dat = datSA)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
plotTr(modOuts[[which(modsParsS == "SAeviSDOldFull")]], filesSA, countriesSA,
       trcol[2], main = "2000-2005, long-term climate", ann = FALSE, dat = datSA)
axis(1, labels = paste0(c(70,75,80), "°"), at = c(70,75,80))
par(new = TRUE, fig = c(0,1,0,1), mai = c(0,0,0.2,0))
plot(0,0, type = "n", axes = FALSE, main = "Growing season EVI", xlab = "", ylab = "", cex.main = 1.5)
par(new = TRUE, fig = c(0,1,0,0.5))
plot(0,0, type = "n", axes = FALSE, main = "EVI variability", xlab = "", ylab = "", cex.main = 1.5)

# Figure S4 - wet season
# png("figures/FigS4_wetSeason.png", height = 480*20, width = 480*15, res = 72*17)
par(mfrow = c(2,2), mai = c(0.4,0.7,0.5,0.1))
plotTr(modOuts[[which(modsParsS == "SAeviMSuppNewSuppNew")]], filesSA, countriesSA,
       trcol[2], main = "Growing season EVI", ann = FALSE, dat = datSA)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
axis(2, labels = paste0(-3:0*10, "°"), at = -3:0*10, las = 2)
plotTr(modOuts[[which(modsParsS == "SAeviSDSuppNewSuppNew")]], filesSA, countriesSA,
       trcol[2], main = "EVI variability", ann = FALSE, dat = datSA)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
ylimAux <- boxplotCHM(modOuts[[which(modsParsS == "SAeviMSuppNewSuppNew")]], datSA, tit = "")
boxplotCHM(modOuts[[which(modsParsS == "SAeviSDSuppNewSuppNew")]], datSA, tit = "", yaxis = FALSE, ylim = ylimAux)

# Figure S5
chPTRU <- ((datRU$datSel$precNew - datRU$datSel$precOld)/100 > 0) + 1 + ((datRU$datSel$tempNew - datRU$datSel$tempOld)/10 > 0) * 2
chPTSA <- ((datSA$datSel$precNew - datSA$datSel$precOld)/100 > 0) + 1 + ((datSA$datSel$tempNew - datSA$datSel$tempOld)/10 > 0) * 2
chEVIRU <- ((datRU$datSel$eviMNew - datRU$datSel$eviMOld)/100 > 0) + 1 + ((datRU$datSel$eviSDNew - datRU$datSel$eviSDOld)/10 > 0) * 2
chEVISA <- ((datSA$datSel$eviMNew - datSA$datSel$eviMOld)/100 > 0) + 1 + ((datSA$datSel$eviSDNew - datSA$datSel$eviSDOld)/10 > 0) * 2
getCoor <- function(transect){
  sampID <- get(paste0("dat", transect))$sampID
  files <- get(paste0("files", transect))
  data.frame(xcoor = xFromCol(files)[(sampID %% dim(files)[2]) + 1], ycoor = yFromRow(files)[1 + trunc(sampID/dim(files)[2])])
}
cols2 <- c("brown", "lightblue", "red", "purple")
cols3 <- c("brown", "#CEDF88", "blue", "#28A15A")

# png("figures/FigS5_change.png", height = 4800*2, width = 4800, res = 72*14)
par(mfrow = c(2, 2), mai = c(1.4,0.6,0.1,0.1))
plotMap(getCoor("RU"), cols2[chPTRU], filesRU, countriesRU, trcol[1])
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotMap(getCoor("SA"), cols2[chPTSA], filesSA, countriesSA, trcol[2])
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
axis(2, labels = paste0(-3:0*10, "°"), at = -3:0*10, las = 2)

plotMap(getCoor("RU"), cols3[chEVIRU], filesRU, countriesRU, trcol[1])
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotMap(getCoor("SA"), cols3[chEVISA], filesSA, countriesSA, trcol[2])
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
axis(2, labels = paste0(-3:0*10, "°"), at = -3:0*10, las = 2)

ys <- 5.15
par(mai = c(0,0,0,0), new = TRUE, mfrow = c(1, 1))
plot(c(0, 1), c(0.45,10), type = "n", ann = FALSE, axes = FALSE)
rect(0.5,0.2+ys,0.6,0.7+ys)
rect(0.5,0.2+ys,0.55,0.45+ys, col = cols2[3])
rect(0.55,0.2+ys,0.6,0.45+ys, col = cols2[4])
rect(0.5,0.45+ys,0.55,0.7+ys, col = cols2[1])
rect(0.55,0.45+ys,0.6,0.7+ys, col = cols2[2])
text(0.45,0.45+ys, "Temperature", adj = c(1, 0.5))
text(0.48,0.575+ys, "-", cex = 1.5)
text(0.48,0.325+ys, "+", cex = 1.5)
text(0.55,1.05+ys, "Precipitation")
text(0.525,0.85+ys, "-", cex = 1.5)
text(0.575,0.85+ys, "+", cex = 1.5)

par(mai = c(0,0,0,0), new = TRUE, mfrow = c(1, 1))
plot(c(0, 1), c(0.45,10), type = "n", ann = FALSE, axes = FALSE)
rect(0.5,0.2,0.6,0.7)
rect(0.5,0.2,0.55,0.45, col = cols3[3])
rect(0.55,0.2,0.6,0.45, col = cols3[4])
rect(0.5,0.45,0.55,0.7, col = cols3[1])
rect(0.55,0.45,0.6,0.7, col = cols3[2])
text(0.45,0.45, "EVI variability", adj = c(1,0.5))
text(0.48,0.575, "-", cex = 1.5)
text(0.48,0.325, "+", cex = 1.5)
text(0.55,1.05, "Growing season EVI")
text(0.525,0.85, "-", cex = 1.5)
text(0.575,0.85, "+", cex = 1.5)

# Figure S6
# png("figures/FigS6_biomes.png", height = 4800*2, width = 4800, res = 72*14)
layout(matrix(c(1,1,1,3,3,3,5,2,2,2,4,4,4,6), ncol = 2))
par(mai = c(0.1,0.3,0.3,0.1))
plotTr(modOuts[[which(modsParsS == "RUeviMNewNew")]], filesRU, countriesRU,
       trcol[1], dat = datRU,
       biomes = selBiomesRU, onlyTipPoints = TRUE, ann = FALSE)
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotTr(modOuts[[which(modsParsS == "SAeviMNewNew")]], filesSA, countriesSA,
       trcol[2], dat = datSA,
       biomes = selBiomesSA, onlyTipPoints = TRUE, ann = FALSE)
axis(2, labels = paste0(-3:0*10, "°"), at = -3:0*10, las = 2)
plotTr(modOuts[[which(modsParsS == "RUeviSDNewNew")]], filesRU, countriesRU,
       trcol[1], dat = datRU,
       biomes = selBiomesRU, onlyTipPoints = TRUE, ann = FALSE)
axis(1, labels = paste0(c(65,70,75,80), "°"), at = c(65,70,75,80))
axis(2, labels = paste0(4:8*10, "°"), at = 4:8*10, las = 2)
plotTr(modOuts[[which(modsParsS == "SAeviSDNewNew")]], filesSA, countriesSA,
       trcol[2], dat = datSA,
       biomes = selBiomesSA, onlyTipPoints = TRUE, ann = FALSE)
axis(1, labels = paste0(1:4*5-70, "°"), at = 1:4*5-70)
axis(2, labels = paste0(-3:0*10, "°"), at = -3:0*10, las = 2)
par(mai = c(0.0,0.3,0.3,0.1))
plot(0,0, type = "n", ann = FALSE, axes = FALSE)
bOrd <- c(4,2,1,3,5,6)
legend("topleft", fill = cols[bOrd], legend = c("Temperate broadleaf and mixed forests", "Boreal forests/taiga", "Temp. grasslands, savannas, shrublands", "Tundra", "Deserts, xeric shrublands", "Water")[bOrd], bty = "n", title = "Boreal transect", cex = 0.8)
plot(0,0, type = "n", ann = FALSE, axes = FALSE)
legend("topleft", fill = cols[1:4], legend = c("(Sub)tropical moist broadleaf forests", "(Sub)tropical dry broadleaf forests", "(Sub)tropical grasslands, savannas, shrub.", "Flooded grasslands, savannas"), bty = "n", title = "Tropical transect", cex = 0.8)
par(new = TRUE, fig = c(0,1,0,1), mai = c(0,0,0.2,0))
plot(0,0, type = "n", axes = FALSE, main = "Growing season EVI", xlab = "", ylab = "")
par(new = TRUE, fig = c(0,1,0,0.57))
plot(0,0, type = "n", axes = FALSE, main = "EVI variability", xlab = "", ylab = "")

# Figure S7
# png("figures/FigS7_worldmap.png", height = 480*6, width = 480*13, res = 72*13)
maps::map("world", fill=TRUE, col="white", bg="#81B5FF", mar=c(0,0,0,0))
rect(ext(filesRU)[1], ext(filesRU)[3], ext(filesRU)[2], ext(filesRU)[4], col = trcol[1])
rect(ext(filesSA)[1], ext(filesSA)[3], ext(filesSA)[2], ext(filesSA)[4], col = trcol[2])

#_
