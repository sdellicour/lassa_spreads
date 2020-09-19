library(blockCV)
library(diagram)
library(dismo)
library(gbm)
library(geosphere)
library(lubridate)
library(ncdf4)
library(ncf)
library(pgirmess)
library(rgdal)
library(seraphim)
library(seqinr)
library(vioplot)

# A. SPECIES DISTRIBUTION MODEL ANALYSES
# B. PHYLOGEOGRAPHIC AND PHYLODYNAMIC ANALYSES
# B1. Preparation of LASV sequences
# B2. Preliminary phylogenetic analyses
# B3. Temporal signal analyses
# B4. Continuous phylogeographic analyses
# B5. Estimating dispersal statistics
# B6. Post hoc landscape phylogeography

wd = getwd()
wda1 = "A1_ENM_analyses_host_virus"
wdb1 = "B1_sequences_preparation"
wdb2 = "B2_preliminary_BEAST_runs"
wdb3 = "B3_temporal_signal_analyses"
wdb4 = "B4_RRW_phylogeography"
wdb5 = "B5_dispersal_stat_estimations"
wdb6 = "B6_all_seraphim_analyses"

study_areas = c("MRU","NGA")
segments = c("L","S"); analyses = c("2","3")
clades = c("MRU","NIG1","NIG2","NIG3")
e_NGA = extent(2, 15.5, 3.3, 14.4)
e_MRU = extent(-15, -4, 4, 13.04)

# A. SPECIES DISTRIBUTION MODEL ANALYSES

setwd(paste(wd,wda1,sep="/"))
analyses = c("M_natalensis","Lassa_virus")
mask = raster("Environmental_rasters/Mask.nc4")
continents = shapefile("Continents_shapefile/Continents.shp")

# A.1. Preparation of land cover and climatic environmental rasters

		# Source: http://luh.umd.edu/data.shtml; see also Lawrence et al. (2016, Geosci. Model Dev.)
		# Units: population (# of people), temperature (Kelvin, °C+273.15), precipitation (kg/m2/second)

	# A.1.1. Preparation of the African shapefile that will be used as a mask

africa1 = subset(continents, continents$CONTINENT=="Africa"); polygons = list(); c = 0
for (i in 1:length(africa1@polygons[[1]]@Polygons))
	{
		if (africa1@polygons[[1]]@Polygons[[i]]@area > 1)
			{
				c = c+1; polygons[[c]] = africa1@polygons[[1]]@Polygons[[i]]
			}
	}
pols = Polygons(polygons, 1); pols_list = list(); pols_list[[1]] = pols
africa2 = SpatialPolygons(pols_list); africa3 = gSimplify(africa2, 0.1)

	# A.1.2. Loading the human population, temperature, precipitation, and land cover rasters

population = raster("Environmental_rasters/Calibration/GSWP3_EWEMBI/population_histsoc_0p5deg_annual_1986_2005_timmean.nc4")
temperature = raster("Environmental_rasters/Calibration/GSWP3_EWEMBI/tas_day_GSWP3+EWEMBI_HistObs_r1i1p1_EWEMBI_1986_2005_timmean.nc4")
precipitation = raster("Environmental_rasters/Calibration/GSWP3_EWEMBI/pr_day_GSWP3+EWEMBI_HistObs_r1i1p1_EWEMBI_1986_2005_timmean.nc4")
land_cover = nc_open("Environmental_rasters/Calibration/GSWP3_EWEMBI/landcover_HistObs_annual_1986_2005_timmean.nc4") 
names(population) = "population"; names(temperature) = "temperature"; names(precipitation) = "precipitation"

	# A.1.3. Preparation of distinct land cover rasters from the initial ".nc" object

landCoverVariableIDs = names(land_cover$var); cols = list()
landCoverVariableNames = as.character(read.csv("Environmental_rasters/LC_vars.csv")[1:12,2])
land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
for (i in 2:13) land_covers1[[i-1]] = brick("Environmental_rasters/Calibration/GSWP3_EWEMBI/landcover_HistObs_annual_1986_2005_timmean.nc4", varname=landCoverVariableIDs[i])
for (i in 1:length(land_covers1))
	{
		names(land_covers1[[i]]) = landCoverVariableNames[i]
		cols[[i]] = colorRampPalette(brewer.pal(9,"YlGn"))(120)[11:110]
		if (i == 1)
			{
				r_global = land_covers1[[1]]
			}	else		{
				r_global[] = r_global[]+land_covers1[[i]][]	
			}
	}
if (!file.exists("Land_cover_rasters.pdf"))	
	{
		pdf("Land_cover_rasters_1.pdf", width=7.5, height=5.0); par(mfrow=c(4,3), oma=c(1.5,2.0,1,0.5), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
		for (i in 1:12)
			{
				plot(land_covers1[[i]], bty="n", box=F, axes=F, legend=F, col=c("gray90",cols[[i]]), colNA="white")
				plot(land_covers1[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.09,0.105,0.14,0.49), adj=3,
					 axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
				mtext(gsub("\\."," ",names(land_covers1[[i]])), side=1, adj=0.5, line=-1.8, at=40, cex=0.45, font=1, col="gray30")
			}
		dev.off()
	}

	# A.1.4. Saving distinct land cover rasters for the "seraphim" analyses

setwd(paste(wd,wdb6,sep="/"))
variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
				   "potentially forested secondary land","potentially non-forested secondary land")
for (i in 1:length(variable_names))
	{
		names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[i])
		if (length(indices) == 0) indices = which(grepl(variable_names[i],names))
		if (variable_names[i] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
		land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[i]
		if (length(indices) > 1)
			{
				for (j in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[j]]][]
			}
		if (!file.exists(paste0("Environmental_files/LUH2_",variable_codes[i],".asc")))
			{
				writeRaster(land_cover[[1]], paste0("Environmental_files/LUH2_",variable_codes[i],".asc"))
			}
		land_covers2[[i]] = land_cover[[1]]; land_covers3[[i]] = raster::aggregate(land_cover[[1]],2)
	}
setwd(paste(wd,wda1,sep="/"))

	# A.1.5. Selecting and treating the environmental rasters for the ENM analyses

envVariables = list()
envVariables[[1]] = temperature; envVariables[[2]] = precipitation
envVariables[[3]] = land_covers3[[4]] # primary forest areas
envVariables[[4]] = land_covers3[[5]] # primary non-forest areas
envVariables[[5]] = land_covers3[[6]] # secondary forest areas
envVariables[[6]] = land_covers3[[7]] # secondary non-forest areas
envVariables[[7]] = land_covers3[[1]] # croplands (all catergories)
envVariables[[8]] = land_covers3[[2]] # managed pasture + rangeland
pLog = population; pLog[] = log10(pLog[]+1); envVariables[[9]] = pLog
for (i in 1:length(envVariables)) envVariables[[i]][is.na(mask[])] = NA
for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], africa2, snap="out")
for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], africa2)
envVariables[[1]][] = envVariables[[1]][]-273.15 # conversion to Celcius degrees
envVariables[[2]][] = envVariables[[2]][]*60*60*24 # conversion to kg/m2/day
for (i in c(1,2,9)) envVariables[[i]][is.na(envVariables[[3]][])] = NA
rasters_stack = stack(envVariables)

	# A.1.6. Plotting the different environmental rasters used for the ENM analyses

showingPlots = FALSE; if (showingPlots == TRUE) {
dev.new(width=8, height=3); par(mfrow=c(2,7), oma=c(0,0,1.5,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
plot(envVariables[[1]], col=colorRampPalette(brewer.pal(9,"YlOrRd"))(150)[1:100], ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
mtext("Air temperature", side=3, line=0.3, cex=0.65, col="gray30"); mtext("near surface (°C)", side=3, line=-0.7, cex=0.65, col="gray30")
plot(envVariables[[1]], col=colorRampPalette(brewer.pal(9,"YlOrRd"))(150)[1:100], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T)
plot(envVariables[[2]], col=colorRampPalette(brewer.pal(9,"YlGnBu"))(100), ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
mtext("Precipitation", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0("(kg/m2/day)"), side=3, line=-0.7, cex=0.65, col="gray30")
plot(envVariables[[2]], col=colorRampPalette(brewer.pal(9,"YlGnBu"))(100), legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T)
plot(envVariables[[3]], col=colorRampPalette(c("gray97","chartreuse4"),bias=1)(100), ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
mtext("Forested", side=3, line=0.3, cex=0.65, col="gray30"); mtext("primary land", side=3, line=-0.7, cex=0.65, col="gray30")
plot(envVariables[[3]], col=colorRampPalette(c("gray97","chartreuse4"),bias=1)(100), legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.0,0), at=c(0,0.3,0.6,0.9)), alpha=1, side=3, horizontal=T)
plot(envVariables[[5]], col=colorRampPalette(c("gray97","olivedrab3"),bias=1)(100), ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
mtext("Forested", side=3, line=0.3, cex=0.65, col="gray30"); mtext("secondary land", side=3, line=-0.7, cex=0.65, col="gray30")
plot(envVariables[[5]], col=colorRampPalette(c("gray97","olivedrab3"),bias=1)(100), legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.0,0), at=c(0,0.3,0.6,0.9)), alpha=1, side=3, horizontal=T)
plot(envVariables[[7]], col=colorRampPalette(c("gray97","navajowhite4"),bias=1)(100), ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
mtext("Croplands", side=3, line=0.3, cex=0.65, col="gray30"); mtext("(all categories)", side=3, line=-0.7, cex=0.65, col="gray30")
plot(envVariables[[7]], col=colorRampPalette(c("gray97","navajowhite4"),bias=1)(100), legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T)
plot(envVariables[[8]], col=colorRampPalette(c("gray97","burlywood3"),bias=1)(100), ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
mtext("Pastures", side=3, line=0.3, cex=0.65, col="gray30"); mtext("and rangeland", side=3, line=-0.7, cex=0.65, col="gray30")
plot(envVariables[[8]], col=colorRampPalette(c("gray97","burlywood3"),bias=1)(100), legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.0,0), at=c(0,0.3,0.6,0.9)), alpha=1, side=3, horizontal=T)
plot(envVariables[[9]], col=colorRampPalette(brewer.pal(9,"BuPu"))(150)[1:100], ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
mtext("Human", side=3, line=0.3, cex=0.65, col="gray30"); mtext("population (log10)", side=3, line=-0.7, cex=0.65, col="gray30")
plot(envVariables[[9]], col=colorRampPalette(brewer.pal(9,"BuPu"))(150)[1:100], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T) }

# A.2. Defining the background areas for the host and for the virus

nullRaster = envVariables[[1]]; nullRaster[!is.na(nullRaster[])] = 1; names(nullRaster) = "nullRaster"
muridaeFam = read.csv("Occurrence_data_sets/Muridae_family_allData_RK220819.csv", header=T)[,c("longitude","latitude")]
natalensis = read.csv("Occurrence_data_sets/Mastomys_natalensis_RK260819.csv", header=T)[,c("longitude","latitude")]
backgroundCells1 = unique(raster::extract(nullRaster, muridaeFam, cellnumbers=T))
backgroundCells2 = unique(raster::extract(nullRaster, natalensis, cellnumbers=T))
background1 = nullRaster; background1[!(1:length(background1[]))%in%backgroundCells1] = NA
background2 = nullRaster; background2[!(1:length(background2[]))%in%backgroundCells2] = NA

# A.3. Boosted regression trees (BRT) analyses with spatial cross-validation

samplingPtsMinDist = function(observations, minDist=500, nberOfPoints=5)
	{
		# function written by Jean Artois (source: Dhingra, Artois, et al. 2016, eLife)
		indices = rep(NA, nberOfPoints)
		selection_list = list(1:nrow(observations)) 
  		indices[1] = sample(1:dim(observations)[1], 1)
		dists = list(spDistsN1(as.matrix(observations), as.matrix(observations[indices[1],]), longlat=T))
		for (i in 2:nberOfPoints)
			{
    				selection = which(dists[[(i-1)]] > minDist)
    				if (length(selection) == 0)
    					{
    						stop("Restarts the function with a smaller minimum distance")
					}
    				selection_list[[i]] = selection
    				test = table(unlist(selection_list))
    				indices_minDist = as.numeric(names(which(test==i)))
    				indices[i] = sample(indices_minDist, 1)   
				dists[[i]] = spDistsN1(as.matrix(observations), as.matrix(observations[indices[i],]), longlat=T)
			}
		return(indices)
	}
foldSelection = function(observations, selectedPoints)
	{
		# function written by Jean Artois (source: Dhingra, Artois, et al. 2016, eLife)
		fold_selection = sapply(1:nrow(observations), function(i) which.min(spDistsN1(as.matrix(selectedPoints), as.matrix(observations[i,]), longlat=T)))
		return(fold_selection)
	}

natalensis = read.csv("Occurrence_data_sets/Mastomys_natalensis_RK050820.csv", header=T)[,c("longitude","latitude")]
lassaVirus = read.csv("Occurrence_data_sets/Lassa_virus_cases_RK070820.csv", header=T)[,c("longitude","latitude")]
backgrounds = list(); backgrounds[[1]] = background1; backgrounds[[2]] = background2
datasets = list(); datasets[[1]] = natalensis; datasets[[2]] = lassaVirus

datas = list()
for (i in 1:length(datasets))
	{
		probas = values(backgrounds[[i]])[!is.na(values(backgrounds[[i]]))]; n = 1000
		if (n > sum(!is.na(values(backgrounds[[i]])))) n = sum(!is.na(values(backgrounds[[i]])))
		pseudo_absences = xyFromCell(backgrounds[[i]], sample(which(!is.na(values(backgrounds[[i]]))), n, prob=probas, replace=F))
		presences = cbind(datasets[[i]], rep(1,dim(datasets[[i]])[1]), raster::extract(rasters_stack, datasets[[i]]))
		absences = cbind(pseudo_absences, rep(0,dim(pseudo_absences)[1]), raster::extract(rasters_stack, pseudo_absences))
		colnames(absences)[1] = "longitude"; colnames(absences)[2] = "latitude"; colnames(absences)[3] = "response"
		colnames(presences)[3] = "response"; data = rbind(presences,absences); data_to_discard = c()
		for (j in 1:length(rasters_stack@layers))
			{
				data_to_discard = c(data_to_discard, which(is.na(raster::extract(rasters_stack[[j]],data[,1:2]))))
			}
		data_to_discard = data_to_discard[order(	data_to_discard)]
		data = data[which(!c(1:dim(data)[1])%in%data_to_discard),]
		cellIDs = cellFromXY(rasters_stack[[1]], data[,1:2]); buffer = c()
		for (j in 1:length(unique(cellIDs)))
			{	# Keeping only one presence or pseudo-absence point per raster cell (priority = presence points):
				if (sum(cellIDs==unique(cellIDs)[j]) > 1)
					{
						tmp = data[which(cellIDs==unique(cellIDs)[j]),]
						if (sum(tmp[,"response"]==1) == 0)
							{
								buffer = rbind(buffer, tmp[sample(1:dim(tmp)[1],1),])
							}
						if (sum(tmp[,"response"]==1) == 1)
							{
								buffer = rbind(buffer, tmp[which(tmp[,"response"]==1),])
							}
						if (sum(tmp[,"response"]==1) >= 2)
							{
								indices = which(tmp[,"response"]==1)
								buffer = rbind(buffer, tmp[sample(indices,1),])
							}
					}	else	{
						buffer = rbind(buffer, data[which(cellIDs==unique(cellIDs)[j]),])
					}
			}
		data = buffer; datas[[i]] = data
		plottingCorrelogram = FALSE
		if (plottingCorrelogram == TRUE)
			{
				correlogram = ncf::correlog(data[,"longitude"], data[,"latitude"], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)
				dev.new(width=4.5, height=3); par(mar=c(2.2,2.2,1.5,1.5))
				plot(correlogram$mean.of.class[-1], correlogram$correlation[-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.4,1.0), xlim=c(0,8500))
				abline(h=0, lwd=0.5, col="red", lty=2)
				points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.35, col="gray30")
				lines(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, col="gray30")
				axis(side=1, pos=-0.4, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,9000,1000))
				axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.4,1,0.2))
				title(xlab="distance (km2)", cex.lab=0.7, mgp=c(0.3,0,0), col.lab="gray30")
				title(ylab="correlation", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30")
			}
		theRanges = c(2000,2000)*1000 # distance in meters
		nberOfReplicates = 10 # one replicate = one folds partition
		gbm.x = names(rasters_stack)
		gbm.y = colnames(data)[3]
		offset = NULL
		tree.complexity = 5 # "tc" = number of nodes in the trees
		learning.rate = 0.005 # "lr" = contribution of each tree to the growing model
		bag.fraction = 0.80 # proportion of data used to train a given tree
		site.weights = rep(1, dim(data)[1])
		var.monotone = rep(0, length(gbm.x))
		n.folds = 5
		prev.stratify = TRUE
		family = "bernoulli"
		n.trees = 100 # initial number of trees
		step.size = 10 # interval at which the predictive deviance is computed and logged
					   # (at each interval, the folds are successively used as test data set
					   # and the remaining folds as training data sets to compute the deviance)
		max.trees = 10000 # maximum number of trees that will be considered
		tolerance.method = "auto"
		tolerance = 0.001
		plot.main = TRUE
		plot.folds = FALSE
		verbose = TRUE
		silent = FALSE
		keep.fold.models = FALSE
		keep.fold.vector = FALSE
		keep.fold.fit = FALSE
		showingFoldsPlot = FALSE		
		brt_model_ccvs = list() # classic cross-validations (CCVs)
		brt_model_scv1 = list() # spatial cross-validations 1 (SCV1)
		brt_model_scv2 = list() # spatial cross-validations 2 (SCV2)
		AUCs = matrix(nrow=nberOfReplicates, ncol=3)
		colnames(AUCs) = c("CCV_AUC","SCV1_AUC","SCV2_AUC")
		for (j in 1:nberOfReplicates)
			{
				# BRT with classic (standard) cross-validation (CCV):
				n.trees = 100; fold.vector = NULL
				brt_model_ccvs[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
					var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
					verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
				dev.copy2pdf(file=paste0("BRT_prediction_files/BRT_models/",analyses[i],"_CCV_replicate_",j,".pdf")); dev.off()
				AUCs[j,"CCV_AUC"] = brt_model_ccvs[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the CCV)
				object = brt_model_ccvs[[j]]; df = as.data.frame(rasters_stack)
				not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
				n.trees = brt_model_ccvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
				prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
				rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
				
				# BRT with spatial (geographic) cross-validation (SCV) based on the folds generation of Dhingra, Artois et al. (2016, eLife):
				folds_with_similar_sizes = FALSE; c = 0
				while (folds_with_similar_sizes == FALSE) # while loop to select a partition where the x folds gather at least
					{									  # a proportion = (1/(x+1)) of the total number of presence points
						data_presence = data[which(data[,3]==1),]; c = c+1; # print(c)
						fivePoints = samplingPtsMinDist(data_presence[,1:2], minDist=200, nberOfPoints=n.folds)
						fold.vector = foldSelection(data[,1:2], selectedPoints=data_presence[fivePoints,1:2])
						fold.vector_presences = fold.vector[which(data[,3]==1)]
						counts = hist(fold.vector_presences, plot=F)$counts
						props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
						if (min(props) > (1/(n.folds+1))) folds_with_similar_sizes = TRUE
					}
				if (showingFoldsPlot == TRUE)
					{
						par(mar=c(0,0,0,0), oma=c(0.0,3.6,0.0,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
						cols = c("olivedrab3","tan3","steelblue3","orange1","tomato2","mediumseagreen")[fold.vector]
						plot(backgrounds[[i]], col="gray90", useRaster=T, colNA=NA, box=F, axes=F, legend=F)
						pchs = c(16,3)[data[,3]+1]; cexs = c(0.25,0.5)[data[,3]+1]
						points(data[,1:2], col=cols, pch=pchs, cex=cexs, lwd=0.7)
					}
				n.trees = 100
				brt_model_scv1[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
					var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
					verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
				dev.copy2pdf(file=paste0("BRT_prediction_files/BRT_models/",analyses[i],"_SCV1_replicate_",j,".pdf")); dev.off()
				# AUCs[j,"Full_AUC"] = brt_model_scv1[[j]]$self.statistics$discrimination # AUC on the complete data set
				AUCs[j,"SCV1_AUC"] = brt_model_scv1[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)
				object = brt_model_scv1[[j]]; df = as.data.frame(rasters_stack)
				not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
				n.trees = brt_model_scv1[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
				prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
				rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
				
				# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
				spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,3:dim(data)[2]], proj4string=crs(nullRaster))
				myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRaster, k=n.folds, theRange=theRanges[1], selection="random")
				fold.vector = myblocks$foldID; n.trees = 100
				brt_model_scv2[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
					var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
					verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
				dev.copy2pdf(file=paste0("BRT_prediction_files/BRT_models/",analyses[i],"_SCV2_replicate_",j,".pdf")); dev.off()
				AUCs[j,"SCV2_AUC"] = brt_model_scv2[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)
				object = brt_model_scv2[[j]]; df = as.data.frame(rasters_stack)
				not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
				n.trees = brt_model_scv2[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
				prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
				rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
			}
		saveRDS(brt_model_ccvs, paste0("BRT_prediction_files/BRT_models/",analyses[i],"_models_CCV.rds"))
		saveRDS(brt_model_scv1, paste0("BRT_prediction_files/BRT_models/",analyses[i],"_models_SCV1.rds"))
		saveRDS(brt_model_scv2, paste0("BRT_prediction_files/BRT_models/",analyses[i],"_models_SCV2.rds"))
		write.csv(AUCs, paste0("BRT_prediction_files/BRT_models/",analyses[i],"_CCV_SCV_AUCs.csv"), row.names=F, quote=F)
	}

# A.4. Estimation of spatial sorting bias (SSB)

	# SSB = Dp/Da (Hijsmans 2012, Ecology), where:
		# Dp = mean distance between testing presence sites and nearest training-presence sites
		# Da = mean distance between testing absence sites and nearest training-presence sites
		# --> SSB = 1 suggests there is no spatial sorting bias
		# --> SSB = 0 suggests extreme spatial sorting bias

SSB_list = list()
for (i in 1:length(datasets))
	{
		data = datas[[i]]; n.folds = 5
		SSBs = matrix(nrow=n.folds, ncol=3)
		colnames(SSBs) = c("SSB_CCV","SSB_SCV1","SSB_SCV2")
		fold.vector = rep(NA, dim(data)[1])
		for (j in 1:dim(data)[1]) fold.vector[j] = sample(1:n.folds,1)
		for (j in 1:n.folds)
			{
				p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
				a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
				reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
				SSB = ssb(p, a, reference); SSBs[j,1] = SSB[1,"p"]/SSB[1,"a"]
			}
		folds_with_similar_sizes = FALSE; c = 0
		while (folds_with_similar_sizes == FALSE)
			{
				data_presence = data[which(data[,3]==1),]; c = c+1; # print(c)
				fivePoints = samplingPtsMinDist(data_presence[,1:2], minDist=200, nberOfPoints=n.folds)
				fold.vector = foldSelection(data[,1:2], selectedPoints=data_presence[fivePoints,1:2])
				fold.vector_presences = fold.vector[which(data[,3]==1)]
				counts = hist(fold.vector_presences, plot=F)$counts
				props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
				if (min(props) > (1/(n.folds+1))) folds_with_similar_sizes = TRUE
			}
		for (j in 1:n.folds)
			{
				p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
				a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
				reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
				SSB = ssb(p, a, reference); SSBs[j,2] = SSB[1,"p"]/SSB[1,"a"]
			}
		spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,3:dim(data)[2]], proj4string=crs(nullRaster))
		myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRaster, k=n.folds, theRange=theRanges[i], selection="random")
		for (j in 1:n.folds)
			{
				fold.vector = myblocks$foldID
				p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
				a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
				reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
				if (dim(reference)[1]>0)
					{
						SSB = ssb(p, a, reference); SSBs[j,3] = SSB[1,"p"]/SSB[1,"a"]
					}
			}
		write.table(round(SSBs,5), paste0("BRT_prediction_files/BRT_models/",analyses[i],"_CCV_SCV_SSBs.csv"), sep=",", row.names=F, quote=F)
		SSB_list[[i]] = SSBs
	}

# A.5. Comparison of the relative influence of each environmental factor

relativeInfluences = matrix(0, nrow=length(envVariables), ncol=length(analyses))
colnames(relativeInfluences) = analyses
for (i in 1:length(analyses))
	{
		brt_model_scv = readRDS(paste0("BRT_prediction_files/BRT_models/",analyses[i],"_models_SCV2.rds"))
		for (j in 1:length(brt_model_scv))
			{
				if ((i == 1)&(j == 1)) envVariableNames = rep(NA, length(envVariables))
				for (k in 1:length(envVariables))
					{
						if ((i == 1)&(j == 1)) envVariableNames[k] = names(envVariables[[k]])
						relativeInfluences[k,i] = relativeInfluences[k,i] + summary(brt_model_scv[[j]])[names(envVariables[[k]]),"rel.inf"]
					}
			}
		if (i == 1) row.names(relativeInfluences) = envVariableNames
		relativeInfluences[,i] = relativeInfluences[,i]/length(brt_model_scv)
	}
write.table(round(relativeInfluences,1), "Relative_influences.csv", quote=F, sep=",")

# A.6. Comparison of the response curves for each environmental factor

envVariableValues_list = list()
for (i in 1:length(analyses))
	{
		sp = SpatialPoints(datasets[[i]]); envVariableValues = matrix(nrow=3, ncol=length(envVariables))
		colnames(envVariableValues) = envVariableNames; row.names(envVariableValues) = c("median","minV","maxV")
		for (j in 1:length(envVariables))
			{
				points = rasterize(sp, envVariables[[j]])
				rast = envVariables[[j]]; rast[is.na(points)] = NA
				minV = min(rast[], na.rm=T); maxV = max(rast[], na.rm=T)
				envVariableValues[,j] = cbind(median(rast[], na.rm=T), minV, maxV)
			}
		envVariableValues_list[[i]] = envVariableValues
	}
	
envVariableValues_list = list()
for (i in 1:length(analyses))
	{
		data = datas[[i]]; data = data[which(data[,"response"]==1),]
		envVariableValues = matrix(nrow=3, ncol=length(envVariables))
		row.names(envVariableValues) = c("median","minV","maxV")
		colnames(envVariableValues) = envVariableNames
		for (j in 1:length(envVariables))
			{
				minV = min(data[,envVariableNames[j]], na.rm=T)
				maxV = max(data[,envVariableNames[j]], na.rm=T)
				medianV = median(data[,envVariableNames[j]], na.rm=T)
				envVariableValues[,j] = cbind(medianV, minV, maxV)
			}
		envVariableValues_list[[i]] = envVariableValues
	}

dev.new(width=8, height=1.5); par(mfrow=c(1,7), oma=c(1.5,1.5,1,1), mar=c(2,1,0.5,1), lwd=0.2, col="gray30")
for (i in c(1,2,3,5,7,8,9))
	{
		predictions_list = list(); dfs = list()
		for (j in 1:length(analyses))
			{
				valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[j]]["maxV",i]-envVariableValues_list[[j]]["minV",i])/100
				df = data.frame(matrix(nrow=length(seq(envVariableValues_list[[j]]["minV",i],envVariableValues_list[[j]]["maxV",i],valuesInterval)),ncol=length(envVariables)))
				colnames(df) = envVariableNames
				for (k in 1:length(envVariables))
					{
						valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[j]]["maxV",k]-envVariableValues_list[[j]]["minV",k])/100
						if (i == k) df[,envVariableNames[k]] = seq(envVariableValues_list[[j]]["minV",k],envVariableValues_list[[j]]["maxV",k],valuesInterval)
						if (i != k) df[,envVariableNames[k]] = rep(envVariableValues_list[[j]]["median",k],dim(df)[1])
					}
				dfs[[j]] = df; predictions = list()
				brt_model_scv = readRDS(paste0("BRT_prediction_files/BRT_models/",analyses[j],"_models_SCV2.rds"))
				for (k in 1:length(brt_model_scv))
					{
						n.trees = brt_model_scv[[k]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model_scv[[k]], newdata=df, n.trees, type, single.tree)
						if ((j == 1)&(k == 1))
							{
								minX = min(df[,envVariableNames[i]]); maxX = max(df[,envVariableNames[i]])
								minY = min(prediction); maxY = max(prediction)
							}	else	{
								if (minX > min(df[,envVariableNames[i]])) minX = min(df[,envVariableNames[i]])
								if (maxX < max(df[,envVariableNames[i]])) maxX = max(df[,envVariableNames[i]])
								if (minY > min(prediction)) minY = min(prediction)
								if (maxY < max(prediction)) maxY = max(prediction)
							}
						predictions[[k]] = prediction
					}
				predictions_list[[j]] = predictions
			}
		cols = c("red","chartreuse3")
		for (k in 1:length(brt_model_scv))
			{
				for (j in 1:length(analyses))
					{
						if ((j == 1)&(k == 1))
							{
								plot(dfs[[j]][,envVariableNames[i]],predictions_list[[j]][[k]],col=cols[j],ann=F,axes=F,lwd=0.2,type="l",xlim=c(minX,maxX),ylim=c(minY,maxY))
							}	else	{
								lines(dfs[[j]][,envVariableNames[i]],predictions_list[[j]][[k]],col=cols[j],lwd=0.2)
							}
					}
			}
		axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.07,0))
		axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.2,0))
		title(ylab="predicted values", cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
		title(xlab=envVariableNames[i], cex.lab=0.9, mgp=c(0.9,0,0), col.lab="gray30")
		text = paste0("(",round(relativeInfluences[i,1],1),"%, ",round(relativeInfluences[i,2],1),"%)")
		mtext(text, side=1, line=1.9, col="gray30", cex=0.6)
		box(lwd=0.2, col="gray30")
	}

# A.7. BRT predictions on historical and future scenarios

years = c(2030, 2050, 2070)
year_intervals = c("2021_2040","2041_2060","2061_2080") 
scenarios = c("Historical","RCP_26","RCP_60","RCP_85")
models = c("GFDL-ESM2M","HadGEM2-ES","IPSL-CM5A-LR","MIROC5")
variables = c("population","pr","tas","landcover")
rasters_stacks = list(); envVariables_list = list(); c = 0
for (s in 1:length(scenarios))
	{
		for (y in 1:length(years))
			{
				for (m in 1:length(models))
					{
						files = list.files(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/"))
						if (s != 1) files = files[grepl(year_intervals[y],files)]
						files = files[!grepl("5min",files)]
						index_population = which(grepl("population",files)); index_temperature = which(grepl("tas_day",files))
						index_precipitation = which(grepl("pr_day",files)); index_land_cover = which(grepl("landcover",files))
						population = raster(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_population]))
						temperature = raster(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_temperature]))
						precipitation = raster(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_precipitation]))
						land_cover = nc_open(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_land_cover]))
						names(population) = "population"; names(temperature) = "temperature"; names(precipitation) = "precipitation"
						landCoverVariableIDs = names(land_cover$var); cols = list()
						landCoverVariableNames = as.character(read.csv("Environmental_rasters/LC_vars.csv")[1:12,2])
						land_covers1 = list(); land_covers2 = list(); land_covers3 = list()
						for (i in 2:13)
							{
								land_covers1[[i-1]] = brick(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",
															files[index_land_cover]), varname=landCoverVariableIDs[i])
							}
						for (i in 1:length(land_covers1)) names(land_covers1[[i]]) = landCoverVariableNames[i]
						variable_codes = c("croplands","pastures","urbanAreas","primaryForest","primaryNonF","secondaryForest","secondaryNonF")
						variable_names = c("crops","pasture","urban land","forested primary land","non-forested primary land",
										   "potentially forested secondary land","potentially non-forested secondary land")
						for (i in 1:length(variable_names))
							{
								names = gsub("\\."," ",landCoverVariableNames); indices = which(landCoverVariableNames==variable_names[i])
								if (length(indices) == 0) indices = which(grepl(variable_names[i],names))
								if (variable_names[i] == "pasture") indices = c(indices, which(grepl("rangeland",names)))
								land_cover = land_covers1[[indices[1]]]; names(land_cover) = variable_codes[i]
								if (length(indices) > 1)
									{
										for (j in 2:length(indices)) land_cover[] = land_cover[]+land_covers1[[indices[j]]][]
									}
								land_covers2[[i]] = land_cover[[1]]; land_covers3[[i]] = raster::aggregate(land_cover[[1]],2)
							}
						envVariables = list()
						envVariables[[1]] = temperature; envVariables[[2]] = precipitation
						envVariables[[3]] = land_covers3[[4]] # primary forest areas
						envVariables[[4]] = land_covers3[[5]] # primary non-forest areas
						envVariables[[5]] = land_covers3[[6]] # secondary forest areas
						envVariables[[6]] = land_covers3[[7]] # secondary non-forest areas
						envVariables[[7]] = land_covers3[[1]] # croplands (all catergories)
						envVariables[[8]] = land_covers3[[2]] # managed pasture + rangeland
						pLog = population; pLog[] = log10(pLog[]+1); envVariables[[9]] = pLog
						for (i in 1:length(envVariables)) envVariables[[i]][is.na(mask[])] = NA
						for (i in 1:length(envVariables)) envVariables[[i]] = crop(envVariables[[i]], africa2, snap="out")
						for (i in 1:length(envVariables)) envVariables[[i]] = mask(envVariables[[i]], africa2)
						envVariables[[1]][] = envVariables[[1]][]-273.15 # conversion to Celcius degrees
						envVariables[[2]][] = envVariables[[2]][]*60*60*24 # conversion to kg/m2/day
						for (i in c(1,2,9)) envVariables[[i]][is.na(envVariables[[3]][])] = NA
						c = c+1; rasters_stacks[[c]] = stack(envVariables); envVariables_list[[c]] = envVariables
					}
			}
	}

predictions_list1 = list(); std_deviations_list1 = list(); human_population_list1 = list()
for (i in 1:length(analyses))
	{
		predictions_list2 = list(); std_deviations_list2 = list(); human_population_list2 = list(); c = 0
		for (s in 1:length(scenarios))
			{
				predictions1 = list(); std_deviations = list(); human_populations = list()
				brt_model_scv = readRDS(paste0("BRT_prediction_files/BRT_models/",analyses[i],"_models_SCV2.rds"))
				for (y in 1:length(years))
					{
						predictions2 = list(); human_population = list()
						for (m in 1:length(models))
							{
								c = c+1; rasters_stack = rasters_stacks[[c]]; replicates = list()
								for (j in 1:length(brt_model_scv))
									{
										df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
										n.trees = brt_model_scv[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
										prediction = predict.gbm(brt_model_scv[[j]], newdata, n.trees, type, single.tree)
										rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction; replicates[[j]] = rast
									}
								rasts = stack(replicates); prediction = mean(rasts); predictions2[[m]] = prediction
							}
						prediction = mean(stack(predictions2)); std_deviation = calc(stack(predictions2), fun=sd)
						names(prediction) = paste0(analyses[i],"_",years[y]); predictions1[[y]] = prediction
						names(std_deviation) = paste0(analyses[i],"_",years[y]); std_deviations[[y]] = std_deviation
						if (scenarios[s] == "Historical")
							{
								writeRaster(prediction, paste0("BRT_prediction_files/BRT_predictions/",analyses[i],"_historical.asc"), overwrite=T)
							}	else	{
								writeRaster(prediction, paste0("BRT_prediction_files/BRT_predictions/",analyses[i],"_",scenarios[s],"_",years[y],".asc"), overwrite=T)
							}
						for (m in 1:length(models))
							{
								files = list.files(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/"))
								if (s != 1) files = files[grepl(year_intervals[y],files)]
								files = files[!grepl("5min",files)]; index_population = which(grepl("population",files))
								human_pop = raster(paste0("Environmental_rasters/",scenarios[s],"/",models[m],"/",files[index_population]))
								human_pop[is.na(mask[])] = NA; human_population[[m]] = mask(crop(human_pop,africa2,snap="out"),africa2)
							}
						human_populations[[y]] = mean(stack(human_population))
					}
				predictions_list2[[s]] = predictions1; std_deviations_list2[[s]] = std_deviations; human_population_list2[[s]] = human_populations
			}
		predictions_list1[[i]] = predictions_list2; std_deviations_list1[[i]] = std_deviations_list2; human_population_list1[[i]] = human_population_list2
	}

suffixes = c("standard_deviation","ENM_difference","index_of_human_exposure"); allIHEs = FALSE
for (h in 1:length(suffixes)) {
plotStdDeviations = FALSE; plotENMdifferences = FALSE; plotHumanExposures = FALSE
if (suffixes[h] == "standard_deviation") plotStdDeviations = TRUE
if (suffixes[h] == "ENM_difference") plotENMdifferences = TRUE
if (suffixes[h] == "index_of_human_exposure") plotHumanExposures = TRUE
for (i in 1:length(analyses))
	{
		for (s in 2:length(scenarios))
			{
				rasters = list()
				rasters[[1]] = predictions_list1[[i]][[1]][[1]] # "historical" prediction
				for (j in 1:length(predictions_list1[[i]][[s]])) rasters[[j+1]] = predictions_list1[[i]][[s]][[j]]
				c = length(rasters)
				if (plotStdDeviations == TRUE)
					{
						for (j in 1:length(std_deviations_list1[[i]][[s]]))
							{
								c = c+1; rasters[[c]] = std_deviations_list1[[i]][[s]][[j]]
							}
					}
				if (plotENMdifferences == TRUE)
					{
						for (j in 1:length(predictions_list1[[i]][[s]]))
							{
								rast1 = predictions_list1[[i]][[s]][[1]]; rast2 = predictions_list1[[i]][[s]][[j]]
								rast3 = rast2; rast3[] = rast3[]-rast1[]; c = c+1; rasters[[c]] = rast3
								# cat(c(min(rasters[[c]][],na.rm=T),max(rasters[[c]][],na.rm=T),"\n"))
							}
					}
				if (plotHumanExposures == TRUE)
					{
						for (j in 1:length(predictions_list1[[i]][[s]]))
							{
								rast1 = predictions_list1[[i]][[s]][[j]]; rast1[rast1[]>=0.5] = 1; rast1[rast1[]<0.5] = 0; rast2 = human_population_list1[[i]][[s]][[j]]
								rast3 = rast1; rast3[] = rast3[]*rast2[]; rast3[] = log10(rast3[]+1); c = c+1; rasters[[c]] = rast3
							}
						if (allIHEs == TRUE)
							{
								rast1 = predictions_list1[[i]][[1]][[1]]; rast1[rast1[]>=0.5] = 1; rast1[rast1[]<0.5] = 0; rast2 = human_population_list1[[i]][[1]][[1]]
								rast3 = rast1; rast3[] = rast3[]*rast2[]; rast3[] = log10(rast3[]+1); c = c+1; rasters[[c]] = rast3 # for the current period
							}
					}
				if (scenarios[[s]] == "RCP_26") scenario = "(RCP 2.6)"
				if (scenarios[[s]] == "RCP_60") scenario = "(RCP 6.0)"
				if (scenarios[[s]] == "RCP_85") scenario = "(RCP 8.5)"
				legend1 = raster(as.matrix(c(0,1)))
				showingPlots = TRUE; if (showingPlots == TRUE) {
				pdf(paste0(analyses[i],"_",scenarios[s],"_",suffixes[h],".pdf"), width=8, height=3) # dev.new(width=8, height=3)
				par(mfrow=c(2,7), oma=c(0,0,1.5,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
				if (allIHEs == FALSE)
					{
						plot(rasters[[1]], col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(rasters[[1]][],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext("Current period (t0)", side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
					 		smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", 
					 		line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3, horizontal=T)
					}
				plot(rasters[[2]], col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(rasters[[2]][],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
				plot(africa3, add=T, border="gray50", lwd=0.5)
				mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[1]," ",scenario), side=3, line=-0.7, cex=0.65, col="gray30")
				plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
					 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", 
					 line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3, horizontal=T)
				plot(rasters[[3]], col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(rasters[[3]][],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
				plot(africa3, add=T, border="gray50", lwd=0.5)
				mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[2]," ",scenario), side=3, line=-0.7, cex=0.65, col="gray30")
				plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
					 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", 
					 line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3, horizontal=T)
				plot(rasters[[4]], col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(rasters[[4]][],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
				plot(africa3, add=T, border="gray50", lwd=0.5)
				mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[3]," ",scenario), side=3, line=-0.7, cex=0.65, col="gray30")
				plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
					 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", 
					 col.axis="gray30", line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3, horizontal=T)
				if (plotStdDeviations == TRUE)
					{
						plot(rasters[[5]], col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(rasters[[2]][],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[1]," (SD)"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
							 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30",
							 col.axis="gray30", line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3, horizontal=T)
						plot(rasters[[6]], col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(rasters[[3]][],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[2]," (SD)"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
							 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", 
							 col.axis="gray30", line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3, horizontal=T)
						plot(rasters[[7]], col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131][1:(max(rasters[[4]][],na.rm=T)*100)], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[3]," (SD)"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend1, col=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(131))[21:131], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
							 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", 
							 col.axis="gray30", line=0, mgp=c(0,0.0,0), at=seq(0,1,0.25), labels=c("0","0.25","0.5","0.75","1")), alpha=1, side=3, horizontal=T)
					}
				if (plotENMdifferences == TRUE)
					{
						legend2 = raster(as.matrix(c(-0.45,0.45)))
						if (i == 1)
							{
								index1 = 45+(min(rasters[[5]][],na.rm=T)*90); index2 = 45+(max(rasters[[5]][],na.rm=T)*90)
								ats = seq(-0.25,0.25,0.25); labels = c("-0.25","0","0.25")
							} 
						if (i == 2)
							{
								index1 = 45+(min(rasters[[5]][],na.rm=T)*45); index2 = 45+(max(rasters[[5]][],na.rm=T)*45)
								ats = seq(-0.25,0.25,0.25); labels = c("-0.50","0","0.50")
							} 
						cols = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(91)); cols[43:47] = "#F7F7F7"
						plot(rasters[[5]], col=cols[index1:index2], ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[1]," ",scenario," - t0"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend2, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3,
							 axis.args=list(cex.axis=0.6, lwd=0,lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0,
							 mgp=c(0,0.0,0), at=ats, labels=labels), alpha=1, side=3, horizontal=T)
						if (i == 1)
							{
								index1 = 45+(min(rasters[[6]][],na.rm=T)*90); index2 = 45+(max(rasters[[6]][],na.rm=T)*90)
								ats = seq(-0.25,0.25,0.25); labels = c("-0.25","0","0.25")
							}
						if (i == 2)
							{
								index1 = 45+(min(rasters[[6]][],na.rm=T)*45); index2 = 45+(max(rasters[[6]][],na.rm=T)*45)
								ats = seq(-0.25,0.25,0.25); labels = c("-0.50","0","0.50")
							}
						cols = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(91)); cols[43:47] = "#F7F7F7"
						plot(rasters[[6]], col=cols[index1:index2], ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[2]," ",scenario," - t0"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend2, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3,
							 axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0,
							 mgp=c(0,0.0,0), at=ats, labels=labels), alpha=1, side=3, horizontal=T)
						if (i == 1)
							{
								index1 = 45+(min(rasters[[7]][],na.rm=T)*90); index2 = 45+(max(rasters[[7]][],na.rm=T)*90)
								ats = seq(-0.25,0.25,0.25); labels = c("-0.25","0","0.25")
							}
						if (i == 2)
							{
								index1 = 45+(min(rasters[[7]][],na.rm=T)*45); index2 = 45+(max(rasters[[7]][],na.rm=T)*45)
								ats = seq(-0.25,0.25,0.25); labels = c("-0.50","0","0.50")
							}
						cols = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(91)); cols[43:47] = "#F7F7F7"
						plot(rasters[[7]], col=cols[index1:index2], ann=F, legend=F, axes=F, box=F); plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[3]," ",scenario," - t0"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend2, col=cols, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3,
							 axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", col.axis="gray30", line=0,
							 mgp=c(0,0.0,0), at=ats, labels=labels), alpha=1, side=3, horizontal=T) }
					}
				if (plotHumanExposures == TRUE)
					{
						legend3 = raster(as.matrix(c(0,max(rasters[[7]][],na.rm=T))))
						if (allIHEs == TRUE)
							{
								plot(rasters[[8]], col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81][1:(max(rasters[[8]][],na.rm=T)*10)], ann=F, legend=F, axes=F, box=F)
								plot(africa3, add=T, border="gray50", lwd=0.5)
								mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext("Current period (t0)", side=3, line=-0.7, cex=0.65, col="gray30")
								plot(legend3, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
							 		 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30",
							 		 col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T)
							}
						plot(rasters[[5]], col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81][1:(max(rasters[[7]][],na.rm=T)*10)], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[1]," (IHE)"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend3, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
							 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30",
							 col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T)
						plot(rasters[[6]], col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81][1:(max(rasters[[7]][],na.rm=T)*10)], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[2]," (IHE)"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend3, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
							 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", 
							 col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T)
						plot(rasters[[7]], col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81][1:(max(rasters[[7]][],na.rm=T)*10)], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext("", side=3, line=0.3, cex=0.65, col="gray30"); mtext(paste0(years[3]," (IHE)"), side=3, line=-0.7, cex=0.65, col="gray30")
						plot(legend3, col=colorRampPalette(brewer.pal(9,"YlOrBr"))(111)[1:81], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, 
							 smallplot=c(0.10,0.80,0.03,0.06), adj=3, axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.6, col="gray30", 
							 col.axis="gray30", line=0, mgp=c(0,0.0,0)), alpha=1, side=3, horizontal=T)
						writeRaster(rasters[[5]], paste0("BRT_prediction_files/BRT_predictions/LASV_IHE05_",scenarios[s],"_",years[1],".asc"), overwrite=T)
						writeRaster(rasters[[6]], paste0("BRT_prediction_files/BRT_predictions/LASV_IHE05_",scenarios[s],"_",years[2],".asc"), overwrite=T)
						writeRaster(rasters[[7]], paste0("BRT_prediction_files/BRT_predictions/LASV_IHE05_",scenarios[s],"_",years[3],".asc"), overwrite=T)
						if (allIHEs == TRUE) writeRaster(rasters[[8]], paste0("BRT_prediction_files/BRT_predictions/LASV_IHE05_",scenarios[s],"_current.asc"), overwrite=T)
					}
				dev.off() # dev.copy2pdf(file=paste0(analyses[i],"_",years[y],".pdf"))
			}
	}}

# A.8. Plotting specific environmental factor predictions

selectedVariables = c(1,2,8); cols = list()
variableNames = c("Temperature","Precipitation","Pastures")
cols[[1]] = colorRampPalette(brewer.pal(9,"YlOrRd"))(150)[1:100]
cols[[2]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
cols[[3]] = colorRampPalette(c("gray97","burlywood3"),bias=1)(100)
cols[[4]] = rev(colorRampPalette(brewer.pal(11,"RdYlGn"))(100))
scenarioNames = c("Historical","RCP 2.6","RCP 6.0","RCP 8.5")
for (i in 1:length(selectedVariables))
	{
		c = 0; rast_list1 = list(); rast_list2 = list()
		for (s in 1:length(scenarios))
			{
				rasts = list()
				for (y in 1:length(years))
					{
						for (m in 1:length(models))
							{
								c = c+1
								if (m == 1) r = envVariables_list[[c]][[selectedVariables[i]]]
								if (m >= 2) r[] = r[] + envVariables_list[[c]][[selectedVariables[i]]][]
							}
						r[!is.na(r[])] = r[!is.na(r[])]/length(models); rasts[[y]] = r
						if ((s == 1) & (y == 1))
							{
								minV1 = min(rasts[[y]][],na.rm=T); maxV1 = max(rasts[[y]][],na.rm=T)
							}	else	{
								if (minV1 > min(rasts[[y]][],na.rm=T)) minV1 = min(rasts[[y]][],na.rm=T)
								if (maxV1 < max(rasts[[y]][],na.rm=T)) maxV1 = max(rasts[[y]][],na.rm=T)
							}
					}
				rast_list1[[s]] = rasts
			}
		for (s in 2:length(scenarios))
			{
				rasts = list(); rasts[[1]] = rast_list1[[1]][[1]]
				for (y in 1:length(years)) rasts[[1+y]] = rast_list1[[s]][[y]]
				for (y in 1:length(years))
					{
						r = rast_list1[[s]][[y]]; r[] = r[]-rast_list1[[1]][[1]][]; rasts[[4+y]] = r
						if ((s == 2) & (y == 1))
							{
								minV2 = min(rasts[[4+y]][],na.rm=T); maxV2 = max(rasts[[4+y]][],na.rm=T)
							}	else		{
								if (minV2 > min(rasts[[4+y]][],na.rm=T)) minV2 = min(rasts[[4+y]][],na.rm=T)
								if (maxV2 < max(rasts[[4+y]][],na.rm=T)) maxV2 = max(rasts[[4+y]][],na.rm=T)
							}
					}
				rast_list2[[s]] = rasts
			}
		for (s in 2:length(scenarios))
			{
				pdf(paste0(variableNames[i],"_",scenarios[s],".pdf"), width=8, height=3) # dev.new(width=8, height=3)
				par(mfrow=c(2,7), oma=c(0,0,1.5,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
				for (j in 1:4)
					{
						index1 = round(((min(rast_list2[[s]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
						index2 = round(((max(rast_list2[[s]][[j]][],na.rm=T)-minV1)/(maxV1-minV1))*100)
						plot(rast_list2[[s]][[j]], col=cols[[i]][index1:index2], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						if (j == 1) mtext("Current period (t0)", side=3, line=-0.7, cex=0.65, col="gray30")
						if (j >= 2) mtext(paste0(years[j-1]," (",scenarioNames[s],")"), side=3, line=-0.7, cex=0.65, col="gray30")
						rastLegend = raster(t(as.matrix(c(minV1,maxV1))))
						plot(rastLegend, col=cols[[i]], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3,
							 axis.args=list(cex.axis=0.6,lwd=0,lwd.tick=0.2,col.tick="gray30",tck=-0.6,col="gray30",col.axis="gray30",line=0,mgp=c(0,0.0,0)),
							 alpha=1,side=3,horizontal=T)
					}
				for (j in 5:7)
					{
						index1 = round(((min(rast_list2[[s]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
						index2 = round(((max(rast_list2[[s]][[j]][],na.rm=T)-minV2)/(maxV2-minV2))*100)
						plot(rast_list2[[s]][[j]], col=cols[[4]][index1:index2], ann=F, legend=F, axes=F, box=F)
						plot(africa3, add=T, border="gray50", lwd=0.5)
						mtext(paste0(years[j-4]," (",scenarioNames[s],") - t0"), side=3, line=-0.7, cex=0.65, col="gray30")
						rastLegend = raster(t(as.matrix(c(minV2,maxV2))))
						plot(rastLegend, col=cols[[4]], legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.80,0.03,0.06), adj=3,
							  axis.args=list(cex.axis=0.6,lwd=0,lwd.tick=0.2,col.tick="gray30",tck=-0.6,col="gray30",col.axis="gray30",line=0,mgp=c(0,0.0,0)),
							  alpha=1,side=3,horizontal=T)
					}
				dev.off()
			}
	}

# B. PHYLOGEOGRAPHIC AND PHYLODYNAMIC ANALYSES

# B1. Preparation of LASV sequences

setwd(paste(wd,wdb1,sep="/"))

	# B1.1. Renaming sequence IDs to match with metadata files

metadata = read.csv("Original_sequence_files/LASV_all_the_metadata.csv", head=T, colClasses="character")
for (i in 1:length(segments))
	{
		for (j in 1:dim(metadata)[1])
			{
				if ((metadata[j,paste0("accession_",segments[i])]!="")
					&(grepl("\\|",metadata[j,paste0("segment_",segments[i])]))
					&(!grepl("LASV|NGA|2019",metadata[j,paste0("segment_",segments[i])])))
					{
						accession = metadata[j,paste0("segment_",segments[i])]
						accession = unlist(strsplit(accession,"\\|"))[1]
						if (accession != metadata[j,paste0("accession_",segments[i])])
							{
								cat(segments[i],j,metadata[j,paste0("accession_",segments[i])],metadata[j,paste0("segment_",segments[i])],"\n",sep=" ")
							}
					}
			}
	}
for (i in 1:length(segments))
	{
		fasta = scan(paste0("Original_sequence_files/LASV_",segments[i],"_alignment_SD.fasta"), what="", sep="\n", quiet=T)
		colNames = c(paste0("segment_",segments[i]),paste0("accession_",segments[i]),"country","admin1","location","latitude","longitude",paste0("kept_in_",segments[i],"_aln"))
		tab1 = metadata[,colNames]; tab1 = tab1[which(tab1[,paste0("segment_",segments[i])]!=""),]
		sequenceIDs = gsub(">","",fasta[which(grepl(">",fasta))])
		missingSequences = which((!tab1[,paste0("accession_",segments[i])]%in%sequenceIDs)&(tab1[,paste0("kept_in_",segments[i],"_aln")]!="0"))
		if (length(tab1[missingSequences,paste0("accession_",segments[i])]) > 0)
			{
				cat(tab1[missingSequences,paste0("accession_",segments[i])],sep="\n")
			}
		if (length(sequenceIDs[which(!sequenceIDs%in%tab1[,paste0("accession_",segments[i])])]) > 0)
			{
				cat(sequenceIDs[which(!sequenceIDs%in%tab1[,paste0("accession_",segments[i])])],sep="\n")
			}
		tab1 = tab1[which(tab1[,paste0("accession_",segments[i])]%in%sequenceIDs),1:7]
		row.names(tab1) = tab1[,paste0("accession_",segments[i])]; tab1 = tab1[sequenceIDs,]
		write.csv(tab1, paste0("LASV_",segments[i],"_alignment_1.csv"), quote=F, row.names=F)
		for (j in 1:length(fasta))
			{
				if (grepl(">",fasta[j]))
					{
						seqID = gsub(">","",fasta[j])
						index = which(tab1[,paste0("accession_",segments[i])]==seqID)
						if (length(index) != 1) print(c(i,j)) 
						fasta[j] = paste0(">",tab1[index,paste0("segment_",segments[i])])
					}
			}
		write(fasta, paste0("LASV_",segments[i],"_alignment_1.fasta"))
	}

	# B1.2. Identifying and discarding suspiciouly duplicated sequences

mismatchesMatrix = function(sequences)
	{
		mmd = matrix(nrow=length(sequences), ncol=length(sequences))
		for (i in 1:length(sequences))
			{
				for (j in 1:i)
					{
						if (i == j)
							{
								mmd[j,i] = 0
							}	else	{
								seq1 = unlist(strsplit(sequences[j],""))
								seq2 = unlist(strsplit(sequences[i],""))
								seq2 = seq2[which((seq1=="A")|(seq1=="C")|(seq1=="G")|(seq1=="T")|(seq1=="a")|(seq1=="c")|(seq1=="g")|(seq1=="t"))]
								seq1 = seq1[which((seq1=="A")|(seq1=="C")|(seq1=="G")|(seq1=="T")|(seq1=="a")|(seq1=="c")|(seq1=="g")|(seq1=="t"))]
								seq1 = seq1[which((seq2=="A")|(seq2=="C")|(seq2=="G")|(seq2=="T")|(seq2=="a")|(seq2=="c")|(seq2=="g")|(seq2=="t"))]
								seq2 = seq2[which((seq2=="A")|(seq2=="C")|(seq2=="G")|(seq2=="T")|(seq2=="a")|(seq2=="c")|(seq2=="g")|(seq2=="t"))]
								if (length(seq1) > 0)
									{
										mm = sum(seq1!=seq2)
										mmd[j,i] = mm; mmd[i,j] = mm
									}	else	{
										mmd[j,i] = NA; mmd[i,j] = NA
									}
							}
					}
			}
		return(mmd)
	}
for (i in 1:length(segments))
	{
		fasta = scan(paste0("LASV_",segments[i],"_alignment_1.fasta"), what="", sep="\n", quiet=T)
		seqIDs = fasta[which(grepl(">",fasta))]; sequences = fasta[which(!grepl(">",fasta))]
		mismatches = mismatchesMatrix(sequences); row.names(mismatches) = gsub(">","",seqIDs)
		write.csv(mismatches, paste0("LASV_",segments[i],"_mismatches.csv"), quote=F)
	}
mismatches_L = read.csv("LASV_L_mismatches.csv", header=T)
mismatches_S = read.csv("LASV_S_mismatches.csv", header=T)
pdf("Mismatch_distributions_NEW.pdf", width=8, height=4); datasets = c(); # dev.new(width=8, height=4)
par(mfrow=c(1,2), mgp=c(1,0.35,0), oma=c(1,0.5,1,2), mar=c(2.0,3,0,0))
hist(mismatches_L[lower.tri(mismatches_L)], breaks=1000, col="gray50", border=NA, axes=F, ann=F)
axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.0,0), lwd=0.2, tck=-0.013, col.axis="gray30")
axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.014, col.axis="gray30")
title(xlab="mismatches", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
title(ylab="frequency", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
title(main="Alignment for L segment", cex.main=0.7, col.main="gray30", line=-1.7)
hist(mismatches_S[lower.tri(mismatches_S)], breaks=1000, col="gray50", border=NA, axes=F, ann=F)
axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.0,0), lwd=0.2, tck=-0.013, col.axis="gray30")
axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.014, col.axis="gray30")
title(xlab="mismatches", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
title(ylab="frequency", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
title(main="Alignment for S segment", cex.main=0.7, col.main="gray30", line=-1.7)
dev.off()

for (s in 1:length(segments))
	{
		fasta = scan(paste0("LASV_",segments[s],"_alignment_1.fasta"), what="", sep="\n", quiet=T)
		tab = read.csv(paste0("LASV_",segments[s],"_alignment_1.csv"), header=T)
		remainingSuspiciousPairs = TRUE; N = 0
		while (remainingSuspiciousPairs == TRUE)
			{
				seqIDs = fasta[which(grepl(">",fasta))]; sequences = fasta[which(!grepl(">",fasta))]
				mismatches = mismatchesMatrix(sequences); row.names(mismatches) = gsub(">","",seqIDs)
				suspicous_pairs_list = list(); suspicous_pairs = c(); nS = list(); N = N+1
				for (m in 0:5)
					{
						n = 0
						if (m == 0) suspicous_pairs = c(suspicous_pairs, paste0(m," NUCLEOTIDE MISMATCH ##############################"), " ")
						if (m == 1) suspicous_pairs = c(suspicous_pairs, " ", paste0(m," NUCLEOTIDE MISMATCHE #############################"), " ")
						if (m >= 2) suspicous_pairs = c(suspicous_pairs, " ", paste0(m," NUCLEOTIDE MISMATCHES ############################"), " ")
						for (i in 2:dim(mismatches)[1])
							{
								for (j in 1:(i-1))
									{
										if ((!is.na(mismatches[i,j]))&&(mismatches[i,j] == m))
											{
												if (m <= 5) suspicous_pairs_list[[length(suspicous_pairs_list)+1]] = c(row.names(mismatches)[i],row.names(mismatches)[j])
						suspicous_pairs = c(suspicous_pairs, paste(row.names(mismatches)[i],as.character(tab[which(tab[,1]==row.names(mismatches)[i]),"location"]),sep=" - "))
						suspicous_pairs = c(suspicous_pairs, paste(row.names(mismatches)[j],as.character(tab[which(tab[,1]==row.names(mismatches)[j]),"location"]),sep=" - "))
												suspicous_pairs = c(suspicous_pairs, " "); n = n + 1; # cat("\n")
											}
									}
							}
						nS[[m+1]] = n
					}
				if (N == 1) write(suspicous_pairs, paste0("LASV_",segments[s],"_alignment_1.txt"))
				if (length(suspicous_pairs_list) == 0)
					{
						remainingSuspiciousPairs = FALSE
					}	else	{
						sequences_to_remove = c(); sequences_to_keep = c()
						for (i in 1:length(suspicous_pairs_list))
							{
								seqIDs = c(); index = NA
								if (sum(grepl("_NGA-",suspicous_pairs_list[[i]])) == 1) # 1 sequence from the "Science" data set
									{
										index = which(grepl("_NGA-",suspicous_pairs_list[[i]]))
									}
								if (sum(grepl("_NGA-",suspicous_pairs_list[[i]])) == 2) # 2 sequences from the "Science" data set
									{
										date1 = unlist(strsplit(suspicous_pairs_list[[i]][1],"_")); date1 = decimal_date(ymd(date1[length(date1)]))
										date2 = unlist(strsplit(suspicous_pairs_list[[i]][2],"_")); date2 = decimal_date(ymd(date2[length(date2)]))
										if (date1 == date2)
											{
												accession1 = unlist(strsplit(suspicous_pairs_list[[i]][1],"\\|"))[1]
												accession2 = unlist(strsplit(suspicous_pairs_list[[i]][2],"\\|"))[1]
												accessions = c(accession1, accession2); accessions = accessions[order(accessions)]
												if (accession1 == accessions[1]) index = 1
												if (accession2 == accessions[1]) index = 2
											}	
										if (date1 < date2) index = 1
										if (date2 < date1) index = 2
									}
								if (sum(grepl("_NGA-",suspicous_pairs_list[[i]])) == 0) # no sequence from the "Science" data set
									{
										accession1 = unlist(strsplit(suspicous_pairs_list[[i]][1],"\\|"))[1]
										accession2 = unlist(strsplit(suspicous_pairs_list[[i]][2],"\\|"))[1]
										accessions = c(accession1, accession2); accessions = accessions[order(accessions)]
										if (accession1 == accessions[1]) index = 1
										if (accession2 == accessions[1]) index = 2
									}
								if (is.na(index)) print(i)
								sequences_to_remove = c(sequences_to_remove, suspicous_pairs_list[[i]][-index])
								sequences_to_keep = c(sequences_to_keep, suspicous_pairs_list[[i]][index])
							}
						buffer = c()
						for (i in 1:length(fasta))
							{
								if (grepl(">",fasta[i]))
									{
										seqID = gsub(">","",fasta[i])
										if (!seqID%in%sequences_to_remove)
											{
												buffer = c(buffer,fasta[i],fasta[i+1])
											}
									}	
							}
						fasta = buffer
						tab = tab[which(!tab[,1]%in%sequences_to_remove),]
					}
			}
		write(fasta, paste0("LASV_",segments[s],"_alignment_2.fasta"))
		write.csv(tab, paste0("LASV_",segments[s],"_alignment_2.csv"), quote=F, row.names=F)
	}
for (i in 1:length(segments)) # to remove CIV and TGO sequences (third, non-analysed, study area)
	{
		fasta2 = scan(paste0("LASV_",segments[i],"_alignment_2.fasta"), what="", sep="\n", quiet=T)
		tab2 = read.csv(paste0("LASV_",segments[i],"_alignment_2.csv"), header=T); buffer = c()
		for (j in 1:length(fasta2))
			{
				if (grepl(">",fasta2[j]))
					{
						seqID = gsub(">","",fasta2[j])
						if (seqID != "KT833204|Sil03|GIN|Guinea|2013-12-6") # wrong coordinates for Guinea
							{
								if ((!grepl("\\|CIV\\|",seqID))&(!grepl("\\|GHA\\|",seqID))
									&(!grepl("\\|TGO\\|",seqID))&(!grepl("Togo",seqID)))
									{
										buffer = c(buffer,fasta2[j],fasta2[j+1])
									}
							}
					}	
			}
		tab2 = tab2[which(tab2[,1]!="KT833204|Sil03|GIN|Guinea|2013-12-6"),]
		tab2 = tab2[which((!grepl("\\|CIV\\|",tab2[,1]))&(!grepl("\\|TGO\\|",tab2[,1]))&(!grepl("Togo",tab2[,1]))),]
		write(buffer, paste0("LASV_",segments[i],"_alignment_2.fasta"))
		write.csv(tab2, paste0("LASV_",segments[i],"_alignment_2.csv"), quote=F, row.names=F)
	}

	# B1.3. Preparing alignments for the continuous phylogeographic analyses

for (i in 1:length(segments))
	{
		fasta2 = scan(paste0("LASV_",segments[i],"_alignment_2.fasta"), what="", sep="\n", quiet=T)
		tab2 = read.csv(paste0("LASV_",segments[i],"_alignment_2.csv"), header=T)
		tab3 = tab2; tab3 = tab3[which(tab3[,"longitude"]!=""),]; fasta3 = c()
		for (j in 1:length(fasta2))
			{
				if (grepl(">",fasta2[j]))
					{
						seqID = gsub(">","",fasta2[j])
						if (seqID%in%tab3[,1])
							{
								fasta3 = c(fasta3,fasta2[j],fasta2[j+1])
							}
					}	
			}		
		write(fasta3, paste0("LASV_",segments[i],"_alignment_3.fasta"))
		write.csv(tab3, paste0("LASV_",segments[i],"_alignment_3.csv"), quote=F, row.names=F)
		tab4 = tab3[,c(paste0("segment_",segments[i]),"latitude","longitude")]
		colnames(tab4) = c("trait","latitude","longitude")
		write.table(tab4, paste0("LASV_",segments[i],"_alignment_3.txt"), sep="	", row.names=F, quote=F)
	}

	# B1.4. Preparing alignments that only contain Nigerian sequences

for (i in 1:length(segments))
	{
		fasta1 = scan(paste0("LASV_",segments[i],"_alignment_1.fasta"), what="", sep="\n", quiet=T)
		tab1 = read.csv(paste0("LASV_",segments[i],"_alignment_1.csv"), header=T)
		lasv1 = read.csv(paste0("LASV_",segments[i],"_seqs_LASV1.csv"), header=F)
		tab4a = tab1[which(grepl("\\|2019\\|",as.character(tab1[,1]))),]
		tab4b = tab1[which(tab1[,1]%in%lasv1[,3]),]; missing = c()
		tab4c = c(); temp = lasv1[which(!lasv1[,3]%in%tab4b[,1]),3]
		for (j in 1:length(temp)) missing = c(missing,unlist(strsplit(as.character(temp[j]),"_NGA-"))[1])
		for (j in 1:length(missing))
			{
				if (sum(grepl(missing[j],tab1[,1])) == 1)
					{
						id = which(grepl(missing[j],tab1[,1])); tab4c = rbind(tab4c, tab1[id,])
					}	else	{
						print(as.character(temp[j]))
					}
			}
		tab4 = rbind(tab4a,tab4b,tab4c); fasta4 = c()
		for (j in 1:length(fasta1))
			{
				if (grepl(">",fasta1[j]))
					{
						seqID = gsub(">","",fasta1[j])
						if (seqID%in%tab4[,1])
							{
								fasta4 = c(fasta4,fasta1[j],fasta1[j+1])
							}
					}	
			}		
		write(fasta4, paste0("LASV_",segments[i],"_NGA_2019.fasta"))
		write.csv(tab4, paste0("LASV_",segments[i],"_NGA_2019.csv"), quote=F, row.names=F)
	}

# B2. Preliminary phylogenetic analyses

setwd(paste(wd,wdb2,sep="/"))

	# B2.1. Preparing and running analyses in BEAST
	
		# Models used: GTR+G substitution model, relaxed lognormal molecular clock model, constant population size coalescent model
		# Burn-in's: 60000000 states for segment L and 30000000 states for segment S
		# Note for segment S: three sequences were discarded because in an isolated clade
		# 	(MF990888|201600568|Human|TGO|2016-03-12, KU961971|Togo/2016/7082|Human|Togo|2016-02-26, pinneo_NGA-Borno_Lassa_LS_1969)

	# B2.2. Generating alignements and metadata for each clade

cladeSpecificRRWanalyses = FALSE
for (i in 1:length(segments))
	{
		setwd(paste(wd,wdb1,sep="/"))
		fastaAll = scan(paste0("LASV_",segments[i],"_alignment_2.fasta"), what="", sep="\n", quiet=T)
		metadata = read.csv(paste0("LASV_",segments[i],"_alignment_2.csv"), colClasses="character", header=T)
		for (j in 1:length(clades))
			{
				setwd(paste(wd,wdb2,sep="/"))
				tree = read.nexus(paste0("LASV_",segments[i],"_alignment_2_",clades[j],".tree"))
				tips = gsub("'","",tree$tip.label); fasta2 = c(); tab2 = c()
				for (k in 1:length(tips))
					{
						index = which(gsub(" ","",fastaAll)==paste0(">",tips[k]))
						if (length(index) != 1) print(c(i,j,k,"fasta"))
						fasta2 = c(fasta2, fastaAll[index], fastaAll[index+1])
						index = which(gsub(" ","",metadata[,1])==tips[k])
						if (length(index) != 1) print(c(i,j,k,"metadata"))
						tab2 = rbind(tab2, metadata[index,])
					}
				setwd(paste(wd,wdb8,sep="/"))
				write(fasta2, paste0("Alignments_all_sequences/LASV_",segments[i],"_alignment_2_",clades[j],".fasta"))
				write.csv(tab2, paste0("Alignments_all_sequences/LASV_",segments[i],"_alignment_2_",clades[j],".csv"), row.names=F, quote=F)
				tab3 = tab2; tab3 = tab3[which(tab3[,"longitude"]!=""),]; fasta3 = c()
				for (k in 1:length(fasta2))
					{
						if (grepl(">",fasta2[k]))
							{
								seqID = gsub(">","",fasta2[k])
								if (seqID%in%tab3[,1])
									{
										fasta3 = c(fasta3,fasta2[k],fasta2[k+1])
									}
							}
					}
				setwd(paste(wd,wdb4,sep="/"))
				write.csv(tab3, paste0("LASV_",segments[i],"_alignment_3_",clades[j],".csv"), row.names=F, quote=F)
				if (cladeSpecificRRWanalyses == TRUE)
					{
						tab4 = tab3[,c(paste0("segment_",segments[i]),"latitude","longitude")]; colnames(tab4) = c("trait","latitude","longitude")
						write.table(tab4, paste0("LASV_",segments[i],"_alignment_3_",clades[j],".txt"), sep="	", row.names=F, quote=F)
						write(fasta3, paste0("LASV_",segments[i],"_alignment_3_",clades[j],".fasta"))			
					}
			}
		samplingDates = rep(NA,dim(metadata)[1])
		for (j in 1:dim(metadata)[1])
			{
				if (grepl("_",metadata[j,paste0("segment_",segments[i])])) separator = "_"
				if (grepl("\\|",metadata[j,paste0("segment_",segments[i])])) separator = "\\|"
				index = length(unlist(strsplit(metadata[j,paste0("segment_",segments[i])],separator)))
				samplingDate = unlist(strsplit(metadata[j,paste0("segment_",segments[i])],separator))[index]
				if (length(unlist(strsplit(samplingDate,"-"))) == 1) samplingDates[j] = as.numeric(samplingDate)
				if (length(unlist(strsplit(samplingDate,"-"))) == 2) samplingDates[j] = decimal_date(ymd(as.character(paste0(samplingDate,"-15"))))
				if (length(unlist(strsplit(samplingDate,"-"))) == 3) samplingDates[j] = decimal_date(ymd(as.character(samplingDate)))				
			}
		cat("Segment ",segments[i],", most recent sampling date = ",max(samplingDates),sep="","\n") # 2019.178 for both segments
	}

# B3. Temporal signal analyses

	# B3.1. Maximum likelihood trees were inferred using PhyML as implemented in SeaView (with default settings)

	# B3.2. Temporal signal evaluation using regressions of root-to-tip genetic distances against sequence sampling times
	
		# B3.2.1. Estimation of regression R2 using the program TempEst 1.5.3 (cfr. screenshot for the results)
		
		# B3.2.2. Estimation of regression p-values were calculated using the approach of Murray et al. (2016)

setwd(paste(wd,wdb3,sep="/"))
samplingDates_list = list(); c = 0
for (i in 1:length(segments))
	{
		for (j in 1:length(analyses))
			{
				tree = read.tree(paste0("LASV_",segments[i],"_alignment_",analyses[j],".tree"))
				tipLabels = tree$tip.label; samplingDates = rep(NA,length(tree$tip.label))
				for (k in 1:length(tipLabels))
					{
						if (grepl("_",tipLabels[k])) separator = "_"
						if (grepl("\\|",tipLabels[k])) separator = "\\|"
						samplingDate = unlist(strsplit(tipLabels[k],separator))[length(unlist(strsplit(tipLabels[k],separator)))]
						if (length(unlist(strsplit(samplingDate,"-"))) == 1) samplingDates[k] = as.numeric(samplingDate)
						if (length(unlist(strsplit(samplingDate,"-"))) == 2) samplingDates[k] = decimal_date(ymd(as.character(paste0(samplingDate,"-15"))))
						if (length(unlist(strsplit(samplingDate,"-"))) == 3) samplingDates[k] = decimal_date(ymd(as.character(samplingDate)))				
					}
				c = c+1; samplingDates_list[[c]] = samplingDates
			}
	}
source("Temporal_signal_functions/mantelCounding.r")
source("Temporal_signal_functions/randRegression.r")
source("Temporal_signal_functions/tempSignalFunctions.r")
pValues1 = c(); c = 0
for (i in 1:length(segments))
	{
		for (j in 1:length(analyses))
			{
				tree = read.tree(paste0("LASV_",segments[i],"_alignment_",analyses[j],".tree")); c = c+1; tipDates = samplingDates_list[[c]]
				test = pathogen.permutation.test(phy=tree, dates=tipDates, use.clusters=F, auto.round.dates=F, nreps=1000)
					#  (approach of Murray et al. based on 1,000 random permutations of the sequence sampling dates)
				pValues = c(pValues, test$p_value)
			}	# p-values: < 0.001 (L, alignment 2), 0.006 (L, alignment 3), < 0.001 (S, alignment 2), < 0.001 (S, alignment 3)
	}

# B4. Continuous phylogeographic analyses

setwd(paste(wd,wdb4,sep="/"))
nberOfExtractionFiles = 1000

	# 4.1. Preparing and running analyses in BEAST
	
		# Models used: GTR+G substitution model, relaxed lognormal molecular clock model, skygrid coalescent model (100 grid points, cut-off = 1500), RRW (gamma, jitter = 0.01)
		# N.B.: productStatistics "treeLengthPrecision" (1, 2, and 3) have to be manually added and logged in the XML files

	# 4.2. Extraction of 1,000 tree files

options(scipen=9)
mostRecentSamplingDatum = 2019.178
burnIns = c(100000000,250000000)

		# 4.2.1. Selecting 1,000 trees from the post-burn-in distribution

for (i in 1:length(segments))
	{
		allTrees = scan(file=paste0("LASV_",segments[i],"_align_3_all_RRW.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		lineIndex1 = which(grepl(paste0("tree STATE_0 "),allTrees))
		lineIndex2 = which(grepl(paste0("tree STATE_",burnIns[i]," "),allTrees))
		lineIndex3 = which(grepl(paste0("End;"),allTrees))[2]
		samplingFrequency = floor(((lineIndex3-1)-(lineIndex2))/nberOfExtractionFiles)
		treesToSample = seq(((lineIndex3-1)-((nberOfExtractionFiles-1)*samplingFrequency)),(lineIndex3-1),samplingFrequency)
		selectedTrees = c(allTrees[1:(lineIndex1-1)],allTrees[treesToSample],"End;")
		write(selectedTrees, paste0("LASV_",segments[i],"_align_3_all_",nberOfExtractionFiles,".trees"))
		log = scan(file=paste0("LASV_",segments[i],"_align_3_all_RRW.log"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		logStates = rep(NA, length(log))
		for (j in 1:length(logStates))
			{
				logStates[j] = unlist(strsplit(log[j],"\t"))[1]
			}
		selectedTrees = scan(file=paste0("LASV_",segments[i],"_align_3_all_1000.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		lines = which(grepl("tree STATE_",selectedTrees)); selectedTreeIDs = rep(NA, length(lines))
		for (j in 1:length(selectedTreeIDs))
			{
				selectedTreeIDs[j] = gsub("STATE_","",unlist(strsplit(selectedTrees[lines[j]]," "))[2])
			}
		selectedLogStates = logStates[which(logStates%in%selectedTreeIDs)]
		indices = c(1:5,which(logStates%in%selectedLogStates))
		write(log[indices], paste0("LASV_",segments[i],"_align_3_all_",nberOfExtractionFiles,".log"))		
	}

		# 4.2.2. Extractions of the spatio-temporal information contained in posterior trees

for (i in 1:length(segments))
	{	
		localTreesDirectory = paste0("LASV_",segments[i],"_align_3_all_1_extractions")
		allTrees = scan(file=paste0("LASV_",segments[i],"_align_3_all_",nberOfExtractionFiles,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		burnIn = 0; randomSampling = FALSE
		nberOfTreesToSample = nberOfExtractionFiles
		coordinateAttributeName = "location"; nberOfCores = 10
		treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
	}

		# 4.2.3. Generating the maximum clade credibility (MCC) trees with TreeAnnotator

for (i in 1:length(segments))
	{
		input = paste0("LASV_",segments[i],"_align_3_all_",nberOfExtractionFiles,".trees")
		output = paste0("LASV_",segments[i],"_align_3_all_MCC.tree")
		system(paste0("BEAST_ver1104_stable_release/bin/treeannotator -burninTrees 0 -heights keep ",input," ",output))
	}

		# 4.2.4. Extractions of the spatio-temporal information contained in MCC trees

source("mccTreeExtractionFunction.r")
for (i in 1:length(segments))
	{
		mcc_tre = readAnnotatedNexus(paste0("LASV_",segments[i],"_align_3_all_MCC.tree"))
		mcc_tab = mccTreeExtractionFunction(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc_tab, paste0("LASV_",segments[i],"_align_3_all_MCC.csv"), row.names=F, quote=F)
	}	# to do: open and export again the MCC trees in FigTree (Nexus format)

		# 4.2.5. Dividing extraction files per lineage (clade)

source("isolateACladeFromTipNodes.r")
for (i in 1:length(segments))
	{
		localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_all_1_extractions")
		metadata_clades = list()
		for (j in 1:length(clades))
			{
				metadata_clades[[j]] = read.csv(paste0("LASV_",segments[i],"_alignment_3_",clades[j],".csv"), header=T)
				metadata_clades[[j]][,1] = gsub(" ","",metadata_clades[[j]][,1])
			}
		mcc_tre = readAnnotatedNexus(paste0("LASV_",segments[i],"_align_3_all_MCC.tree"))
		tipNodes = mcc_tre$edge[which(!mcc_tre$edge[,2]%in%mcc_tre$edge[,1]),2]
		tipClades = rep(NA, length(tipNodes)); samplingCoordinates = rep(NA, length(tipNodes))
		mcc_tre$tip.label = gsub(" ","",gsub("'","",mcc_tre$tip.label))
		for (j in 1:length(tipNodes))
			{
				indices = c(); index = NA
				for (k in 1:length(clades))
					{
						if (mcc_tre$tip.label[tipNodes[j]]%in%metadata_clades[[k]][,paste0("segment_",segments[i])])
							{
								indices = c(indices, k)
							}
					}
				lon = round(mcc_tre$annotations[[which(mcc_tre$edge[,2]==tipNodes[j])]]$location2,7)
				lat = round(mcc_tre$annotations[[which(mcc_tre$edge[,2]==tipNodes[j])]]$location1,7)
				if (length(indices) != 1)
					{
						tipClades[j] = NA; samplingCoordinates[j] = paste0(lon,"_",lat)
					}	else	{
						index = indices
						tipClades[j] = clades[index]; samplingCoordinates[j] = paste0(lon,"_",lat)
					}
			}
		metadata = cbind(mcc_tre$tip.label,samplingCoordinates,tipClades)
		for (j in 1:length(clades))
			{
				localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_",clades[j],"_extractions")
				dir.create(file.path(localTreesDirectory2), showWarnings=F)
				for (k in 0:nberOfExtractionFiles)
					{
						if (k == 0)
							{
								tab1 = read.csv(paste0("LASV_",segments[i],"_align_3_all_MCC.csv"), header=T)
							}	else	{
								tab1 = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",k,".csv"), header=T)
							}
						samplingCoordinates = rep(NA, dim(tab1)[1]); node2s = rep(NA, dim(tab1)[1])
						for (l in 1:dim(tab1)[1])
							{
								lon = round(tab1[l,"endLon"],7); lat = round(tab1[l,"endLat"],7)
								samplingCoordinates[l] = paste0(lon,"_",lat); node2s[l] = tab1[l,"node2"]
							}
						indices = which(!tab1[,"node2"]%in%tab1[,"node1"])
						samplingCoordinates = samplingCoordinates[indices]; node2s = node2s[indices]
						indices = which(metadata[,"tipClades"]==clades[j]); tipNodes = rep(NA, length(indices))
						for (l in 1:length(indices))
							{
								index = which(samplingCoordinates == metadata[indices[l],"samplingCoordinates"])
								if (length(index) != 1) cat(i,j,k,l,"\n")
								tipNodes[l] = node2s[index]
							}
						tab2 = isolateACladeFromTipNodes(tab1, tipNodes)
						if (k == 0)
							{
								write.csv(tab2, paste0("LASV_",segments[i],"_align_3_",clades[j],"_MCC.csv"), row.names=F, quote=F)
							}	else	{
								write.csv(tab2, paste0(localTreesDirectory2,"/TreeExtractions_",k,".csv"), row.names=F, quote=F)
							}
					}
			}
	}
	
		# 4.2.6. Gathering all the extraction files per study area

for (i in 1:length(segments))
	{
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_NGA_extractions")
		dir.create(file.path(localTreesDirectory2), showWarnings=F)
		for (j in 0:nberOfExtractionFiles)
			{
				tab = c()
				for (k in 2:length(clades))
					{
						if (j == 0) tab = rbind(tab, read.csv(paste0("LASV_",segments[i],"_align_3_",clades[k],"_MCC.csv"), head=T))
						localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_",clades[k],"_extractions")
						if (j != 0) tab = rbind(tab, read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",j,".csv"), head=T))
					}
				if (j == 0) write.csv(tab, paste0("LASV_",segments[i],"_align_3_NGA_MCC.csv"), row.names=F, quote=F)
				if (j != 0) write.csv(tab, paste0(localTreesDirectory2,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
			}
	}
	
		# 4.2.7. Gathering all the extraction files per segment

for (i in 1:length(segments))
	{
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_all_2_extractions")
		dir.create(file.path(localTreesDirectory2), showWarnings=F)
		for (j in 1:nberOfExtractionFiles)
			{
				tab = c()
				for (k in 1:length(clades))
					{
						localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_",clades[k],"_extractions")
						tab = rbind(tab, read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",j,".csv"), head=T))
					}
				write.csv(tab, paste0(localTreesDirectory2,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
			}
	}
	
		# 4.2.8. Generating extraction tables only gathering tip branches

for (i in 1:length(segments))
	{
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_tips2_extractions")
		dir.create(file.path(localTreesDirectory2), showWarnings=F)
		for (j in 1:nberOfExtractionFiles)
			{
				localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_all_2_extractions")
				tab = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",j,".csv"), head=T)
				tab = tab[which(!tab[,"node2"]%in%tab[,"node1"]),]
				write.csv(tab, paste0(localTreesDirectory2,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
			}
	}
	
		# 4.2.9. Regurlarly subsampling 100 out of 1,000 tree extractions

for (i in 1:length(segments))
	{
		localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_all_2_extractions")
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_all_2_extract_100")
		dir.create(file.path(localTreesDirectory2), showWarnings=F); n = 0
		for (j in seq(10,nberOfExtractionFiles,10))
			{
				tab = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",j,".csv"), head=T); n = n+1
				write.csv(tab, paste0(localTreesDirectory2,"/TreeExtractions_",n,".csv"), row.names=F, quote=F)
			}
		localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_tips2_extractions")
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_tips2_extract_100")
		dir.create(file.path(localTreesDirectory2), showWarnings=F); n = 0
		for (j in seq(10,nberOfExtractionFiles,10))
			{
				tab = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",j,".csv"), head=T); n = n+1
				write.csv(tab, paste0(localTreesDirectory2,"/TreeExtractions_",n,".csv"), row.names=F, quote=F)
			}
	}

	# 4.3. Generating and saving annual progression polygons (convex hull polygons)

for (i in 1:length(segments))
	{
		for (j in 1:length(study_areas))
			{
				points = c(); years = c(1900:2019)
				localTreesDirectory = paste0("LASV_",segments[i],"_align_3_",study_areas[j],"_extractions")
				for (k in 1:nberOfExtractionFiles)
					{
						tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",k,".csv"), header=T)
						points = rbind(points, tab[,c("endLon","endLat","endYear")])
					}
				colnames(points) = c("lon","lat","time")
				for (y in 1:length(years))
					{
						selectedNodes = points[which(floor(points[,"time"])<=years[y]+1),]
						if (dim(selectedNodes)[1] > 0)
							{
								selectedNodes = unique(selectedNodes[,c("lon","lat")]); hull = chull(selectedNodes); hull = c(hull,hull[1])
								p = Polygon(selectedNodes[hull,]); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								contourPolygons_df = SpatialPolygonsDataFrame(sps, data.frame(ID=1:length(sps))); # plot(sps, lwd=0.05, add=T)
								polygonName = paste0("LASV_segment_",segments[i],"_",study_areas[j],"_",years[y])
								writeOGR(contourPolygons_df, dsn="LASV_minimum_convex_hulls", layer=polygonName, driver="ESRI Shapefile", overwrite_layer=T)
							}
					}
			}
	}
	
	# 4.4. Estimating kernel density polygons

percentage = 80; precision = 20	
for (i in 1:length(segments))
	{
		for (j in 1:length(clades))
			{
				localTreesDirectory = paste0("LASV_",segments[i],"_align_3_",clades[j],"_extractions")
				mcc = read.csv(paste0("LASV_",segments[i],"_align_3_",clades[j],"_MCC.csv"))
				startDatum = min(mcc[,"startYear"]); prob = percentage/100
				polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
				for (k in 1:length(polygons))
					{
						year = as.numeric(names(polygons[[k]]))
						layerName = paste0("LASV_",segments[i],"_align_3_",clades[j],"_",percentage,"_",precision,"_",year)
						writeOGR(polygons[[k]], dsn=paste0("LASV_80HPD_polygons_20yrs"), layer=layerName, driver="ESRI Shapefile")
					}
			}
	}

	# 4.5. Plotting the global dispersal history graphs

setwd(wd)
background1 = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/Elevation_MRU_0.008.asc"))
background1 = crop(background1, e_MRU); background1[!is.na(background1[])] = 1; r = background1
background_cols1 = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
borders1 = crop(readOGR(dsn=paste0(wdb6,"/Environmental_files/Natural_Earth/"), layer="International_borders"), e_MRU)
lakes1 = crop(readOGR(dsn=paste0(wdb6,"/Environmental_files/Natural_Earth/"), layer="Natural_Earth_lakes"), e_MRU)
rivers1 = crop(readOGR(dsn=paste0(wdb6,"/Environmental_files/LASV_rivers/"), layer="FAO_GeoNetwork"), e_MRU)
rivers1@data[which(rivers1@data[,"Strahler"]<3),"Strahler"] = 0
riverCols1 = rep("#AFBED5",dim(rivers1@data)[1]); riverCols1[which(rivers1@data[,"Strahler"]<3)] = NA
background2 = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/Elevation_NGA_0.008.asc"))
background2 = crop(background2, e_NGA); background2[!is.na(background2[])] = 1; r = background2
background_cols2 = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
borders2 = crop(readOGR(dsn=paste0(wdb6,"/Environmental_files/Natural_Earth/"), layer="International_borders"), e_NGA)
lakes2 = crop(readOGR(dsn=paste0(wdb6,"/Environmental_files/Natural_Earth/"), layer="Natural_Earth_lakes"), e_NGA)
rivers2 = crop(readOGR(dsn=paste0(wdb6,"/Environmental_files/LASV_rivers/"), layer="FAO_GeoNetwork"), e_NGA)
rivers2@data[which(rivers2@data[,"Strahler"]<3),"Strahler"] = 0
riverCols2 = rep("#AFBED5",dim(rivers2@data)[1]); riverCols2[which(rivers2@data[,"Strahler"]<3)] = NA
titles = c("MRU lineage","NGA lineage IV","NGA lineage II","NGA lineage III")
cols = colorRampPalette(brewer.pal(11,'RdYlGn'))(131)[21:121]
for (i in 1:length(segments))
	{
		pdf(paste0("Dispersal_graphs_seg",segments[i],".pdf"), width=12, height=10.2)
		par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(1.5,2.7,1.0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		for (j in 1:length(clades))
			{
				if (clades[j] == "MRU")
					{
						background = background1; background_cols = background_cols1
						borders = borders1; lakes = lakes1; rivers = rivers1; riverCols = riverCols1
					}	else	{
						background = background2; background_cols = background_cols2
						borders = borders2; lakes = lakes2; rivers = rivers2; riverCols = riverCols2
					}
				mcc = read.csv(paste0(wdb4,"/LASV_",segments[i],"_align_3_",clades[j],"_MCC.csv"))
				mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
				mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
				minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
				ancestralC = (((mcc[1,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				startYearsC = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				endYearsC = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
				col_ancestral = cols[ancestralC]; cols_endYears = cols[endYearsC]
				files = list.files(paste0(wdb4,"/LASV_80HPD_polygons_20yrs")); pols = list(); cols_pol = list()
				files = files[which((grepl(paste0("_",segments[i],"_"),files))&(grepl(clades[j],files))&(grepl(".shp",files)))]
				for (k in 1:length(files))
					{
						pol = shapefile(paste0(wdb4,"/LASV_80HPD_polygons_20yrs/",files[k]))
						pol@proj4string = CRS("+init=epsg:4326"); pols[[k]] = pol
						date = gsub(".shp","",unlist(strsplit(files[k],"_"))[length(unlist(strsplit(files[k],"_")))])
						precision = 20; date = as.numeric(date)+1-(precision/2)
						yearM = (((date-minYear)/(maxYear-minYear))*100)+1
						cols_pol[[k]] = cols[yearM]
					}
				for (k in 1:length(cols_pol)) cols_pol[[k]] = paste0(cols_pol[[k]],"70")
				plot(background, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="white", colNA="gray90")
				for (k in 1:length(pols)) plot(pols[[k]], axes=F, col=cols_pol[[k]], add=T, border=NA)
				plot(lakes, add=T, lwd=0.7, col=unique(riverCols)[!is.na(unique(riverCols))], border=NA)
				plot(rivers, add=T, lwd=rivers@data[,"Strahler"]/3.5, col=riverCols)
				if (j == 2)
					{
						mtext(expression(bold(Nigeria)), side=3, line=-8.8, at=6.85, cex=0.9, font=1, col="gray30")
					}
				mtext(titles[j], side=3, line=-27.3, at=3.7, cex=0.9, font=1, col="gray30")
				plot(borders, add=T, lwd=2, col="gray60", lty=2)
				rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.2)
				axis(1,c(ceiling(xmin(background)), floor(xmax(background))), pos=ymin(background), mgp=c(0,0.7,0), 
					 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
				axis(2, c(ceiling(ymin(background)), floor(ymax(background))), pos=xmin(background), mgp=c(0,0.8,0),
					 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
				for (k in 1:dim(mcc)[1])
					{
						curvedarrow(cbind(mcc[k,"startLon"],mcc[k,"startLat"]), cbind(mcc[k,"endLon"],mcc[k,"endLat"]), arr.length=0,
						    		arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
					}
				for (k in dim(mcc)[1]:1)
					{
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=16, col=cols_endYears[k], cex=1.0)
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=1, col="gray30", cex=1.0, lwd=0.1)
						if (k == 1)
							{
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=16, col=col_ancestral, cex=1.0)
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=1, col="gray30", cex=1.0, lwd=0.1)
							}
					}
				legendRast = raster(as.matrix(c(minYear,maxYear)))
				plot(legendRast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.810,0.825,0.09,0.56), alpha=1.0,
					 legend.args=list(text="", cex=0.9, line=0.5, col="gray30"), axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.2, tck=-0.8, 
					 col.axis="gray30", line=0, mgp=c(0,0.6,0)))
				if (j == 1)
					{
						legend(x=-14.5, y=6.3, legend=c("S = 3 ","S = 4 ","S = 5 ","S = 6 ","S = 7 "), lwd=c(3:7)/3.5, cex=0.9, lty=1,
							   col=unique(riverCols)[!is.na(unique(riverCols))], x.intersp=0.5, border=NA, bg="white", box.lwd=0.2, box.col="gray30")
			   		}
			}
		dev.off()
	}
titles = c("MRU lineages","NGA lineages")
cols = list(); cols1 = list(); cols2 = list()
cols[[1]] = rev(colorRampPalette(brewer.pal(9,'YlGn'))(161)[11:111])
cols[[2]] = colorRampPalette(brewer.pal(9,'Reds'))(151)[11:111]
cols[[3]] = colorRampPalette(brewer.pal(9,'Blues'))(151)[11:111]
cols[[4]] = colorRampPalette(brewer.pal(9,'Oranges'))(191)[11:111]
cols1[[1]] = rgb(222,67,39,255,maxColorValue=255); cols2[[1]] = rgb(222,67,39,100,maxColorValue=255) # red
cols1[[2]] = rgb(70,118,187,255,maxColorValue=255); cols2[[2]] = rgb(70,118,187,100,maxColorValue=255) # blue
cols1[[3]] = rgb(250,165,33,255,maxColorValue=255); cols2[[3]] = rgb(250,165,33,100,maxColorValue=255) # orange
pdf(paste0("Dispersal_graphs_all_NEW.pdf"), width=12, height=10.2); # dev.new(width=12, height=10.2)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(1.5,2.7,1.0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 1:length(segments))
	{	
		for (j in 1:length(clades))
			{
				if (clades[j] == "MRU")
					{
						background = background1; background_cols = background_cols1
						borders = borders1; lakes = lakes1; rivers = rivers1; riverCols = riverCols1
					}	else	{
						background = background2; background_cols = background_cols2
						borders = borders2; lakes = lakes2; rivers = rivers2; riverCols = riverCols2
					}
				mcc = read.csv(paste0(wdb4,"/LASV_",segments[i],"_align_3_",clades[j],"_MCC.csv"))
				mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
				mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
				if (clades[j] == "MRU")
					{
						minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
					}
				if (clades[j] == "NIG1")
					{
						minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
						for (k in 3:length(clades))
							{
								mcc = read.csv(paste0(wdb4,"/LASV_",segments[i],"_align_3_",clades[k],"_MCC.csv"))
								mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
								mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
								if (minYear > min(mcc[,"startYear"])) minYear = min(mcc[,"startYear"])
								if (maxYear < min(mcc[,"endYear"])) maxYear = max(mcc[,"endYear"])
							}
						mcc = read.csv(paste0(wdb4,"/LASV_",segments[i],"_align_3_",clades[j],"_MCC.csv"))
						mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
						mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
					}
				ancestralC = (((mcc[1,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				startYearsC = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				endYearsC = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
				col_ancestral = cols[[j]][ancestralC]; cols_endYears = cols[[j]][endYearsC]
				files = list.files(paste0(wdb4,"/LASV_80HPD_polygons_20yrs")); pols = list(); cols_pol = list()
				files = files[which((grepl(paste0("_",segments[i],"_"),files))&(grepl(clades[j],files))&(grepl(".shp",files)))]
				for (k in 1:length(files))
					{
						pol = shapefile(paste0(wdb4,"/LASV_80HPD_polygons_20yrs/",files[k]))
						pol@proj4string = CRS("+init=epsg:4326"); pols[[k]] = pol
						date = gsub(".shp","",unlist(strsplit(files[k],"_"))[length(unlist(strsplit(files[k],"_")))])
						precision = 20; date = as.numeric(date)+1-(precision/2)
						yearM = (((date-minYear)/(maxYear-minYear))*100)+1
						cols_pol[[k]] = cols[[j]][yearM]
					}
				for (k in 1:length(cols_pol)) cols_pol[[k]] = paste0(cols_pol[[k]],"70")
				if ((clades[j] == "MRU")|(clades[j] == "NIG1"))
					{
						plot(background, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="white", colNA="gray90")
					 }
				for (k in 1:length(pols)) plot(pols[[k]], axes=F, col=cols_pol[[k]], add=T, border=NA)
				if (clades[j] == "NIG1") mtext(paste0("SEGMENT ",segments[i]), side=3, line=-24.0, at=3.55, cex=0.8, font=1, col="gray30")
				if (clades[j] == "NIG1") mtext("NGA lineage IV", side=3, line=-25.0, at=3.65, cex=0.8, font=1, col=cols1[[1]])
				if (clades[j] == "NIG2") mtext("NGA lineage II", side=3, line=-26.0, at=3.65, cex=0.8, font=1, col=cols1[[2]])
				if (clades[j] == "NIG3") mtext("NGA lineage III", side=3, line=-27.0, at=3.65, cex=0.8, font=1, col=cols1[[3]])
				if ((clades[j] == "MRU")|(clades[j] == "NIG1"))
					{
						plot(lakes, add=T, lwd=0.7, col=unique(riverCols)[!is.na(unique(riverCols))], border=NA)
						plot(rivers, add=T, lwd=rivers@data[,"Strahler"]/3.5, col=riverCols)
						plot(borders, add=T, lwd=2, col="gray60", lty=2)
						rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.2)
						axis(1,c(ceiling(xmin(background)), floor(xmax(background))), pos=ymin(background), mgp=c(0,0.7,0), 
							 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
						axis(2, c(ceiling(ymin(background)), floor(ymax(background))), pos=xmin(background), mgp=c(0,0.8,0),
							 cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
					}	
				if ((clades[j] == "MRU")|(clades[j] == "NIG3"))
					{	
						if (clades[j] == "MRU")
							{
								mtext(paste0("SEGMENT ",segments[i]), side=3, line=-26, at=-11.7, cex=0.8, font=1, col="gray30")
								mtext(titles[j], side=3, line=-27.0, at=-11.65, cex=0.8, font=1, col="gray30")
							}
						if (clades[j] == "NIG3")
							{
								mtext(expression(bold(Nigeria)), side=3, line=-8.8, at=6.85, cex=0.9, font=1, col="gray30")
							}
					 }
				for (k in 1:dim(mcc)[1])
					{
						curvedarrow(cbind(mcc[k,"startLon"],mcc[k,"startLat"]), cbind(mcc[k,"endLon"],mcc[k,"endLat"]), arr.length=0,
						    		arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
					}
				for (k in dim(mcc)[1]:1)
					{
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=16, col=cols_endYears[k], cex=1.0)
						points(mcc[k,"endLon"], mcc[k,"endLat"], pch=1, col="gray30", cex=1.0, lwd=0.1)
						if (k == 1)
							{
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=16, col=col_ancestral, cex=1.0)
								points(mcc[k,"startLon"], mcc[k,"startLat"], pch=1, col="gray30", cex=1.0, lwd=0.1)
							}
					}
				legendRast = raster(as.matrix(c(minYear,maxYear)))
				if ((clades[j] == "MRU")|(clades[j] == "NIG3"))
					{
						plot(legendRast, legend.only=T, add=T, col=cols[[j]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.810,0.825,0.09,0.56), alpha=1.0,
					 		 legend.args=list(text="", cex=0.9, line=0.5, col="gray30"), axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.2, tck=-0.8, 
					 		 col.axis="gray30", line=0, mgp=c(0,0.6,0)))
					 }
				if (clades[j] == "NIG2")
					{
						plot(legendRast, legend.only=T, add=T, col=cols[[j]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.795,0.810,0.09,0.56), alpha=1.0,
					 		 legend.args=list(text="", cex=0.9, line=0.5, col="gray30"), axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.2, tck=0, 
					 		 col.axis="gray30", line=0, mgp=c(0,0.6,0), at=c(2000), labels=c("")))					 	
					 }
				if (clades[j] == "NIG1")
					{
						plot(legendRast, legend.only=T, add=T, col=cols[[j]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.780,0.795,0.09,0.56), alpha=1.0,
					 		 legend.args=list(text="", cex=0.9, line=0.5, col="gray30"), axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.2, tck=-0.8, 
					 		 col.axis="gray30", line=0, mgp=c(0,0.6,0), at=c(2000), labels=c("")))					 	
					 }
				if (clades[j] == "MRU")
					{
						legend(x=-14.5, y=6.3, legend=c("S = 3 ","S = 4 ","S = 5 ","S = 6 ","S = 7 "), lwd=c(3:7)/3.5, cex=0.9, lty=1,
							   col=unique(riverCols)[!is.na(unique(riverCols))], x.intersp=0.5, border=NA, bg="white", box.lwd=0.2, box.col="gray30")
			   		}
			}
	}
dev.off()

# B5. Estimating dispersal statistics

setwd(wd)
nberOfExtractionFiles = 1000
timSlices = 200; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 10; slidingWindow = 5
for (i in 1:length(segments))
	{
		for (j in 1:length(study_areas))
			{
				localTreesDirectory = paste0(wdb4,"/LASV_",segments[i],"_align_3_",study_areas[j],"_extractions")
				outputName = paste0(wdb4,"/LASV_",segments[i],"_align_3_",study_areas[j])
				spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
			}
	}
for (i in 1:length(segments))
	{
		for (j in 2:length(clades))
			{
				localTreesDirectory = paste0(wdb4,"/LASV_",segments[i],"_align_3_",clades[j],"_extractions")
				outputName = paste0(wdb4,"/LASV_",segments[i],"_align_3_",clades[j])
				spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
			}
	}

analyses = c("LASV_L_align_3_MRU","LASV_L_align_3_NGA","LASV_L_align_3_NIG1","LASV_L_align_3_NIG2","LASV_L_align_3_NIG3",
			 "LASV_S_align_3_MRU","LASV_S_align_3_NGA","LASV_S_align_3_NIG1","LASV_S_align_3_NIG2","LASV_S_align_3_NIG3")
statistics = c("mean_branch_dispersal_velocity","weighted_branch_dispersal_velocity")
results = matrix(nrow=length(analyses), ncol=length(statistics)); colNames = c(); rounds = c(2,1)
for (i in 1:length(analyses))
	{
		tab = read.table(paste0(wdb5,"/",analyses[i],"_estimated_dispersal_statistics.txt"), header=T)
		for (j in 1:length(statistics))
			{
				median = round(median(tab[,statistics[j]]),rounds[j])
				qtiles = round(quantile(tab[,statistics[j]],c(0.025,0.975)),rounds[j])
				results[i,j] = as.character(paste0(median," [",qtiles[1],"-",qtiles[2],"]"))
			}
	}
colNames = c("Mean branch dispersal velocity (km/year)","Weighted branch dispersal velocity (km/year)")
colnames(results) = colNames; row.names(results) = analyses
write.table(results, "LASV_dispersal_statistics_NEW.txt", sep="	", quote=F)

pdf("Dispersal_velocities_NEW.pdf", width=8.5, height=2.7); # dev.new(width=8.5, height=2.7)
par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(1,1.3,0,0.75), mar=c(1.5,1.2,1,0)); cols1 = list(); cols2 = list()
lineage_names = c("MRU lineage","NGA lineage II","NGA lineage III")
analyses = c("LASV_L_align_3_MRU","LASV_L_align_3_NIG2","LASV_L_align_3_NIG3",
			 "LASV_S_align_3_MRU","LASV_S_align_3_NIG2","LASV_S_align_3_NIG3")
cols1[[1]] = rgb(69,139,0,255,maxColorValue=255); cols2[[1]] = rgb(69,139,0,100,maxColorValue=255) # green
cols1[[2]] = rgb(70,118,187,255,maxColorValue=255); cols2[[2]] = rgb(70,118,187,100,maxColorValue=255) # blue
cols1[[3]] = rgb(250,165,33,255,maxColorValue=255); cols2[[3]] = rgb(250,165,33,100,maxColorValue=255) # orange
cols1[[4]] = cols1[[1]]; cols1[[5]] = cols1[[2]]; cols1[[6]] = cols1[[3]]
cols2[[4]] = cols2[[1]]; cols2[[5]] = cols2[[2]]; cols2[[6]] = cols2[[3]]
for (i in 1:length(analyses))
	{
		segment = unlist(strsplit(analyses[i],"_"))[2]
		tab = read.table(paste0(wdb5,"/",analyses[i],"_estimated_dispersal_statistics.txt"), header=T)[,statistics[s]]
		if (i%in%c(1,4)) plot(density(tab), type="l", axes=F, ann=F, col=NA, xlim=c(0.3,1.8), ylim=c(0,10.9))
		polygon(density(tab), lwd=1, col=cols2[[i]], border=NA); lines(density(tab), lwd=0.7, col=cols1[[i]])
		if (i%in%c(1,4)) 
			{
				axis(side=1, pos=0.0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0, tck=-0.020, col.axis="gray30")
				axis(side=2, pos=0.3, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.1,0), lwd=0, tck=-0.015, col.axis="gray30") 
				title(xlab="Weighted lineage dispersal velocity (km/year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
				rect(0.3, 0, 1.8, 10.9, lwd=0.2, border="gray30"); # box(lwd=0.2, col="gray30")
				mtext(side=3, paste0("SEGMENT ",segment), line=-2.2, at=1.52, cex=0.75, col="gray30")
			}
		if (i == 1)
			{
				title(ylab="Density", cex.lab=0.65, mgp=c(0.5,0,0), col.lab="gray30")
				legend(1.34, 9, legend=lineage_names, pch=16, col=unlist(cols2)[1:3], border="gray30", text.col="gray30", 
					   bty="n", pt.cex=1, cex=0.6, x.intersp=0.9, y.intersp=1.2)
			}
	}
dev.off()

pdf("Spatial_wavefronts_NEW.pdf", width=8.5, height=2.7); # dev.new(width=8.5, height=2.7)
par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(1,1.7,0,0.75), mar=c(1.5,1.2,1,0)); cols1 = list(); cols2 = list()
lineage_names = c("MRU lineage","NGA lineage II","NGA lineage III")
analyses = c("LASV_L_align_3_MRU","LASV_L_align_3_NIG2","LASV_L_align_3_NIG3",
			 "LASV_S_align_3_MRU","LASV_S_align_3_NIG2","LASV_S_align_3_NIG3")
cols1[[1]] = rgb(69,139,0,255,maxColorValue=255); cols2[[1]] = rgb(69,139,0,100,maxColorValue=255) # green
cols1[[2]] = rgb(70,118,187,255,maxColorValue=255); cols2[[2]] = rgb(70,118,187,100,maxColorValue=255) # blue
cols1[[3]] = rgb(250,165,33,255,maxColorValue=255); cols2[[3]] = rgb(250,165,33,100,maxColorValue=255) # orange
cols1[[4]] = cols1[[1]]; cols1[[5]] = cols1[[2]]; cols1[[6]] = cols1[[3]]
cols2[[4]] = cols2[[1]]; cols2[[5]] = cols2[[2]]; cols2[[6]] = cols2[[3]]
for (i in 1:length(analyses))
	{
		segment = unlist(strsplit(analyses[i],"_"))[2]
		tab1 = read.table(paste0(wdb5,"/",analyses[i],"_median_spatial_wavefront_distance.txt"), header=T)
		tab2 = read.table(paste0(wdb5,"/",analyses[i],"_95%HPD_spatial_wavefront_distance.txt"), header=T)
		if (i%in%c(1,4)) plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,750), xlim=c(1300,2018), col=NA)
		slicedTimes = tab1[,1]; waveFrontDistances1MedianValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
		xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=cols2[[i]], border=0)
		lines(slicedTimes, waveFrontDistances1MedianValue, lwd=1, col=cols1[[i]], xlim=c(1400,2018))
		if (i%in%c(1,4)) 
			{
				axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.05,0), lwd=0, tck=-0.020, col.axis="gray30")
				axis(side=2, pos=1300, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.1,0), lwd=0, tck=-0.015, col.axis="gray30", 
				     at=seq(0,600,200), labels=c("0","200","400","      600 km"))
				title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
				rect(1300, 0, 2018, 750, lwd=0.2, border="gray30"); # box(lwd=0.2, col="gray30")
				mtext(side=3, paste0("SEGMENT ",segment), line=-2.0, at=1420, cex=0.75, col="gray30")
				legend(1337, 630, legend=lineage_names, pch=16, col=unlist(cols2)[1:3], border="gray30", text.col="gray30", 
					   bty="n", pt.cex=1, cex=0.6, x.intersp=0.9, y.intersp=1.2)
			}
		if (i == 1) title(ylab="Furthest extent of epidemic wavefront", cex.lab=0.65, mgp=c(0.5,0,0), col.lab="gray30")
	}
dev.off()

# B6. Post hoc landscape phylogeography

setwd(wd)
nberOfExtractionFiles = 1000

	# 6.1. Generating null models of lineages dispersal

		# 6.1.1. Randomisation of tree branches position (not used)

analyses = c("LASV_L_align_3_MRU","LASV_L_align_3_NIG1","LASV_L_align_3_NIG2","LASV_L_align_3_NIG3",
			 "LASV_S_align_3_MRU","LASV_S_align_3_NIG1","LASV_S_align_3_NIG2","LASV_S_align_3_NIG3")
for (i in 1:length(analyses))
	{
		studyArea = unlist(strsplit(analyses[i],"_"))[length(unlist(strsplit(analyses[i],"_")))]
		if (studyArea != "MRU") studyArea = "NGA"
		rast = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/Elevation_",studyArea,"_0.008.asc")); rast[!is.na(rast[])] = 0
		localTreesDirectory = paste0(wdb4,"/",analyses[i],"_extractions")
		envVariables = list(rast); randomProcedure = 3; nberOfCores = 1; showingPlots = FALSE
		treesRandomisation(localTreesDirectory, nberOfExtractionFiles, envVariables, randomProcedure, nberOfCores, showingPlots)
	}

		# 6.1.2. Relaxed random walk simulations along posterior trees

analyses = c("LASV_L_align_3_MRU","LASV_L_align_3_NIG1","LASV_L_align_3_NIG2","LASV_L_align_3_NIG3",
			 "LASV_S_align_3_MRU","LASV_S_align_3_NIG1","LASV_S_align_3_NIG2","LASV_S_align_3_NIG3")
backgroundRasters_list = list()
for (i in 1:length(segments))
	{
		backgroundRasters = list()
		for (j in 1:length(study_areas))
			{
				studyArea = study_areas[j]; points1 = c(); points2 = c()
				if (studyArea != "MRU") studyArea = "NIG"
				indices = which(grepl(paste0("_",segments[i],"_"),analyses)&grepl(studyArea,analyses))
				for (k in 1:length(indices))
					{
						localTreesDirectory = paste0(wdb4,"/",analyses[indices[k]],"_extractions")
						for (l in 1:nberOfExtractionFiles)
							{
								tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",l,".csv"), header=T)
								points1 = rbind(points1, tab[,c("startLon","startLat")])
								points2 = rbind(points2, tab[,c("endLon","endLat")])
							}
					}
				colnames(points1) = c("lon","lat"); colnames(points2) = c("lon","lat")
				points = rbind(points1, points2); jitter = 0
				if (jitter > 0)
					{
						points1 = points; points2 = points; points3 = points; points4 = points
						points1[,1] = points1[,1]-jitter; points2[,1] = points2[,1]+jitter
						points3[,2] = points3[,2]-jitter; points4[,2] = points4[,2]+jitter
						points = rbind(points, points1, points2, points3, points4)
					}
				hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,])
				ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
				rast1 = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/Elevation_",study_areas[j],"_0.008.asc"))
				rast1[!is.na(rast1[])] = 0
				rast2 = rasterize(sps, rast1, getCover=T)
				rast2[rast2[]==0] = NA; rast2[is.na(rast1[])] = NA
				if (sum(is.na(raster::extract(rast2, tab[,c("startLon","startLat")]))) != 0) print(c(i,j))
				if (sum(is.na(raster::extract(rast2, tab[,c("endLon","endLat")]))) != 0) print(c(i,j))
				backgroundRasters[[j]] = rast2
			}
		backgroundRaster = merge(backgroundRasters[[1]],backgroundRasters[[2]])
		backgroundRasters_list[[i]] = backgroundRaster
	}
analyses = c("LASV_L_align_3_all","LASV_S_align_3_all")
model = "gamma"; showingPlots = FALSE; nodesOnly = FALSE; newPlot = FALSE; n1 = 100; n2 = 10
for (i in 1:length(analyses))
	{
		segment = unlist(strsplit(analyses[i],"_"))[2]
		envVariables = list(backgroundRasters_list[[i]])
		localTreesDirectory = paste0(wdb4,"/",analyses[i],"_extractions")
		mostRecentSamplingDatum = 0
		for (j in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"))
				if (mostRecentSamplingDatum < max(tab[,"endYear"])) mostRecentSamplingDatum = max(tab[,"endYear"])
			}
		trees = readAnnotatedNexus(paste0(wdb4,"/LASV_",segment,"_align_3_all_",nberOfExtractionFiles,".trees"))
			# (n.b.: before reading, remove "'" and spaces in tip labels)
		log = read.table(paste0(wdb4,"/LASV_",segment,"_align_3_all_",nberOfExtractionFiles,".log"), header=T)
		cladeTipLabels = list()
		for (j in 1:length(clades))
			{
				tab = read.csv(paste0(wdb4,"/LASV_",segment,"_alignment_3_",clades[j],".csv"), header=T)
				cladeTipLabels[[j]] = gsub(" ","_",tab[,paste0("segment_",segment)])
			}
		for (j in 1:nberOfExtractionFiles)
			{
				tree = trees[[j]]; MRCAs = rep(NA, length(cladeTipLabels))
				tree$tip.label = gsub("'","",tree$tip.label)
				for (k in 1:length(cladeTipLabels))
					{
						MRCAs[k] = findMRCA(tree, as.character(cladeTipLabels[[k]]), type="node")
					}
				if (segment == "L")
					{
						MRCAs[length(MRCAs)+1] = which(tree$tip.label=="pinneo_NGA-Borno_Lassa_LS_1969")
					}
				internalNodesToExclude = c(MRCAs); buffer = c()
				for (k in 1:length(MRCAs))
					{
						internalNodesToExclude = c(internalNodesToExclude, tree$edge[which(tree$edge[,2]==MRCAs[k]),1])
					}
				internalNodesToExclude = unique(internalNodesToExclude)
				for (k in 1:length(internalNodesToExclude))
					{
						reachingRoot = FALSE
						internalNode = internalNodesToExclude[k]
						while (reachingRoot == FALSE)
							{
								index = which(tree$edge[,2]==internalNode)
								if (length(index) == 1)
									{
										buffer = c(buffer, tree$edge[index,1])
										internalNode = tree$edge[index,1]
									}	else		{
										reachingRoot = TRUE
									}
							}
					}
				internalNodesToExclude = unique(c(internalNodesToExclude, buffer))				
				tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T)
				ancestID = which(!tab[,"node1"]%in%tab[,"node2"])[1]
				ancestPosition = c(tab[ancestID,"startLon"], tab[ancestID,"startLat"])
				rates = c(); geoDists = matrix(nrow=dim(tab)[1], ncol=1)
				for (k in 1:length(tree$annotations))
					{
						rates = c(rates, tree$annotations[[k]]$location.rate)
					}
				for (k in 1:dim(tab)[1])
					{
						x1 = cbind(tab[k,"startLon"], tab[k,"startLat"])
						x2 = cbind(tab[k,"endLon"], tab[k,"endLat"])
						geoDists[k,1] = rdist.earth(x1, x2, miles=F, R=NULL)
					}
				col11 = log[j,"treeLengthPrecision1"]
				col12 = log[j,"treeLengthPrecision3"]
				col22 = log[j,"treeLengthPrecision2"]
				my_prec = c(col11, col12, col12, col22)
				halfDF = log[j,"location.halfDF"]
				if (model == "cauchy") reciprocalRates = TRUE
				if (model == "gamma") reciprocalRates = TRUE
				if (model == "logN") reciprocalRates = FALSE
				tab = tab[order(tab[,"startYear"]),]
				cor = cor((tab[,"endLon"]-tab[,"startLon"]),(tab[,"endLat"]-tab[,"startLat"]))
				my_var = solve(matrix(my_prec,nrow=2))		
				sigma1 = sqrt(my_var[1,1]); sigma2 = sqrt(my_var[2,2])
				sigmas = c(sigma1, sigma2)
				cor = my_var[1,2]/(sqrt(my_var[1,1])*sqrt(my_var[2,2]))
				simulationsDirectory = localTreesDirectory
				fixedNodes = internalNodesToExclude
				simulation = simulatorRRW1(tree, rates, sigmas, cor, envVariables, mostRecentSamplingDatum,
									       ancestPosition, reciprocalRates, n1, n2, showingPlots, newPlot, fixedNodes)
				if (showingPlots == TRUE)
					{
						dev.new(); plot(backgroundRasters_list[[i]], col="gray90", box=F, ann=F, axes=F, legend=F)
						for (k in 1:dim(tab)[1])
							{
								segments(tab[k,"startLon"], tab[k,"startLat"], tab[k,"endLon"], tab[k,"endLat"], col="red", lwd=0.2)
							}
					}
				if (showingPlots == TRUE)
					{
						dev.new(); plot(backgroundRasters_list[[i]], col="gray90", box=F, ann=F, axes=F, legend=F)
						for (k in 1:dim(simulation)[1])
							{
								segments(simulation[k,"startLon"], simulation[k,"startLat"], simulation[k,"endLon"], simulation[k,"endLat"], col="red", lwd=0.2)
							}
					}
				file = as.character(paste(simulationsDirectory,"/TreeSimulations_",j,".csv",sep=""))
				write.csv(simulation, file, row.names=F, quote=F)
			}
	}
setwd(paste(wd,wdb4,sep="/"))
for (i in 1:length(segments))
	{
		localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_all_1_extractions")
		for (j in 1:length(clades))
			{
				localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_",clades[j],"_extractions")
				for (k in 1:nberOfExtractionFiles)
					{
						tab1 = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",k,".csv"), header=T)
						tab2 = read.csv(paste0(localTreesDirectory2,"/TreeExtractions_",k,".csv"), header=T)
						branches1 = paste(round(tab1[,"startYear"],4),round(tab1[,"endYear"],4),round(tab1[,"length"],4),sep="_")
						branches2 = paste(round(tab2[,"startYear"],4),round(tab2[,"endYear"],4),round(tab2[,"length"],4),sep="_")
						tab3 = tab1[which(branches1%in%branches2),]
						if (dim(tab3)[1]!=dim(tab2)[1]) print(c(i,j,k))
						write.csv(tab3, paste0(localTreesDirectory2,"/TreeSimulations_",k,".csv"), row.names=F, quote=F)
					}			
			}		
	}
for (i in 1:length(segments))
	{
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_NGA_extractions")
		for (j in 1:nberOfExtractionFiles)
			{
				tab = c()
				for (k in 2:length(clades))
					{
						localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_",clades[k],"_extractions")
						tab = rbind(tab, read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",j,".csv"), head=T))
					}
				write.csv(tab, paste0(localTreesDirectory2,"/TreeSimulations_",j,".csv"), row.names=F, quote=F)
			}
	}
for (i in 1:length(segments))
	{
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_all_2_extractions")
		dir.create(file.path(localTreesDirectory2), showWarnings=F)
		for (j in 1:nberOfExtractionFiles)
			{
				tab = c()
				for (k in 1:length(clades))
					{
						localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_",clades[k],"_extractions")
						tab = rbind(tab, read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",j,".csv"), head=T))
					}
				write.csv(tab, paste0(localTreesDirectory2,"/TreeSimulations_",j,".csv"), row.names=F, quote=F)
			}
	}
for (i in 1:length(segments))
	{
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_tips2_extractions")
		dir.create(file.path(localTreesDirectory2), showWarnings=F)
		for (j in 1:nberOfExtractionFiles)
			{
				localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_all_2_extractions")
				tab = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",j,".csv"), head=T)
				tab = tab[which(!tab[,"node2"]%in%tab[,"node1"]),]
				write.csv(tab, paste0(localTreesDirectory2,"/TreeSimulations_",j,".csv"), row.names=F, quote=F)
			}
	}
for (i in 1:length(segments))
	{
		localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_all_2_extractions")
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_all_2_extract_100")
		dir.create(file.path(localTreesDirectory2), showWarnings=F); n = 0
		for (j in seq(10,nberOfExtractionFiles,10))
			{
				tab = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",j,".csv"), head=T); n = n+1
				write.csv(tab, paste0(localTreesDirectory2,"/TreeSimulations_",n,".csv"), row.names=F, quote=F)
			}
		localTreesDirectory1 = paste0("LASV_",segments[i],"_align_3_tips2_extractions")
		localTreesDirectory2 = paste0("LASV_",segments[i],"_align_3_tips2_extract_100")
		dir.create(file.path(localTreesDirectory2), showWarnings=F); n = 0
		for (j in seq(10,nberOfExtractionFiles,10))
			{
				tab = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",j,".csv"), head=T); n = n+1
				write.csv(tab, paste0(localTreesDirectory2,"/TreeSimulations_",n,".csv"), row.names=F, quote=F)
			}
	}

	# 6.2. Plotting the different environmental rasters
	
envVariableNames = c("IGBP_forests","IGBP_grasslands","IGBP_savannas","IGBP_croplands",
					 "Elevation","Annual_mean_temperature","Annual_precipitation","Pop_density")
envVariableNames1 = c("","","","","","Annual","Annual","Human pop.")
envVariableNames2 = c("Forests","Grasslands","Savannas","Croplands","Elevation","mean temp.","precipitation","density (log10)")
rS = list(); cols1 = list(); cols2 = list(); merged = list(); borders = list(); extents = list(e_MRU, e_NGA); colour1 = "white"
for (i in 1:length(study_areas))
	{
		buffer = list(); raster_names = paste0(envVariableNames,"_",study_areas[i],"_0.04")
		for (j in 1:length(raster_names))
			{
				buffer[[j]] = crop(raster(paste0(wdb6,"/Environmental_files/LASV_rasters/",raster_names[j],".asc")), extents[[i]])
				if (j == 5) buffer[[5]][buffer[[5]][]<0] = 0
				if (j == 7) buffer[[7]][] = buffer[[7]][]/100
				if (j == 8) buffer[[8]][] = log10(1+buffer[[8]][])
				# (legend: temperature in °C and precipitation in meters)
			}
		rS[[i]] = buffer
		buffer = readOGR(dsn=paste0(wdb6,"/Environmental_files/Natural_Earth/"), layer="International_borders")
		borders[[i]] = crop(gSimplify(buffer,0.02), extents[[i]])
	}
for (i in 1:length(envVariableNames))
	{
		merged[[i]] = raster::merge(rS[[1]][[i]],rS[[2]][[i]])
	}
for (i in 1:length(study_areas))
	{
		buffer1 = list(); buffer2 = list()
		r1 = merged[[1]]; r2 = rS[[i]][[1]]; colour2 = "chartreuse4"
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[1]] = colorRampPalette(c(colour1,colour2),bias=1)(101)[index1:index2]
		buffer2[[1]] = colorRampPalette(c(colour1,colour2),bias=1)(101)
		r1 = merged[[2]]; r2 = rS[[i]][[2]]; colour2 = "olivedrab3"
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[2]] = colorRampPalette(c(colour1,colour2),bias=1)(101)[index1:index2]
		buffer2[[2]] = colorRampPalette(c(colour1,colour2),bias=1)(101)
		r1 = merged[[3]]; r2 = rS[[i]][[3]]; colour2 = "darkkhaki"
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[3]] = colorRampPalette(c(colour1,colour2),bias=1)(101)[index1:index2]
		buffer2[[3]] = colorRampPalette(c(colour1,colour2),bias=1)(101)
		r1 = merged[[4]]; r2 = rS[[i]][[4]]; colour2 = "navajowhite4"
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[4]] = colorRampPalette(c(colour1,colour2),bias=1)(101)[index1:index2]
		buffer2[[4]] = colorRampPalette(c(colour1,colour2),bias=1)(101)
		r1 = merged[[5]]; r2 = rS[[i]][[5]]
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[5]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(101)[index1:index2]
		buffer2[[5]] = colorRampPalette(brewer.pal(9,"YlOrBr"))(101)
		r1 = merged[[6]]; r2 = rS[[i]][[6]]
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[6]] = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101][index1:index2]
		buffer2[[6]] = colorRampPalette(brewer.pal(9,"YlOrRd"))(151)[1:101]
		r1 = merged[[7]]; r2 = rS[[i]][[7]]
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[7]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)[index1:index2]
		buffer2[[7]] = colorRampPalette(brewer.pal(9,"YlGnBu"))(101)
		r1 = merged[[8]]; r2 = rS[[i]][[8]]
		index1 = (((min(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		index2 = (((max(r2[],na.rm=T)-min(r1[],na.rm=T))/(max(r1[],na.rm=T)-min(r1[],na.rm=T)))*100)+1
		buffer1[[8]] = c("#FFFFFF",colorRampPalette(brewer.pal(9,"BuPu"))(101)[index1:index2])
		buffer2[[8]] = c("#FFFFFF",colorRampPalette(brewer.pal(9,"BuPu"))(101))
		cols1[[i]] = buffer1; cols2[[i]] = buffer2
	}
dev.new(width=9, height=4); xmin = -14; xmax = -5; ymin = 5; ymax = 12
par(mfrow=c(2,4), oma=c(2,2.5,1,0.3), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
labsX = c(expression(14*degree*W), expression(5*degree*W))
labsY = c(expression(5*degree*N), expression(12*degree*N))
for (j in 1:length(rS[[1]]))
	{
		plot(rS[[1]][[j]], bty="n", box=F, axes=F, legend=F, col=cols1[[1]][[j]], colNA="gray90")
		if (j%in%c(5:8))
			{
				axis(1, c(xmin,xmax), labels=labsX, pos=ymin(rS[[1]][[j]]), col="gray30", cex.axis=0.7, col.axis="gray30",
					 lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.02, mgp=c(0,0.15,0))
			}
		if (j%in%c(1,5))
			{
				axis(2, c(ymin,ymax), labels=labsY, pos=xmin(rS[[1]][[j]]), col="gray30", cex.axis=0.7, col.axis="gray30",
					 lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.02, mgp=c(0,0.27,0))
			}
		legendRast = raster(as.matrix(c(min(merged[[j]][],na.rm=T), max(merged[[j]][],na.rm=T))))
		plot(legendRast, legend.only=T, add=T, col=cols2[[1]][[j]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.050,0.065,0.11,0.4), adj=3,
			 axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
		plot(borders[[1]], add=T, lwd=0.2, col="gray30", lty=1)
		if (nchar(envVariableNames1[j] > 0)) mtext(envVariableNames1[j], side=1, adj=0, line=-3.1, at=-13.2, cex=0.55, font=1, col="gray30")
		if (nchar(envVariableNames2[j] > 0)) mtext(envVariableNames2[j], side=1, adj=0, line=-2.3, at=-13.2, cex=0.55, font=1, col="gray30")
		rect(xmin(rS[[1]][[j]]), ymin(rS[[1]][[j]]), xmax(rS[[1]][[j]]), ymax(rS[[1]][[j]]), lwd=0.2, border="gray30")
	}
dev.new(width=9, height=4); xmin = 3; xmax = 14; ymin = 4; ymax = 13
par(mfrow=c(2,4), oma=c(2,2.5,1,0.3), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
labsX = c(expression(3*degree*W), expression(14*degree*W))
labsY = c(expression(4*degree*N), expression(13*degree*N))
for (j in 1:length(rS[[2]]))
	{
		plot(rS[[2]][[j]], bty="n", box=F, axes=F, legend=F, col=cols1[[2]][[j]], colNA="gray90")
		if (j%in%c(5:8))
			{
				axis(1, c(xmin,xmax), labels=labsX, pos=ymin(rS[[2]][[j]]), col="gray30", cex.axis=0.7, col.axis="gray30",
					 lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.02, mgp=c(0,0.15,0))
			}
		if (j%in%c(1,5))
			{
				axis(2, c(ymin,ymax), labels=labsY, pos=xmin(rS[[2]][[j]]), col="gray30", cex.axis=0.7, col.axis="gray30",
					 lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.02, mgp=c(0,0.27,0))
			}
		plot(borders[[2]], add=T, lwd=0.2, col="gray30", lty=1)
		rect(xmin(rS[[2]][[j]]), ymin(rS[[2]][[j]]), xmax(rS[[2]][[j]]), ymax(rS[[2]][[j]]), lwd=0.2, border="gray30")
	}
		
	# 6.3. Analysing the impact of environmental factors on the dispersal velocity of lineages

setwd(wd)

		# 6.3.1. Analyses based on continuous environmental rasters

analyses = c("LASV_L_align_3_all_2","LASV_S_align_3_all_2","LASV_L_align_3_tips2","LASV_S_align_3_tips2")
nberOfExtractionFiles = 100; nberOfRandomisations = 0; randomProcedure = 3; fourCells = FALSE
showingPlots = FALSE; nberOfCores = 10; nberOfCores_CS = 1; OS = "Unix"
envVariableNames = c("IGBP_croplands","IGBP_forests","IGBP_grasslands","IGBP_savannas",
					 "Elevation","Annual_mean_temperature","Annual_precipitation","Pop_density")
for (i in 1:length(analyses))
	{
		studyArea = unlist(strsplit(gsub("_2","",analyses[i]),"_"))[length(unlist(strsplit(gsub("_2","",analyses[i]),"_")))]
		if (studyArea == "tips2") studyArea = "all"
		localTreesDirectory = paste0(wdb4,"/",analyses[i],"_extract_1°°")
		raster_names = paste0(envVariableNames,"_",studyArea,"_0.04")
		envVariables = list(); resistances = list(); avgResistances = list(); c = 0
		for (k in c(10,100,1000))
			{
				for (j in 1:length(raster_names))
					{
						c = c+1
						rast = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/",raster_names[j],".asc"))
						rast[rast[]<0] = 0
						if (grepl("Pop_density",raster_names[j]))
							{
								rast[] = rast[]+1; rast[] = log10(rast[])
							}
						M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(raster_names[j], "_k", k, sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(raster_names[j],"_k",k,sep="")
						resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
					}
				for (j in 1:length(raster_names))
					{
						c = c+1
						rast = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/",raster_names[j],".asc"))
						rast[rast[]<0] = 0
						if (grepl("Pop_density",raster_names[j]))
							{
								rast[] = rast[]+1; rast[] = log10(rast[])
							}
						M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(raster_names[j], "_k", k, sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(raster_names[j],"_k",k,sep="")
						resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
					}
			}
		pathModel = 2; simulations = FALSE; outputName = paste0(analyses[i],"_seraphim_LC_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F)
		pathModel = 2; simulations = TRUE; outputName = paste0(analyses[i],"_seraphim_LC_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=T)
		pathModel = 3; simulations = FALSE; outputName = paste0(analyses[i],"_seraphim_CS_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F)
		pathModel = 3; simulations = TRUE; outputName = paste0(analyses[i],"_seraphim_CS_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=T)
	}

analyses = c("LASV_L_align_3_all_2","LASV_S_align_3_all_2","LASV_L_align_3_tips2","LASV_S_align_3_tips2")
nberOfExtractionFiles = 100; pathModels = c("Least-cost path model","Circuitscape path model")
envVariableNames = c("IGBP_forests","IGBP_grasslands","IGBP_savannas","IGBP_croplands",
					 "Elevation","Annual_mean_temperature","Annual_precipitation","Pop_density")
for (a in 1:length(analyses))
	{
		studyArea = unlist(strsplit(gsub("_2","",analyses[a]),"_"))[length(unlist(strsplit(gsub("_2","",analyses[a]),"_")))]
		raster_names = paste0(envVariableNames,"_",studyArea,"_0.04"); extractions = list(); simulations = list()
		extractions[[1]] = read.table(paste0(wdb6,"/All_seraphim_results/",analyses[a],"_seraphim_LC_extractions_linear_regression_results.txt"), header=T)
		extractions[[2]] = read.table(paste0(wdb6,"/All_seraphim_results/",analyses[a],"_seraphim_CS_extractions_linear_regression_results.txt"), header=T)
		simulations[[1]] = read.table(paste0(wdb6,"/All_seraphim_results/",analyses[a],"_seraphim_LC_simulations_linear_regression_results.txt"), header=T)
		simulations[[2]] = read.table(paste0(wdb6,"/All_seraphim_results/",analyses[a],"_seraphim_CS_simulations_linear_regression_results.txt"), header=T)
		allResults = matrix(nrow=length(envVariableNames)*2*2*3, ncol=4); kS = c(10,100,1000); CR = c("R","C"); L = 0
		colnames(allResults) = c("Regression coefficient","Q statistic","p(Q) > 0","BF")
		for (i in 1:length(pathModels))
			{
				for (j in 1:length(envVariableNames))
					{
						for (k in 1:length(CR))
							{
								for (l in 1:length(kS))
									{
										L = L+1; c1 = 0; c2 = 0; # allResults[L,1] = pathModels[i]
										# allResults[L,2] = paste0(envVariableNames[j]," (",CR[k],")"); allResults[L,3] = kS[l]
										index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
										index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
										index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(envVariableNames[j],colnames(simulations[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
										R2 = extractions[[i]][,index1]; Qe = extractions[[i]][,index2]; Qs = simulations[[i]][,index3]
										for (m in 1:length(Qe))
											{
												if (Qs[m] < Qe[m]) c1 = c1+1
											}
										p = c1/length(Qe); BF1 = p/(1-p)
										allResults[L,1] = paste0(round(median(R2),3)," [",round(quantile(R2,0.025),3)," - ",round(quantile(R2,0.975),3),"]")
										allResults[L,2] = paste0(round(median(Qe),3)," [",round(quantile(Qe,0.025),3)," - ",round(quantile(Qe,0.975),3),"]")
										allResults[L,3] = sum(Qe>0)/nberOfExtractionFiles
										if (as.numeric(allResults[L,3]) >= 0.9)
											{
												allResults[L,4] = round(BF1,1)
											}	else	{
												allResults[L,4] = "-"
											}
									}
							}
					}
			}
		write.csv(allResults, paste0("SERAPHIM_",analyses[a],"_TEMP.csv"), row.names=F, quote=F)
	}
	
		# 6.3.2. Analyses based on rivers network rasters

minS = 3; maxS = c(5,7)
for (i in 1:length(study_areas))
	{
		template = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/Elevation_",study_areas[i],"_0.008.asc"))
		rivers = crop(readOGR(dsn=paste0(wdb6,"/Environmental_files/LASV_rivers/"), layer="FAO_GeoNetwork"), extent(template))
		for (s in minS:maxS[i])
			{
				selection = subset(rivers, rivers@data[,"Strahler"]>=s)
				rast = rasterize(selection, template)
				rast[!is.na(rast[])] = 1; rast[(is.na(rast[]))&(!is.na(template[]))] = 0
				writeRaster(rast, paste0(wdb6,"/Environmental_files/LASV_rasters/Main_rivers_S",s,"_",study_areas[i],"_0.008.asc"), overwrite=T)
			}
	}
	
analyses = c("LASV_L_align_3_all_2","LASV_S_align_3_all_2")
nberOfExtractionFiles = 100; nberOfRandomisations = 0; randomProcedure = 3; fourCells = FALSE
showingPlots = FALSE; nberOfCores = 10; nberOfCores_CS = 1; OS = "Unix"
for (i in 1:length(analyses))
	{
		studyArea = unlist(strsplit(gsub("_2","",analyses[i]),"_"))[length(unlist(strsplit(gsub("_2","",analyses[i]),"_")))]
		if (studyArea == "tips2") studyArea = "all"
		localTreesDirectory = paste0(wdb4,"/",analyses[i],"_extract_100")
		if (studyArea == "all")
			{
				raster_names = c("Main_rivers_S3_all_0.008","Main_rivers_S4_all_0.008","Main_rivers_S5_all_0.008",
								 "Main_rivers_S6_all_0.008","Main_rivers_S7_all_0.008")
			}
		if (studyArea == "MRU")
			{
				raster_names = c("Main_rivers_S3_MRU_0.008","Main_rivers_S4_MRU_0.008","Main_rivers_S5_MRU_0.008")
			}
		if (studyArea == "NGA")
			{
				raster_names = c("Main_rivers_S3_NGA_0.008","Main_rivers_S4_NGA_0.008","Main_rivers_S5_NGA_0.008",
								 "Main_rivers_S6_NGA_0.008","Main_rivers_S7_NGA_0.008")
			}
		envVariables = list(); resistances = list(); avgResistances = list(); c = 0
		for (k in c(10,100,1000,10000))
			{
				for (j in 1:length(raster_names))
					{
						c = c+1
						rast = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/",raster_names[j],".asc"))
						rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(raster_names[j], "_k", k, sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(raster_names[j],"_k",k,sep="")
						resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
					}
			}
		pathModel = 2; simulations = FALSE; outputName = paste0(analyses[i],"_seraphim_LC_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=F)
		pathModel = 2; simulations = TRUE; outputName = paste0(analyses[i],"_seraphim_LC_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations=T)
	}

analyses = c("LASV_L_align_3_all_2","LASV_S_align_3_all_2")
nberOfExtractionFiles = 100; pathModels = c("Least-cost path model")
raster_names = c("Main_rivers_S3_all_0.008","Main_rivers_S4_all_0.008","Main_rivers_S5_all_0.008","Main_rivers_S6_all_0.008","Main_rivers_S7_all_0.008")
for (a in 1:length(analyses))
	{
		extractions = list(); simulations = list()
		extractions[[1]] = read.table(paste0(wdb6,"/All_seraphim_results/",analyses[a],"_seraphim_LC_extractions_linear_regression_rivers.txt"), header=T)
		simulations[[1]] = read.table(paste0(wdb6,"/All_seraphim_results/",analyses[a],"_seraphim_LC_simulations_linear_regression_rivers.txt"), header=T)
		allResults = matrix(nrow=length(raster_names)*3, ncol=7); kS = c(10,100,1000); CR = c("R"); L = 0
		colnames(allResults) = c("Path model","Environmental factor","k","Regression coefficient","Q statistic","p(Q) > 0","BF")
		for (i in 1:length(pathModels))
			{
				for (j in 1:length(raster_names))
					{
						for (k in 1:length(CR))
							{
								for (l in 1:length(kS))
									{
										L = L+1; c1 = 0; c2 = 0; allResults[L,1] = pathModels[i]
										allResults[L,2] = paste0(raster_names[j]," (",CR[k],")"); allResults[L,3] = kS[l]
										index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(raster_names[j],colnames(extractions[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
										index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(raster_names[j],colnames(extractions[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
										index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(raster_names[j],colnames(simulations[[i]]))
													   &grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
										R2 = extractions[[i]][,index1]; Qe = extractions[[i]][,index2]; Qs = simulations[[i]][,index3]
										for (m in 1:length(Qe))
											{
												if (Qs[m] < Qe[m]) c1 = c1+1
											}
										p = c1/length(Qe); BF1 = p/(1-p)
										allResults[L,4] = paste0(round(median(R2),3)," [",round(quantile(R2,0.025),3)," - ",round(quantile(R2,0.975),3),"]")
										allResults[L,5] = paste0(round(median(Qe),3)," [",round(quantile(Qe,0.025),3)," - ",round(quantile(Qe,0.975),3),"]")
										allResults[L,6] = sum(Qe>0)/nberOfExtractionFiles
										if (as.numeric(allResults[L,6]) >= 0.9)
											{
												allResults[L,7] = round(BF1,1)
											}	else	{
												allResults[L,7] = "-"
											}
									}
							}
					}
			}
		write.csv(allResults, paste0("RIVERS_",analyses[a],"_TEMP.csv"), row.names=F, quote=F)
	}

	# 6.4. Analysing the impact of rivers on the dispersal frequency of lineages

nberOfExtractionFiles = 100; interval = 10; registerDoMC(cores=5)
analyses = c("LASV_L_align_3_all_2","LASV_S_align_3_all_2")
for (i in 1:length(analyses))
	{
		if (!file.exists(paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_BayesFactors_least-cost_weights_rivers.csv")))
			{
				BFs = matrix(nrow=5*4, ncol=1); colnames(BFs) = analyses[i]
				write.csv(round(BFs,1), paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_BayesFactors_least-cost_weights_rivers.csv"), row.names=F)				
			}	else	{
				BFs = read.csv(paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_BayesFactors_least-cost_weights_rivers.csv"), header=T)
			}
		studyArea = unlist(strsplit(gsub("_2","",analyses[i]),"_"))[length(unlist(strsplit(gsub("_2","",analyses[i]),"_")))]
		if (studyArea == "tips2") studyArea = "all"
		localTreesDirectory = paste0(wdb4,"/",analyses[i],"_extract_100")
		if (studyArea == "all")
			{
				raster_names = c("Main_rivers_S3_all_0.008","Main_rivers_S4_all_0.008","Main_rivers_S5_all_0.008",
								 "Main_rivers_S6_all_0.008","Main_rivers_S7_all_0.008")
			}
		if (studyArea == "MRU")
			{
				raster_names = c("Main_rivers_S3_MRU_0.008","Main_rivers_S4_MRU_0.008","Main_rivers_S5_MRU_0.008")
			}
		if (studyArea == "NGA")
			{
				raster_names = c("Main_rivers_S6_NGA_0.008","Main_rivers_S7_NGA_0.008")
			}
		envVariables = list(); resistances = list(); avgResistances = list(); rowNames = c(); c = 0
		for (j in 1:length(raster_names))
			{
				for (k in c(10,100,1000,10000))
					{
						c = c+1
						rast = raster(paste0(wdb6,"/Environmental_files/LASV_rasters/",raster_names[j],".asc"))
						rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(raster_names[j],"_k",k,sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(raster_names[j],"_k",k,sep="")
						resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
						rowNames = c(rowNames, paste0("Main_rivers_",unlist(strsplit(raster_names[j],"_"))[3],"_k",k))
					}
			}
		if (!file.exists(paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_extractions_least-cost_weights_sum_rivers.csv")))
			{
				LC_weights_extractions = matrix(nrow=(nberOfExtractionFiles), ncol=length(envVariables)); colnames(LC_weights_extractions) = rowNames
				LC_weights_simulations = matrix(nrow=(nberOfExtractionFiles), ncol=length(envVariables)); colnames(LC_weights_simulations) = rowNames
				write.csv(LC_weights_extractions, paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_extractions_least-cost_weights_sum_rivers.csv"), row.names=F, quote=F)
				write.csv(LC_weights_simulations, paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_simulations_least-cost_weights_sum_rivers.csv"), row.names=F, quote=F)	
			}	else	{
				LC_weights_extractions = read.csv(paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_extractions_least-cost_weights_sum_rivers.csv"), header=T)
				LC_weights_simulations = read.csv(paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_simulations_least-cost_weights_sum_rivers.csv"), header=T)	
			}
		for (j in 1:length(envVariables))
			{
				trEnvVariable = transition(envVariables[[j]], function(x) 1/mean(x), directions=8)
				# trEnvVariable = transition(aggregate(envVariables[[j]],10), function(x) 1/mean(x), directions=8)
				trEnvVariableCorr = geoCorrection(trEnvVariable, type="c", multpl=F, scl=T); n = 0; buffer = list()
				buffer = foreach(k = 1:nberOfExtractionFiles) %dopar% {
						obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",k,".csv"), header=T)
						sim = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",k,".csv"), header=T)
						LC1 = diag(costDistance(trEnvVariableCorr, as.matrix(obs[,c("startLon","startLat")]), as.matrix(obs[,c("endLon","endLat")])))
						LC2 = diag(costDistance(trEnvVariableCorr, as.matrix(sim[,c("startLon","startLat")]), as.matrix(sim[,c("endLon","endLat")])))
						tmp1 = LC1; tmp2 = LC2; LC1[(!is.finite(tmp1))|(!is.finite(tmp2))] = NA; LC2[(!is.finite(tmp1))|(!is.finite(tmp2))] = NA
						cbind(sum(LC1, na.rm=T), sum(LC2, na.rm=T))
					}
				for (k in 1:length(buffer))
					{	
						LC_weights_extractions[k,j] = buffer[[k]][1,1]
						LC_weights_simulations[k,j] = buffer[[k]][1,2]
						if (buffer[[k]][1,1] < buffer[[k]][1,2]) n = n+1
					}
				p = n/nberOfExtractionFiles; BFs[j,1] = p/(1-p)
				write.csv(LC_weights_extractions, paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_extractions_least-cost_weights_sum_rivers.csv"), row.names=F, quote=F)
				write.csv(LC_weights_simulations, paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_simulations_least-cost_weights_sum_rivers.csv"), row.names=F, quote=F)
				write.csv(BFs, paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_BayesFactors_least-cost_weights_rivers.csv"), row.names=F, quote=F)
			}
	}

table = matrix(nrow=5*4, ncol=length(analyses)*2); colnames(BFs) = analyses
for (i in 1:length(analyses))
	{
		tab1 = read.csv(paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_extractions_least-cost_weights_sum_rivers.csv"), header=T)
		tab2 = read.csv(paste0(wdb6,"/All_seraphim_results/",analyses[i],"_seraphim_simulations_least-cost_weights_sum_rivers.csv"), header=T)
		for (j in 1:dim(tab1)[2])
			{
				table[j,((i-1)*2)+1] = paste0(round(median(tab1[,j]))," [",round(quantile(tab1[,j],0.025)),"-",round(quantile(tab1[,j],0.975)),"]"); n = 0
				for (k in 1:dim(tab1)[1])
					{
						if (tab1[k,j] < tab2[k,j])
							{
								if (k == 1) n = 1
								if (k >= 2) n = n+1
							}
					}
				p = n/dim(tab1)[1]; BF = p/(1-p); table[j,((i-1)*2)+2] = round(BF,1)
			}
	}
write.table(table, "Table_1_without_layout.txt", row.names=F, quote=F, sep="\t")

