############################################################################################################
### Script to run a root-to-tip against sampling date regression on a dated phylogeny            														        ###
### (c) 2015 Gemma GR Murray, John J Welch                                                              					 													 ###
############################################################################################################
############################################################################################################

############################################################################################################
## Randomise dates for permutation tests
create.randomised.dates <- function(phy, dates, auto.round.dates=F, rounded.dates=NULL, nreps=1000, use.clusters=T, clusters=NULL, do.plot=F,stat='nve-rms')
{
	udates 		<- unique(dates)
	
	# Generate unique random orderings, either across sequences or clusters 
	if(!use.clusters) { # Just get unique random permutations of the dates across sequences
		# Number of distinct ways of arranging 'udates' dates among n sequences. May limit the number of replicates.
		nposs <- exp(lfactorial(length(dates))-sum(lfactorial(table(dates))))
		if(nposs<nreps)
			nreps <- nposs
		# Generate 'nreps' sets of unique random permutations of the dates among the sequences
		S <- perm.rand(nreps,mm=rbind(match(dates,udates),dates),luc=NULL)
	} else {	
		taxnames <- phy$tip.label
		
		#	If clusters are specified, use these, and create plot describing clusters if 'do.plot'
		if(!is.null(clusters)) {
			clusters <- clusters
			if(do.plot) {plot.clusters(phy,clusters)}
			}	

		#	If clusters are not specified, define clusters as single-date clades using 'dates' (rounded if 'auto.round,dates'=T) or 'rounded.dates' if specified, and create plot describing clusters if 'do.plot'		
		if(is.null(clusters)) {clusters <- get.clusters(phy,dates, auto.round.dates=auto.round.dates, rounded.dates=rounded.dates,  do.plot=do.plot, stat=stat)
		}
						
		luc <- length(unique(clusters))
		if(luc==1) stop("Only one cluster in data. Dates need to be finer grained.")
		cluster_dates_single<-rbind(clusters,dates)[,!duplicated(clusters)]
		cluster_dates<-rbind(clusters,dates)

		# Number of distinct ways of arranging dates among 'luc' clusters. May limit the number of replicates.
		nposs <-exp(lfactorial(luc)-sum(lfactorial(table(cluster_dates_single[2,]))))
		if(nposs<nreps) {nreps <- nposs-1}
			
		# Generate 'nreps' sets of unique random permutations of the dates among the clusters
		S <- perm.rand(nreps,cluster_dates,luc,clusters)

	}
	u <- matrix(S,dim(S))
	return(u)
}
##############################################################################################################	
## Run the regression analysis, using the specified statistic (i.e., maximise r or minimises residual mean squares)
pathogen.regression <- function(phy,dates,stat='nve-rms',do.path.plot=F,reroot=T)
{	
	library(ape)
	
	# Root the tree to optimise model fit based on chosen stat
	if(!is.rooted(phy)| reroot==T) {
		phy.allroots <- get.all.roots(phy)	
		phy <- get.best.root(phy.allroots,dates,stat)
	}
	
	# Root-to-tip distances	
	div <- dist.nodes(phy)[(length(phy$tip.label)+1),c(1:length(phy$tip.label))]
	
	# Fit a linear model
	mylm <- lm(div~dates)
	slm  <- summary(mylm)
	r	 <-sqrt(slm$r.squared)*sign(slm$coefficients[2,1])	
	slope	   <-slm$coefficients[2]
	intercepty <-slm$coefficients[1]
	tmrca<-(-intercepty/slope)
	
	# Plot regression
	if(do.path.plot)
	{
		if(tmrca<=max(dates)){
		plot(div~dates,xlab="Dates",ylab="Root-to-tip distance",main="Regression",pch=19,col="dodgerblue",xlim=c(tmrca,max(dates)),ylim=c(0,max(div)))}
		if(tmrca>max(dates)){
		plot(div~dates,xlab="Dates",ylab="Root-to-tip distance", main="Regression",pch=19,col="dodgerblue")}
		abline(a=intercepty,b=slope,col='darkgrey')
		mtext(bquote(atop(paste(italic("r")," = ", .(round(r,digits=2)),sep=""), paste("Time of MRCA = ",.(round(tmrca)), sep=""))), side=1, line=7, outer=F, cex=0.8)		
	}
	
	regression<-c()
	regression$phy<-phy
	regression$dates<-dates
	regression$date_mrca<-tmrca
	regression$r<-r
	regression$rms<-sum((mean(residuals(mylm))-residuals(mylm))^2)/mylm$df.residual		
	return(regression)
}
############################################################################################################	
## Run the randomisation test on the regression analysis
pathogen.permutation.test <- function(phy, dates, auto.round.dates=F, rounded.dates=NULL, reroot=T, stat='nve-rms', nreps=1000, use.clusters=T, clusters=NULL, print.progress=T, output.rooted.tree=T)
{
	library(ape)
	source('Temporal_signal_functions/tempSignalFunctions.r')
	
	if(any(!is.numeric(dates))) stop("Dates need to be numeric.")
	if(!is.null(rounded.dates)) if(any(!is.numeric(rounded.dates))) stop("Rounded dates need to be numeric.")
	
	if(phy$Nnode < (length(phy$tip.label)-1)) {phy<-multi2di(phy); cat("\n\nPolytomies in tree were resolved randomly.\n\n")}
	
	if(use.clusters)
		{par(mfrow=c(1,3),mar=c(10,4.5,3,1.5))}
	else
		{par(mfrow=c(1,2),mar=c(10,4.5,3,1.5))}
	
	# Root the tree to optimise model fit based on chosen stat
	if(!is.rooted(phy) | reroot==T) {
		phy.allroots <- get.all.roots(phy)	
		phy   <- get.best.root(phy.allroots, dates, stat)
	}
	
	# Linear regression of root-to-tip distance against sampling date
	regression.test <- pathogen.regression(phy=phy, dates=dates, stat=stat, do.path.plot=T, reroot=F)
	true_r <- regression.test$r
	
	# Create sets of randomised dates
	t_random<-create.randomised.dates(phy=phy, dates=dates, rounded.dates=rounded.dates, nreps=nreps, use.clusters=use.clusters, clusters=clusters, do.plot=T, auto.round.dates=auto.round.dates, stat=stat)
	
	# Repeat linear regression with randomised dates
	nperms <- dim(t_random)[2]
	randomised_r<-rep(NA,nperms)	
	for(i in 1:nperms)
	{
		randomised_r[i]<-pathogen.regression(phy,as.numeric(t_random[,i]),stat=stat,do.path.plot=F, reroot=T)$r	
		if(print.progress){print(paste("Randomisation ",i, " done.",sep=""))}	
	}	
	p_value<-length(which(true_r <=c(true_r,randomised_r)))/(nperms+1)
	min_p_value<-1/(nperms+1)
	
	# Create a plot of the permutation test
	hist(randomised_r,xlim=c(min(c(randomised_r,true_r)), max(c(randomised_r,true_r))), xlab=bquote(italic("r")), ylab="Frequency", main="Randomised data sets", col="gray80", border=F)
	box()
	abline(v=true_r,col="red")	
	mtext(bquote(atop(paste(italic("p"),"-value = ",.(round(p_value,digits=3)),sep=""),paste("min ", italic("p"),"-value = ",.(round(min_p_value,digits=3)),sep=""))),side=1,line=7,outer=F,cex=0.8)		
	
	if(output.rooted.tree==T) {
		write.tree(phy,"rooted tree.tre")
	}
	
	if(!is.null(clusters)) {
			clusters <- clusters
			}				
	if(is.null(clusters)) {
		clusters <- get.clusters(phy,dates, auto.round.dates=auto.round.dates, rounded.dates=rounded.dates,  do.plot=F, stat=stat)
		}
						
	if(use.clusters==T) {
		cluster_tip<-cbind(cbind(clusters, phy$tip.label))
		colnames(cluster_tip)[2]<-"sequences"	
		regression.test$permutation_clusters<-cluster_tip}
	if(!is.null(rounded.dates)) {regression.test$rounded_dates<-rounded.dates}

	regression.test$p_value<-p_value
	regression.test$min_p_value<-min_p_value
	regression.test$clustered<-use.clusters
	regression.test$rounded.dates<-rounded.dates
	regression.test$randomised_r<-randomised_r
	
	return(regression.test)
}
############################################################################################################
### Instructions for implementation of function																																	  ###
############################################################################################################

# setwd('') # To implement these functions, please place all scripts and the tree file in the same working directory and run the following commands.

# Read in a phylogeny estimated without use of dates (e.g. a NJ tree)
# library(ape)
# phy <- '12_ST22.tre' # An example of a confounded data set. 
# phy <- read.tree(phy)

#A list of dates are required in numeric format in the same order as the tip labels of the tree (phy$tip.label)	
# dates <- as.numeric(gsub(".*_","",phy$tip.label))		# This provides a list of dates in a numerical format such that the order matches phy$tip.label, this function works only if dates are provided in tip labels after the last '_'. 

#Regression test function
#pathogen.permutation.test(phy, dates, auto.round.dates=T, rounded.dates= NULL,  reroot=T, stat='nve-rms', nreps=100, use.clusters=T, clusters=NULL, print.progress=T, output.rooted.tree=T)

# Read the read_me.pdf document for a description of all the arguments of this function.
############################################################################################################