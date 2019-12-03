############################################################################################################
### Script to run a Mantel test of confounding of genetic distance and sampling dates for both unclustered and clustered data          ###
### (c) 2015 Gemma GR Murray, John J Welch                                                                                                                                     ###
############################################################################################################
############################################################################################################

############################################################################################################
## Read in a nucleotide alignment in fasta format, with dates either specified as a vector in the order of the sequences, or as the last element of the sequence names, after a specified character.
read.data <- function(dnafile, dateschar ="_",dates=NULL)
{
	library(ape)
	library(seqinr)
	d<-read.dna(dnafile,format="fasta",as.matrix=T)
	g<-dist.dna(d,as.matrix=T)
	
	names<-rownames(g)
	if(is.null(dates))
		{
		x<-paste(".*\\", dateschar,sep="")
		dates<-as.numeric(gsub(x,"",names))
		}

	return(list(g = g, dates=dates))	
}
############################################################################################################
## Read in a tree file in nexus or newick format, with dates either specified as a vector in the order of the tip labels, or as the last element of the tip names, after a specified character.
read.data.tree <- function(treefile, dateschar ="_", dates=NULL, type="nexus")
{
	library(ape)
	library(seqinr)
	if(type=="nexus") d<-read.nexus(treefile)
	if(type=="newick") d<- read.tree(treefile)
	g<-cophenetic(d)
		
	names<-rownames(g)
	if(is.null(dates))
		{
			x<-paste(".*\\", dateschar,sep="")
			dates<-as.numeric(gsub(x,"",names))
		}

	return(list(g = g, dates=dates))	
}
############################################################################################################
## Define clusters based on monophyly and shared sampling date.
get.clusters.Mt <- function(inputfile, dateschar="_", dates=NULL, auto.round.dates=F, rounded.dates=NULL, datatype="dna.fasta", do.plot=T,approx.dates=F, reroot=T, stat='nve-rms')
{
	library(ape)
	library(seqinr)
	source('Temporal_signal_functions/tempSignalFunctions.r')
	if(datatype=="dna.fasta") {x <- read.data(inputfile, dateschar); phy <- nj(x$g)}
	if(datatype=="tree.nexus") {x <- read.data.tree(inputfile, dateschar, dates, 'nexus'); phy <- read.nexus(inputfile)}
	if(datatype=="tree.newick") {x <- read.data.tree(inputfile, dateschar, dates, 'newick'); phy <- read.tree(inputfile)}
	if(phy$Nnode < (length(phy$tip.label)-1)) {phy<-multi2di(phy); cat("\n\nPolytomies in tree were resolved randomly.\n\n")}
	
	clusters <- get.clusters(phy,x$dates,auto.round.dates=auto.round.dates, rounded.dates=rounded.dates, do.plot, stat=stat, reroot=reroot)
	return(clusters)	
}
############################################################################################################
##Perform a Mantel test of confounding of temporal and genetic structure based on a dated nucleotide alignment or a dated-tip phylogeny. The test can be applied to individual samples or clusters of samples.
mantel.confounding.test <- function(inputfile, datatype="dna.fasta", dateschar="_", dates=NULL, auto.round.dates=F, rounded.dates=NULL, nreps=1000, pval=0.05, use.clusters=F, clusters=NULL, do.plot=F, reroot=T, stat='nve-rms')
{
	library(ade4)
	source('Temporal_signal_functions/tempSignalFunctions.r')
	if(datatype=="dna.fasta") {x <- read.data(inputfile, dateschar, dates); phy <- nj(x$g)}
	if(datatype=="tree.nexus") {x <- read.data.tree(inputfile, dateschar, dates, 'nexus'); phy <- read.nexus(inputfile)}
	if(datatype=="tree.newick") {x <- read.data.tree(inputfile, dateschar, dates, 'newick'); phy <- read.tree(inputfile)}
	if(any(!is.numeric(x$g))) stop("Dates need to be numeric.")

	# Clustering for test
	if(use.clusters==T)
	{
		if(do.plot) {par(mfrow=c(1,2))}
		# Define clusters
		if(!is.null(clusters)) 
		{
			clusters <- clusters
			if(do.plot) {plot.clusters(phy,clusters)}
		}
		if(is.null(clusters)) 
		{
			clusters <- get.clusters.Mt (inputfile,dateschar= dateschar, dates=dates, auto.round.dates=auto.round.dates, rounded.dates=rounded.dates, datatype=datatype, do.plot=do.plot, reroot=reroot)
		}
		if(length(unique(clusters))==1) stop("Only one cluster in data. Dates need to be finer grained.")
		# Collapse distance matrices into clusters
		newdates<-c()
		cl <- clusters
		for(i in 1:	length(unique(clusters)))
		{
			if(length(which(cl ==i))>1){
		 	newgrow <- colSums(x$g[which(cl ==i),-which(cl ==i)])/length(which(cl ==i))
			x$g <- x$g[-which(cl ==i),-which(cl ==i)]
			x$g <- cbind(x$g,newgrow)
			x$g <- rbind(x$g,c(newgrow,0))
			rownames(x$g) <- c(rownames(x$g)[-length(rownames(x$g))],i)
			colnames(x$g) <- c(colnames(x$g)[-length(colnames(x$g))],i)	
			cl <- c(cl[-which(cl ==i)],NA)
		}	
		if(length(which(cl ==i))==1)
		{
			rownames(x$g)[which(cl ==i)]<-i
			colnames(x$g)[which(cl ==i)]<-i
		}
			newdate <- sum(x$dates[which(clusters==i)])/length(which(clusters==i))
			newdates <- c(newdates,newdate)	
		}	
		x$dates<-newdates[as.numeric(rownames(x$g))]
	}
	
	# Run Mantel test
	m <- suppressWarnings(mantel.rtest(as.dist(x$g),dist(x$dates),nrepet=nreps))
	p <- m$pval
	cat('r =',round(m$obs,3),'p =',round(p,3),'\n')
	if(p > pval) cat('\tNo evidence of counfounding in this dataset.\n\n')
	if(p < pval) cat('\tEvidence of significant counfounding in this dataset. \n\tTry clustering.\n\n')
	if(do.plot) {plot(as.dist(x$g),dist(x$dates),pch=16, col=rgb(0,0,1,1/4),ylab="Temporal distance",xlab="Genetic distance", main="Test of confounding")
	mtext(bquote(atop(paste(italic("p"),"-value = ",.(round(p,digits=3))))),side=3,line=-4,outer=F)}	
}

############################################################################################################
### Instructions for implementation of function															  															              ###
############################################################################################################

#setwd('') # To implement these functions, please place all scripts and the data file in the same working directory and run the following commands.

# Specify the data file, which can be an alignment in fasta format, or a tree in nexus or newick format.
# inputfile <-'12_ST22.fasta' #an example of a confounded data set

# Run the Mantel test
#mantel.confounding.test(inputfile, datatype='dna.fasta', dateschar="_", dates=NULL, auto.round.dates=T, rounded.dates=NULL, nreps=1000, pval=0.05, do.plot=T, use.clusters=T, clusters=NULL, reroot=T, stat='nve-rms')

# Read the read_me.pdf document for a description of all the arguments of this function.
############################################################################################################