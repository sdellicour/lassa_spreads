############################################################################################################
### Script to create replicate BEAST v1.8 xml files with randomly permuted sampling dates               												 ###
### (c) 2015 Gemma GR Murray, John J Welch                                                              																		 ###
############################################################################################################


############################################################################################################
### Get the sequences in an xml
get.seqs <- function(xml)
{
	library(ape)
	ss <- grep('<sequence>',xml)
	se <- grep('</sequence>',xml)
	seqs<-c()
	for(i in 1:length(ss)){
		s <- ss[i]+2
		f <- ''
		while(s < se[i]) {
			f <- paste(f,tolower(gsub('\t','',xml[s],fixed=T)),sep='')
			s <- s+1
		}
		fa<-strsplit(f,'')[[1]]
		seqs<-rbind(seqs,fa)
	}
	rownames(seqs)<-gsub("\"/>","",gsub("\t\t\t<taxon idref=\"","",xml[ss+1],fixed=T),fixed=T)
	seqs
}
############################################################################################################	
## Randomise the dates in an xml, using the "clustered" method if the flag = T
create.randomised.xmls <- function(fname=NULL,nreps=10,use.clusters=T,do.plot=T,auto.round.dates=F, rounded.dates=NULL, clusters=NULL, stat='nve-rms')
{
	source('Temporal_signal_functions/tempSignalFunctions.r')
	if(is.null(fname))
		stop('xml file not specified')
	xml		  	<- readLines(fname)
	file.stem 	<- strsplit(fname,'.xml')[[1]]
	taxnames 	<- gsub('\t','',gsub('\">','',gsub('<taxon id=\"','',xml[grep('<taxon id=',xml,fixed=T)],fixed=T),fixed=T),fixed=T)
	datelines 	<- grep("<date value=",xml)
	udatestrings <- unique(xml[datelines])
	dates 		<- as.numeric(gsub('\" direction.*','',gsub('\t','',gsub('<date value=\"','',xml[datelines],fixed=T),fixed=T)))
	udates 		<- unique(dates)
	n 			<- length(dates)
	
	# (1) Generate unique random orderings, either across sequences or clusters 
	if(!use.clusters) { # Just get unique random permutations of the dates
		add.stem <- '.perm'
		# Number of distinct ways of arranging dates among n sequences
		nposs <- exp(lfactorial(n)-sum(lfactorial(table(dates))))
		if(nposs<nreps)
			nreps <- nposs
		# Second argument here matches the date of each sequence to udates
		S <- perm.rand(nreps,mm=rbind(match(dates,udates),dates),luc=NULL)
	} 
	if(use.clusters)  {	# Check for clusters of sequences sharing a date, and permute across clusters
		add.stem <- '.clusteredperm.'
		seqs <- get.seqs(xml)
		phy <- nj(as.dist(dist.dna(as.DNAbin(seqs))))
		if(!is.null(clusters)) {clusters <- clusters
			clusters <- clusters[match(taxnames,phy$tip.label)]
			if(do.plot) {plot.clusters(phy,clusters)}}	

		#	If clusters are not specified, define clusters as single-date clades using 'dates' (rounded if 'auto.round,dates'=T) or 'rounded.dates' if specified, and create plot describing clusters if 'do.plot'		
		if(is.null(clusters)) 
		{
			clusters <- get.clusters(phy,dates[match(phy$tip.label,taxnames)], auto.round.dates=auto.round.dates, rounded.dates=rounded.dates,  do.plot=do.plot, stat=stat)
			clusters <- clusters[match(taxnames,rownames(seqs))]
		}
			luc <- length(unique(clusters))		
		
		cluster_dates_single<-rbind(clusters,dates)[,!duplicated(clusters)]
		cluster_dates<-rbind(clusters,dates)
		
		# Number of distinct ways of arranging dates among luc clusters

		nposs <-exp(lfactorial(luc)-sum(lfactorial(table(cluster_dates_single[2,]))))
		if(nposs<nreps) {nreps <- nposs-1}

		S <- perm.rand(nreps,cluster_dates,luc,clusters)
		}
	# (2) Write and save the xmls
	for(i in 1:dim(S)[2]) {
		new.xml <- xml
		new.xml <- gsub(file.stem,paste(file.stem,add.stem,i,sep=''),new.xml)
		new.xml[datelines] <- paste("\t\t\t<date value=\"",S[,i],"\" direction=\"forwards\" units=\"years\"/>",sep="")
		writeLines(new.xml,paste(file.stem,add.stem,i,'.xml',sep=''))
	}
}
############################################################################################################
### Instructions for implementation of function																																	 ###
############################################################################################################

#setwd('') # To implement these scripts, please place all scripts and the xml file in the same working directory and run the following commands.

#xml.file.name <- '12_ST22.xml' 	#A BEAST v1.8 xml file with dates specified as a number of time units since some time in the past, i.e. direction = 'forwards'.

#Function to create randomised xmls (note the that number of permuted data sets produced is limited by the number of possible unique permutations)
#create.randomised.xmls(xml.file.name, nreps=10, use.clusters=T, do.plot=T, auto.round.dates=T, rounded.dates=NULL, clusters=NULL, stat='nve-rms')

# Read the read_me.pdf document for a description of all the arguments of this function.

############################################################################################################