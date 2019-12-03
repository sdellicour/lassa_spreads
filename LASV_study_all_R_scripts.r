library(diagram)
library(dismo)
library(gbm)
library(geosphere)
library(lubridate)
library(ncdf4)
library(pgirmess)
library(seraphim)
library(seqinr)
library(vioplot)

# 1. Preparation of LASV sequences
# 2. Preliminary phylogenetic analyses
# 3. Temporal signal analyses
# 4. Continuous phylogeographic analyses
# 5. Generating dispersal history graphs
# 6. Estimating dispersal statistics
# 7. Post-hoc landscape phylogeography
# 8. Prepration of skygrid-GLM covariates
# 9. Analysis of skygrid-GLM analyses
# 10. Species distribution modelling

wd = getwd()
wd1 = "01_sequences_preparation"wd2 = "02_preliminary_BEAST_runs"wd3 = "03_temporal_signal_analyses"wd4 = "04_RRW_phylogeography"wd5 = "05_dispersal_history_graphs"wd6 = "06_dispersal_stat_estimations"wd7 = "07_all_seraphim_analyses"wd8 = "08_covariates_preparation"wd9 = "09_skygrid_GLM_analyses"wd10 = "10_SDM_analyses_host_virus"

# 1. Preparation of LASV sequences

setwd(paste(wd,wd1,sep="/"))

	# 1.1. Renaming sequence IDs to match with metadata files

segments = c("L","S"); test = decimal_date(ymd("1986-01-26"))
metadata = read.csv("Original_sequence_files/LASV_all_the_metadata.csv", head=T, colClasses="character")
for (i in 1:length(segments))
	{
		fasta = scan(paste0("Original_sequence_files/LASV_",segments[i],"_alignment_SD.fasta"), what="", sep="\n", quiet=T)
		colNames = c(paste0("segment_",segments[i]),paste0("accession_",segments[i]),"country","admin1","location","latitude","longitude")
		tab1 = metadata[,colNames]; tab1 = tab1[which(tab1[,paste0("segment_",segments[i])]!=""),]
		sequenceIDs = gsub(">","",fasta[which(grepl(">",fasta))])
		if (length(sequenceIDs[which(sequenceIDs%in%tab1[,paste0("accession_",segments[i])])]) > 0)
			{
				cat(sequenceIDs[which(!sequenceIDs%in%tab1[,paste0("accession_",segments[i])])],sep="\n")
			}
		tab1 = tab1[which(tab1[,paste0("accession_",segments[i])]%in%sequenceIDs),]
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

	# 1.2. Identifying and discarding suspiciouly duplicated sequences

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

	# 1.3. Preparing alignments for the continuous phylogeographic analyses

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
	}

	# 1.4. Preparing alignments that only contain Nigerian sequences

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

# 2. Preliminary phylogenetic analyses

	# 2.1. Preparing and running analyses in BEAST
	
		# Models used: GTR+G substitution model, relaxed lognormal molecular clock model, constant population size coalescent model

	# 2.2. Generating alignements and metadata for each clade
	
segments = c("L","S"); clades = c("NI1","NI2","NI3","MRU")
for (i in 1:length(segments))
	{
		setwd(paste(wd,wd1,sep="/"))
		fastaAll = scan(paste0("LASV_",segments[i],"_alignment_2.fasta"), what="", sep="\n", quiet=T)
		metadata = read.csv(paste0("LASV_",segments[i],"_alignment_2.csv"), colClasses="character", header=T)
		setwd(paste(wd,wd2,sep="/"))
		for (j in 1:length(clades))
			{
				tree = read.nexus(paste0("LASV_",segments[i],"_alignment_2_",clades[j],".tree"))
				tips = gsub("'","",tree$tip.label); fasta2 = c(); tab2 = c()
				for (k in 1:length(tips))
					{
						index = which(fastaAll==paste0(">",tips[k]))
						if (length(index) != 1) print(c(i,j,k,"fasta"))
						fasta2 = c(fasta2, fastaAll[index], fastaAll[index+1])
						index = which(metadata[,1]==tips[k])
						if (length(index) != 1) print(c(i,j,k,"metadata"))
						tab2 = rbind(tab2, metadata[index,])
					}
				write(fasta2, paste0("LASV_",segments[i],"_alignment_2_",clades[j],".fasta"))
				write.csv(tab2, paste0("LASV_",segments[i],"_alignment_2_",clades[j],".csv"), row.names=F, quote=F)
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
				write(fasta3, paste0("LASV_",segments[i],"_alignment_3_",clades[j],".fasta"))
				write.csv(tab3, paste0("LASV_",segments[i],"_alignment_3_",clades[j],".csv"), row.names=F, quote=F)
				tab4 = tab3[,c(paste0("segment_",segments[i]),"latitude","longitude")]; colnames(tab4) = c("trait","latitude","longitude")
				write.table(tab4, paste0("LASV_",segments[i],"_alignment_3_",clades[j],".txt"), sep="	", row.names=F, quote=F)
			}
	}

# 3. Temporal signal analyses

setwd(paste(wd,wd3,sep="/"))
source("Temporal_signal_functions/mantelCounding.r")
source("Temporal_signal_functions/randRegression.r")
source("Temporal_signal_functions/tempSignalFunctions.r")
names = c("LASV_L_alignment","LASV_S_alignment")
for (i in 1:length(names))
	{
		if (!file.exists(paste0(names[i],".txt")))
			{
				tree = read.tree(paste0(names[i],".tree"))
				labs = tree$tip.label; years = c()
				for (j in 1:length(labs))
					{
						date = unlist(strsplit(labs[j],"_"))[length(unlist(strsplit(labs[j],"_")))]
						if (length(unlist(strsplit(date,"-"))) == 3)
							{
								years = c(years, decimal_date(ymd(as.character(date))))
							}
						if (length(unlist(strsplit(date,"-"))) == 2)
							{
								date = paste0(date,"-15")
								years = c(years, decimal_date(ymd(as.character(date))))
							}
						if (length(unlist(strsplit(date,"-"))) == 1)
							{
								years = c(years, as.numeric(date)+0.5)
							}
					}
				write.table(years, paste0(names[i],".txt"), quote=F, row.name=F, col.names=F)
			}
	}
for (i in 1:length(names)) # Mantel tests:
	{
		fasta = paste0(names[i],".fasta")
		mantel.confounding.test(inputfile=fasta, use.clusters=F, auto.round.dates=T, do.plot=T,
                        			nreps=10000, pval=".05", datatype="dna.fasta", dateschar="_")
        # (test with clustering by date within clades:)
		mantel.confounding.test(inputfile=fasta, use.clusters=T, auto.round.dates=T, do.plot=T,
								nreps=10000, pval=".05", datatype="dna.fasta", dateschar="_")
	}
pValues = c()
for (i in 1:length(names)) # regression tests:
	{
		tree = read.tree(paste0(names[i],".tree"))
		tipDates = read.table(paste0(names[i],".txt"), header=F)
		tipDates = c(as.numeric(tipDates$V1))
		test = pathogen.permutation.test(phy=tree, dates=tipDates, use.clusters=F, auto.round.dates=F, nreps=1000)
		pValues = c(pValues, test$p_value)
	}
pValues = c()
for (i in 1:length(names)) # tests with clustering by date within clades (not done here):
	{
		tree = read.tree(paste0(names[i],".tree"))
		tipDates = read.table(paste0(names[i],".txt"), header=F)
		tipDates = c(as.numeric(tipDates$V1))
		test = pathogen.permutation.test(phy=tree, dates=tipDates, use.clusters=T, auto.round.dates=T, nreps=1000)
		pValues = c(pValues, test$p_value)
	}

# 4. Continuous phylogeographic analyses

setwd(paste(wd,wd4,sep="/"))

	# 4.1. Preparing and running analyses in BEAST
	
		# Models used: GTR+G substitution model, relaxed lognormal molecular clock model, skygrid coalescent model (with 100 grid points), RRW model (gamma, jitter = 0.01)
		# Skygrid cutoff: 1600 (L-MRU), 1200 (L-NI1), 1300 (L-NI2), 1450 (L-NI3), 1700 (S-MRU), 1600 (S-NI1), 1800 (S-NI2), 1700 (S-NI3)

	# 4.2. Extraction of 1,000 tree files

analyses = c("L_3_MRU","L_3_NIG1","L_3_NIG2","L_3_NIG3","S_3_MRU","S_3_NIG1","S_3_NIG2","S_3_NIG3")
mostRecentSamplingDates = c(2014.3534246575343, 2016.035519125683, 2019.1780821917807, 2019.1041095890412,
									    2018.0, 2019.043835616438, 2019.1780821917807, 2019.1041095890412)

		# 4.2.1. Selecting 1,000 trees from the post-burn-in distribution

burnIns = c(20000000,15000000,100000000,20000000,20000000,20000000,20000000,20000000); options(scipen=9)
for (i in 1:length(analyses))
	{
		allTrees = scan(file=paste0("LASV_",analyses[i],"_RRW.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		lineIndex1 = which(grepl(paste0("tree STATE_0 "),allTrees))
		lineIndex2 = which(grepl(paste0("tree STATE_",burnIns[i]," "),allTrees))
		lineIndex3 = which(grepl(paste0("End;"),allTrees))[2]
		samplingFrequency = floor(((lineIndex3-1)-(lineIndex2))/nberOfExtractionFiles)
		treesToSample = seq(((lineIndex3-1)-((nberOfExtractionFiles-1)*samplingFrequency)),(lineIndex3-1),samplingFrequency)
		selectedTrees = c(allTrees[1:(lineIndex1-1)],allTrees[treesToSample],"End;")
		write(selectedTrees, paste0("LASV_",analyses[i],"_",nberOfExtractionFiles,".trees"))
	}	

		# 4.2.2. Extractions of the spatio-temporal information contained in posterior trees

for (i in 1:length(analyses))
	{	
		allTrees = scan(file=paste0("LASV_",analyses[i],"_1000.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		segment = unlist(strsplit(analyses[i],"_"))[1]; area = unlist(strsplit(analyses[i],"_"))[3]
		localTreesDirectory = paste0("Tree_extractions/LASV2_",area,"_seg",segment)
		burnIn = 0; randomSampling = FALSE
		nberOfTreesToSample = nberOfExtractionFiles
		mostRecentSamplingDatum = mostRecentSamplingDates[i]
		coordinateAttributeName = "location"; nberOfCores = 1
		treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
	}

		# 4.2.3. Generating the maximum clade credibility (MCC) trees with TreeAnnotator

for (i in 1:length(analyses))
	{
		input = paste0("LASV_",analyses[i],"_1000.trees"); output = paste0("LASV_",analyses[i],"_MCC.tree")
		system(paste0("BEAST_ver_1.10.4/bin/treeannotator -burninTrees 0 -heights keep ",input," ",output))
	}

		# 4.2.4. Extractions of the spatio-temporal information contained in MCC trees

source("Divers_R_functions/mccTreeExtraction.r")
for (i in 1:length(analyses))
	{
		mcc_tre = readAnnotatedNexus(paste0("LASV_",analyses[i],"_MCC.tree"))
		mostRecentSamplingDatum = mostRecentSamplingDates[i]
		mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc_tab, paste0("LASV_",analyses[i],"_MCC.csv"), row.names=F, quote=F)
	}

		# 4.2.5. Gathering all the extraction files per study area

for (i in c(2,6))
	{
		segment1 = unlist(strsplit(analyses[i+0],"_"))[1]; area1 = unlist(strsplit(analyses[i+0],"_"))[3]
		localTreesDirectory1 = paste0("Tree_extractions/LASV2_",area1,"_seg",segment1)
		segment2 = unlist(strsplit(analyses[i+1],"_"))[1]; area2 = unlist(strsplit(analyses[i+1],"_"))[3]
		localTreesDirectory2 = paste0("Tree_extractions/LASV2_",area2,"_seg",segment2)
		segment3 = unlist(strsplit(analyses[i+2],"_"))[1]; area3 = unlist(strsplit(analyses[i+2],"_"))[3]
		localTreesDirectory3 = paste0("Tree_extractions/LASV2_",area3,"_seg",segment3)
		localTreesDirectory = gsub("NIG1","NIG",localTreesDirectory1)
		dir.create(file.path(localTreesDirectory), showWarnings=F)
		for (j in 1:nberOfExtractionFiles)
			{
				tab1 = read.csv(paste0(localTreesDirectory1,"/TreeExtractions_",j,".csv"), header=T)
				tab2 = read.csv(paste0(localTreesDirectory2,"/TreeExtractions_",j,".csv"), header=T)
				tab3 = read.csv(paste0(localTreesDirectory3,"/TreeExtractions_",j,".csv"), header=T)
				genotypes1 = matrix(1, nrow=dim(tab1), ncol=1); colnames(genotypes1) = "genotype"
				genotypes2 = matrix(2, nrow=dim(tab2), ncol=1); colnames(genotypes2) = "genotype"
				genotypes3 = matrix(3, nrow=dim(tab3), ncol=1); colnames(genotypes3) = "genotype"
				tab1 = cbind(tab1, genotypes1); tab2 = cbind(tab2, genotypes2); tab3 = cbind(tab3, genotypes3)
				tab = tab1; maxNodeID = max(tab[,c("node1","node2")])
				tab2[,c("node1","node2")] = tab2[,c("node1","node2")]+maxNodeID
				tab = rbind(tab,tab2); maxNodeID = max(tab[,c("node1","node2")])
				tab3[,c("node1","node2")] = tab3[,c("node1","node2")]+maxNodeID; tab = rbind(tab,tab3)
				write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
			}
	}
	
	# 4.3. Generating and saving annual progression polygons (convex hull polygons)

for (i in 1:length(study_areas))
	{
		for (j in 1:length(segments))
			{
				points = c(); years = c(1900:2019)
				localTreesDirectory = paste0("Tree_extractions/LASV2_",study_areas[i],"_seg",segments[j])
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
								polygonName = paste0("LASV2_",study_areas[i],"_seg",segments[j],"_",years[y])
								writeOGR(contourPolygons_df, dsn="Minimum_convex_hulls", layer=polygonName, driver="ESRI Shapefile", overwrite_layer=T)
							}
					}
			}
	}
	
# 5. Generating dispersal history graphs

setwd(paste(wd,wd5,sep="/"))
percentage = 80; precision = 20
analyses = c("MRU_segL","NIG1_segL","NIG2_segL","NIG3_segL","MRU_segS","NIG1_segS","NIG2_segS","NIG3_segS")

	# 5.1. Estimating kernel density polygons
	
for (a in 1:length(analyses))
	{
		localTreesDirectory = paste0("Tree_extractions/LASV2_",analyses[a])
		area = unlist(strsplit(analyses[a],"_"))[1]
		segment = unlist(strsplit(analyses[a],"_seg"))[2]
		mcc = read.csv(paste0("LASV_",segment,"_3_",area,"_MCC.csv"))
		startDatum = min(mcc[,"startYear"]); prob = percentage/100
		polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
		for (i in 1:length(polygons))
			{
				year = as.numeric(names(polygons[[i]]))
				layerName = paste0("LASV_",area,"_",segment,"_",percentage,"_",precision,"_",year)
				writeOGR(polygons[[i]], dsn=paste0("./LASV_",percentage,"HPD_",precision,"_years"), layer=layerName, driver="ESRI Shapefile")
			}
	}

	# 5.2. Plotting the global dispersal history graphs
	
background1 = raster("Environmental_files/LASV_rasters/Elevation_MRU_0.008.asc")
background1 = crop(background1, e_MRU); background1[!is.na(background1[])] = 1; r = background1
background_cols1 = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
borders1 = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="International_borders"), e_MRU)
lakes1 = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), e_MRU)
rivers1 = crop(readOGR(dsn="Environmental_files/LASV_rivers/", layer="FAO_GeoNetwork"), e_MRU)
rivers1@data[which(rivers1@data[,"Strahler"]<3),"Strahler"] = 0
riverCols1 = rep("#AFBED5",dim(rivers1@data)[1]); riverCols1[which(rivers1@data[,"Strahler"]<3)] = NA
background2 = raster("Environmental_files/LASV_rasters/Elevation_NGA_0.008.asc")
background2 = crop(background2, e_NGA); background2[!is.na(background2[])] = 1; r = background2
background_cols2 = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]
borders2 = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="International_borders"), e_NGA)
lakes2 = crop(readOGR(dsn="Environmental_files/Natural_Earth/", layer="Natural_Earth_lakes"), e_NGA)
rivers2 = crop(readOGR(dsn="Environmental_files/LASV_rivers/", layer="FAO_GeoNetwork"), e_NGA)
rivers2@data[which(rivers2@data[,"Strahler"]<3),"Strahler"] = 0
riverCols2 = rep("#AFBED5",dim(rivers2@data)[1]); riverCols2[which(rivers2@data[,"Strahler"]<3)] = NA
titles = c("MRU lineage","NGA lineage 1","NGA lineage 2","NGA lineage 3")
cols = colorRampPalette(brewer.pal(11,'RdYlGn'))(131)[21:121]
for (h in 1:2)
	{
		pdf(paste0("Dispersal_graphs_",segments[h],".pdf"), width=12, height=10.2) # dev.new(width=12, height=10.2)
		par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(1.5,2.7,1.0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o"); n = 0
		for (i in 1:4)
			{
				if (i == 1)
					{
						background = background1; background_cols = background_cols1
						borders = borders1; lakes = lakes1; rivers = rivers1; riverCols = riverCols1
					}	else	{
						background = background2; background_cols = background_cols2
						borders = borders2; lakes = lakes2; rivers = rivers2; riverCols = riverCols2
					}
				a = ((h-1)*4)+i
				area = unlist(strsplit(analyses[a],"_"))[1]
				segment = unlist(strsplit(analyses[a],"_seg"))[2]
				mcc = read.csv(paste0("LASV_",segment,"_3_",area,"_MCC.csv"))
				minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
				ancestralC = (((mcc[1,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				startYearsC = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				endYearsC = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
				col_ancestral = cols[ancestralC]; cols_endYears = cols[endYearsC]
				files = list.files(paste0("LASV_",percentage,"HPD_",precision,"_years")); pols = list(); cols_pol = list()
				files = files[which((grepl(area,files))&(grepl(paste0("_",segment,"_"),files))&(grepl(".shp",files)))]
				for (j in 1:length(files))
					{
						pol = shapefile(paste0("LASV_",percentage,"HPD_",precision,"yrs/",files[j]))
						pol@proj4string = CRS("+init=epsg:4326"); pols[[j]] = pol
						date = gsub(".shp","",unlist(strsplit(files[j],"_"))[length(unlist(strsplit(files[j],"_")))])
						date = as.numeric(date)+1-(precision/2)
						yearM = (((date-minYear)/(maxYear-minYear))*100)+1
						cols_pol[[j]] = cols[yearM]
					}
				for (j in 1:length(cols_pol)) cols_pol[[j]] = paste0(cols_pol[[j]],"70")
				plot(background, main="", cex.main=1, cex.axis=0.7, bty="n", box=F, axes=F, legend=F, axis.args=list(cex.axis=0.7), col="white", colNA="gray90")
				for (j in 1:length(pols)) plot(pols[[j]], axes=F, col=cols_pol[[j]], add=T, border=NA)
				plot(lakes, add=T, lwd=0.7, col=unique(riverCols)[!is.na(unique(riverCols))], border=NA)
				plot(rivers, add=T, lwd=rivers@data[,"Strahler"]/3.5, col=riverCols)
				if (i != 1)
					{
						mtext(expression(bold(Nigeria)), side=3, line=-8.8, at=6.85, cex=0.9, font=1, col="gray30")
					}
				mtext(titles[i], side=3, line=-27.3, at=3.7, cex=0.9, font=1, col="gray30")
				plot(borders, add=T, lwd=2, col="gray60", lty=2)
				rect(xmin(background), ymin(background), xmax(background), ymax(background), xpd=T, lwd=0.2)
				axis(1, c(ceiling(xmin(background)), floor(xmax(background))), pos=ymin(background), mgp=c(0,0.7,0), cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
				axis(2, c(ceiling(ymin(background)), floor(ymax(background))), pos=xmin(background), mgp=c(0,0.8,0), cex.axis=0.9, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
				for (j in 1:dim(mcc)[1])
					{
						curvedarrow(cbind(mcc[j,"startLon"],mcc[j,"startLat"]), cbind(mcc[j,"endLon"],mcc[j,"endLat"]), arr.length=0,
						    		arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
					}
				for (j in dim(mcc)[1]:1)
					{
						points(mcc[j,"endLon"], mcc[j,"endLat"], pch=16, col=cols_endYears[j], cex=1.0)
						points(mcc[j,"endLon"], mcc[j,"endLat"], pch=1, col="gray30", cex=1.0, lwd=0.1)
						if (j == 1)
							{
								points(mcc[j,"startLon"], mcc[j,"startLat"], pch=16, col=col_ancestral, cex=1.0)
								points(mcc[j,"startLon"], mcc[j,"startLat"], pch=1, col="gray30", cex=1.0, lwd=0.1)
							}
					}
				legendRast = raster(as.matrix(c(minYear,maxYear)))
				plot(legendRast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.810,0.825,0.09,0.56), alpha=1.0,
					 legend.args=list(text="", cex=0.9, line=0.5, col="gray30"), axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.2, tck=-0.8, col.axis="gray30", line=0, mgp=c(0,0.6,0)))
				if (i == 1)
					{
						legend(x=-14.5, y=6.3, legend=c("S = 3 ","S = 4 ","S = 5 ","S = 6 ","S = 7 "), lwd=c(3:7)/3.5, cex=0.9, lty=1,
							   col=unique(riverCols)[!is.na(unique(riverCols))], x.intersp=0.5, border=NA, bg="white", box.lwd=0.2, box.col="gray30")
			   		}
			}
		dev.off()
	}

# 6. Estimating dispersal statistics

setwd(paste(wd,wd6,sep="/"))
analyses = c("MRU_segL","NIG_segL","NIG1_segL","NIG2_segL","NIG3_segL","MRU_segS","NIG_segS","NIG1_segS","NIG2_segS","NIG3_segS")
timSlices = 200; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 1; slidingWindow = 5
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0("Tree_extractions/LASV2_",analyses[i]); outputName = analyses[i]
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	}

rounds = c(2,1); genotypes = c("II","III"); n = 0
results = matrix(nrow=length(analyses), ncol=length(statistics))
colNames = c()
for (i in 1:length(analyses))
	{
		tab = read.table(paste0("Dispersal_stats_LASV_data_sets/",analyses[i],"_estimated_statistics.txt"), header=T)
		for (j in 1:length(statistics))
			{
				median = round(median(tab[,statistics[j]]),rounds[j])
				qtiles = round(quantile(tab[,statistics[j]],c(0.025,0.975)),rounds[j])
				results[i,j] = as.character(paste0(median," [",qtiles[1],"-",qtiles[2],"]"))
			}
	}
colNames = c("Weighted dispersal velocity (km/year)","Weighted diffusion coefficient (km2/year)")
colnames(results) = colNames; row.names(results) = analyses
write.table(results, "Dispersal_statistics.txt", sep="	", quote=F)

pdf("Spatial_wavefronts_NEW.pdf", width=7, height=3); datasets = c(); # dev.new(width=7, height=3)
par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(1,1.2,0,0), mar=c(1.5,1.5,1,1)); cols1 = list(); cols2 = list()
cols1[[1]] = rgb(233, 157, 34, 255, maxColorValue=255); cols2[[1]] = rgb(233, 157, 34, 100, maxColorValue=255); cols1[[3]] = cols1[[1]]; cols2[[3]] = cols2[[1]]
cols1[[2]] = rgb(43, 153, 43, 255, maxColorValue=255); cols2[[2]] = rgb(43, 153, 43, 100, maxColorValue=255); cols1[[4]] = cols1[[2]]; cols2[[4]] = cols2[[2]]
genotype_names = c("genotype II","genotype III"); segment_names = c("Segment L","Segment L","Segment S","Segment S")
for (i in 1:4)
	{
		tab1 = read.table(paste0("Dispersal_stats_LASV_data_sets/",analyses[i],"_median_spatial_wavefront_distance.txt"), header=T)
		tab2 = read.table(paste0("Dispersal_stats_LASV_data_sets/",analyses[i],"_95%HPD_spatial_wavefront_distance.txt"), header=T)
		if (i%in%c(1,3)) plot(tab1[,1], tab1[,2], type="l", axes=F, ann=F, ylim=c(0,700), xlim=c(1200,2018), col=NA)
		slicedTimes = tab1[,1]; waveFrontDistances1MedianValue = tab1[,2]; lower_l_1 = tab2[,2]; upper_l_1 = tab2[,3]
		xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
		getOption("scipen"); opt = options("scipen"=20); polygon(xx_l, yy_l, col=cols2[[i]], border=0)
		lines(slicedTimes, waveFrontDistances1MedianValue, lwd=1, col=cols1[[i]], xlim=c(1400,2018))
		axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.02,0), lwd=0, tck=-0.020, col.axis="gray30")
		axis(side=2, pos=1200, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.2,0), lwd=0, tck=-0.015, col.axis="gray30")
		title(xlab="time (year)", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
		title(ylab="distance (km)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
		# title(main="Furthest extent of epidemic wavefront (spatial distance from epidemic origin)", cex.main=0.55, col.main="gray30", line=0.3)
		rect(1200, 0, 2018, 700, lwd=0.2, border="gray30"); # box(lwd=0.2, col="gray30")
		mtext(side=3, segment_names[i], line=-2.2, at=1360, cex=0.65, col="gray30")
		legend(1265, 600, legend=genotype_names, pch=16, col=unlist(cols1), border="gray30", text.col="gray30", bty="n", pt.cex=1, cex=0.6, x.intersp=0.9, y.intersp=1.2)
	}
dev.off()

# 7. Post-hoc landscape phylogeography

setwd(paste(wd,wd7,sep="/"))

	# 7.1. Generating null models of lineages dispersal

		# 7.1.1. Randomisation of tree branches position (randomisation procedure n°3)

source("Seraphim_functions/treesRandomisation.r")
for (i in 1:length(analyses))
	{
		if (i <= 4) rast = raster("Environmental_files/LASV_rasters/Elevation_NGA_0.008.asc")
		if (i >= 5) rast = raster("Environmental_files/LASV_rasters/Elevation_MRU_0.008.asc")
		rast[!is.na(rast[])] = 0
		localTreesDirectory = localTreesDirectories[i]; envVariables = list(rast); randomProcedure = 3; nberOfCores = 1
		treesRandomisation(localTreesDirectory, nberOfExtractionFiles, envVariables, randomProcedure, nberOfCores)
	}
for (i in c(1,3))
	{
		localTreesDirectory1 = localTreesDirectories[i]
		localTreesDirectory2 = localTreesDirectories[i+1]
		localTreesDirectory = gsub("genII","allG",localTreesDirectories[i])
		for (j in 1:nberOfExtractionFiles)
			{
				tab1 = read.csv(paste0(localTreesDirectory1,"/TreeRandomisation_",j,".csv"), header=T)
				tab2 = read.csv(paste0(localTreesDirectory2,"/TreeRandomisation_",j,".csv"), header=T)
				genotypes2 = matrix(2, nrow=dim(tab1), ncol=1); colnames(genotypes2) = "genotype"
				genotypes3 = matrix(3, nrow=dim(tab2), ncol=1); colnames(genotypes3) = "genotype"
				tab1 = cbind(tab1, genotypes2); tab2 = cbind(tab2, genotypes3)
				maxNodeID = max(tab1[,c("node1","node2")])
				tab2[,c("node1","node2")] = tab2[,c("node1","node2")]+maxNodeID
				tab = rbind(tab1,tab2)
				write.csv(tab, paste0(localTreesDirectory,"/TreeRandomisation_",j,".csv"), row.names=F, quote=F)
			}
	}

		# 7.1.2. Relaxed random walk simulations along posterior trees

model = "gamma"; showingPlots = FALSE; nodesOnly = FALSE; newPlot = FALSE; n1 = 100; n2 = 0
rast1 = raster("Environmental_files/LASV_rasters/Elevation_NGA_0.008.asc")
rast1[!is.na(rast1[])] = 0; test = rast1; test[] = 0
for (i in 1:length(analyses))
	{
		if (i <= 4) rast1 = raster("Environmental_files/LASV_rasters/Elevation_NGA_0.008.asc")
		if (i >= 5) rast1 = raster("Environmental_files/LASV_rasters/Elevation_MRU_0.008.asc")
		rast1[!is.na(rast1[])] = 0; test[] = 0; points = c()
		if (i <= 4) analysesToConsider = c(1:4)
		if (i >= 5) analysesToConsider = c(5:6)
		for (j in analysesToConsider)
			{
				points1 = c(); points2 = c()
				localTreesDirectory = localTreesDirectories[j]
				for (k in 1:nberOfExtractionFiles)
					{
						if (j != 5)
							{
								tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",k,".csv"), header=T)
							}	else	{
								tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",k,"_with_root.csv"), header=T)
							}
						points1 = rbind(points1, tab[,c("startLon","startLat")])
						points2 = rbind(points2, tab[,c("endLon","endLat")])
					}
				colnames(points1) = c("lon","lat"); colnames(points2) = c("lon","lat")
				points = rbind(points, points1, points2)
			}
		hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,])
		ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pointsRaster = rasterize(points, crop(rast1, sps, snap="out"))
		pointsRaster[!is.na(pointsRaster[])] = 0; # plot(mask(simRasters[[h]],sps))
		hullRaster = crop(rast1, sps, snap="out"); bufferRaster = hullRaster
		rast2 = mask(hullRaster, sps, snap="out")
		rast2[!is.na(pointsRaster[])] = bufferRaster[!is.na(pointsRaster[])]
		envVariables = list(rast2); # envVariables = list(test)
		localTreesDirectory = localTreesDirectories[i]
		mostRecentSamplingDatum = mostRecentSamplingDates[i]
		trees = readAnnotatedNexus(gsub(".trees",paste0("_",nberOfExtractionFiles,".trees"),treesFiles[i]))
		log = read.table(gsub(".trees",paste0("_",nberOfExtractionFiles,".log"),treesFiles[i]), header=T)
		for (j in 1:nberOfExtractionFiles) # for (j in 1:100)
			{
				tree = trees[[j]]
				if (i != 5)
					{
						fixingRootBranches = FALSE
						tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T)
						ancestID = which(!tab[,"node1"]%in%tab[,"node2"])[1]
						ancestPosition = c(tab[ancestID,"startLon"], tab[ancestID,"startLat"])
					}	else	{
						fixingRootBranches = TRUE
						tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,"_with_root.csv"), header=T)
						ancestIDs = which(!tab[,"node1"]%in%tab[,"node2"])
						ancestPosition = c(mean(tab[ancestIDs[1],"endLon"],tab[ancestIDs[2],"endLon"]), mean(tab[ancestIDs[1],"endLat"],tab[ancestIDs[2],"endLat"]))
					}
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
				source("Seraphim_functions/simulatorRRW1.r"); # showingPlots = TRUE; newPlot = TRUE
				output = simulatorRRW1(tree, rates, sigmas, cor, envVariables, mostRecentSamplingDatum,
									   ancestPosition, reciprocalRates, n1, n2, showingPlots, newPlot, fixingRootBranches)
				# pdf(paste0(analyses[i],"_treeSimulations_",j,".pdf"))
				# plot(density(geoDists[,1])); lines(density(output[,"greatCircleDist_km"]), col="red")
				# dev.off()
				file = as.character(paste(simulationsDirectory,"/TreeSimulations_",j,".csv",sep=""))
				write.csv(output, file, row.names=F, quote=F); sim = output
				# tab = tab[order(tab[,"startYear"],tab[,"endYear"]),]; sim = sim[order(sim[,"startYear"],sim[,"endYear"]),]
				# indices1 = which(!tab[,"node1"]%in%tab[,"node2"]); indices2 = which(!sim[,"node1"]%in%sim[,"node2"])
				# print(round(tab[indices1[1],"endLon"]-sim[indices2[1],"endLon"],5))
				# print(round(tab[indices1[1],"endLat"]-sim[indices2[1],"endLat"],5))
			}
	}
for (i in c(1,3))
	{
		localTreesDirectory1 = localTreesDirectories[i]
		localTreesDirectory2 = localTreesDirectories[i+1]
		localTreesDirectory = gsub("genII","allG",localTreesDirectories[i])
		for (j in 1:nberOfExtractionFiles)
			{
				tab1 = read.csv(paste0(localTreesDirectory1,"/TreeSimulations_",j,".csv"), header=T)
				tab2 = read.csv(paste0(localTreesDirectory2,"/TreeSimulations_",j,".csv"), header=T)
				genotypes2 = matrix(2, nrow=dim(tab1), ncol=1); colnames(genotypes2) = "genotype"
				genotypes3 = matrix(3, nrow=dim(tab2), ncol=1); colnames(genotypes3) = "genotype"
				tab1 = cbind(tab1, genotypes2); tab2 = cbind(tab2, genotypes3)
				maxNodeID = max(tab1[,c("node1","node2")])
				tab2[,c("node1","node2")] = tab2[,c("node1","node2")]+maxNodeID
				tab = rbind(tab1,tab2)
				write.csv(tab, paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), row.names=F, quote=F)
			}
	}
for (i in c(5))
	{
		localTreesDirectory = localTreesDirectories[i]
		for (j in 1:nberOfExtractionFiles)
			{
				tab1 = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), header=T)
				indices = which(!tab1[,"node1"]%in%tab1[,"node2"]); tab2 = tab1
				if (length(indices) == 2)
					{
						tab2 = tab2[which(tab1[,"node1"]%in%tab1[,"node2"]),]
						write.csv(tab1, paste0(localTreesDirectory,"/TreeSimulations_",j,"_with_root.csv"), row.names=F, quote=F)
						write.csv(tab2, paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), row.names=F, quote=F)
					}
			}
	}	
for (i in c(1,3,5,6)) # check-up:
	{
		if (i <= 4) localTreesDirectory = gsub("genII","allG",localTreesDirectories[i])
		if (i >= 5) localTreesDirectory = localTreesDirectories[i]
		dists = matrix(nrow=nberOfExtractionFiles, ncol=2)
		# localTreesDirectory = localTreesDirectories[i]; dists = matrix(nrow=100, ncol=2)
		for (j in 1:nberOfExtractionFiles) # for (j in 1:100)
			{
				obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T)
				sim = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), header=T)
				obs = obs[order(obs[,"node2"]),]; geoDists_obs = matrix(nrow=dim(obs)[1], ncol=1)
				sim = sim[order(sim[,"node2"]),]; geoDists_sim = matrix(nrow=dim(sim)[1], ncol=1)
				for (k in 1:dim(obs)[1])
					{
						x1 = cbind(obs[k,"startLon"], obs[k,"startLat"])
						x2 = cbind(obs[k,"endLon"], obs[k,"endLat"])
						geoDists_obs[k,1] = rdist.earth(x1, x2, miles=F, R=NULL)
						x1 = cbind(sim[k,"startLon"], sim[k,"startLat"])
						x2 = cbind(sim[k,"endLon"], sim[k,"endLat"])
						geoDists_sim[k,1] = rdist.earth(x1, x2, miles=F, R=NULL)
					}
				dists[j,1:2] = cbind(sum(geoDists_obs[,1]), sum(geoDists_sim[,1]))
				# plot(density(geoDists_obs[,1])); lines(density(geoDists_sim[,1]), col="red")
				# plot(geoDists_obs,geoDists_sim)
			}
		dev.new(width=6, height=2.5); par(mfrow=c(1,1), mgp=c(1,0.35,0), oma=c(0.5,1,1,1), mar=c(2.5,2,0.5,0), lwd=0.2)
		if (i == 1) { xLim = c(0,8000); yLim=c(0,0.0025); X = 1000; Y = 0.0025 }
		if (i == 3) { xLim = c(0,12000); yLim=c(0,0.0030); X = 1500; Y = 0.003 }
		if (i == 5) { xLim = c(0,4000); yLim=c(0,0.0075); X = 500; Y = 0.0075 }
		if (i == 6) { xLim = c(0,6500); yLim=c(0,0.0035); X = 700; Y = 0.0035 }
		cols1 = list(); cols1[[1]] = rgb(204,0,0,255,maxColorValue=255); cols1[[2]] = rgb(120,120,120,255,maxColorValue=255)
		cols2 = list(); cols2[[1]] = rgb(204,0,0,100,maxColorValue=255); cols2[[2]] = rgb(120,120,120,100,maxColorValue=255)
		plot(density(dists[,1]), lwd=0.7, col=cols1[[1]], ylim=yLim, xlim=xLim, axes=F, ann=F)
		lines(density(dists[,2]), lwd=0.7, col=cols1[[2]])
		polygon(density(dists[,1]), col=cols2[[1]], border=NA)
		polygon(density(dists[,2]), col=cols2[[2]], border=NA)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,-0.02,0), lwd=0.2, tck=-0.025, col.axis="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.025, col.axis="gray30")
		title(xlab="total great-circle distance travelled by branches (km)", cex.lab=0.7, mgp=c(0.8,0,0), col.lab="gray30")
		title(ylab="density", cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		legend(x=X, y=Y, legend=c("inferred trees","simulated trees"), lwd=0.7, cex=0.7,
			   col=c(unlist(cols1)), text.col=c(unlist(cols1)), border=NA, x.intersp=0.5, bty="n")
		if (i <= 4) dev.copy2pdf(file=paste0(gsub("genII","RRW1",analyses[i]),".pdf"))
		if (i >= 5) dev.copy2pdf(file=paste0(analyses[i],"_RRW1.pdf"))
	}

	# 7.2. Analysing the impact of rivers on the dispersal frequency of lineages

onInternalBranches = FALSE; # onInternalBranches = TRUE
for (i in c(1,3,5,6))
	{
		if (i <= 4) rast = raster("Environmental_files/LASV_rasters/Elevation_NGA_0.008.asc")
		if (i >= 5) rast = raster("Environmental_files/LASV_rasters/Elevation_MRU_0.008.asc")
		rast[!is.na(rast[])] = 0
		rivers = crop(readOGR(dsn="Environmental_files/LASV_rivers/", layer="FAO_GeoNetwork"), extent(rast))
		rivers_list = list()
		suffix = 1
		if (onInternalBranches == TRUE) suffix = 2
		if (i <= 4) strahlerNumbers = c(1:7)
		if (i >= 5) strahlerNumbers = c(1:5)
		for (j in 1:length(strahlerNumbers))
			{
				indices = which(rivers@data[,"Strahler"] >= strahlerNumbers[j])
				rivers_list[[j]] = SpatialLines(rivers@lines[indices])
				rivers_list[[j]]@proj4string = rivers@proj4string
			}
		for (j in 1:length(strahlerNumbers)) print(length(rivers_list[[j]]))
		if (i <= 4)
			{
				localTreesDirectory = gsub("genII","allG",localTreesDirectories[i])
				if (!file.exists(paste0(gsub("_genII","_obs",analyses[i]),".csv")))
					{
						crossingRiversObs = matrix(nrow=nberOfExtractionFiles, ncol=length(strahlerNumbers)); colnames(crossingRiversObs) = strahlerNumbers
						crossingRiversSim = matrix(nrow=nberOfExtractionFiles, ncol=length(strahlerNumbers)); colnames(crossingRiversSim) = strahlerNumbers
						crossingRiversRan = matrix(nrow=nberOfExtractionFiles, ncol=length(strahlerNumbers)); colnames(crossingRiversRan) = strahlerNumbers
						write.csv(crossingRiversObs, paste0(gsub("_genII","_obs",analyses[i]),suffix,".csv"), row.names=F, quote=F)
						write.csv(crossingRiversSim, paste0(gsub("_genII","_sim",analyses[i]),suffix,".csv"), row.names=F, quote=F)
						write.csv(crossingRiversRan, paste0(gsub("_genII","_ran",analyses[i]),suffix,".csv"), row.names=F, quote=F)
					}	else		{
						crossingRiversObs = read.csv(paste0(gsub("_genII","_obs",analyses[i]),suffix,".csv"), header=T)
						crossingRiversSim = read.csv(paste0(gsub("_genII","_sim",analyses[i]),suffix,".csv"), header=T)
						crossingRiversRan = read.csv(paste0(gsub("_genII","_ran",analyses[i]),suffix,".csv"), header=T)
					}
			}
		if (i >= 5)
			{
				localTreesDirectory = localTreesDirectories[i]
				if (!file.exists(paste0(analyses[i],"_obs.csv")))
					{
						crossingRiversObs = matrix(nrow=nberOfExtractionFiles, ncol=length(strahlerNumbers)); colnames(crossingRiversObs) = strahlerNumbers
						crossingRiversSim = matrix(nrow=nberOfExtractionFiles, ncol=length(strahlerNumbers)); colnames(crossingRiversSim) = strahlerNumbers
						crossingRiversRan = matrix(nrow=nberOfExtractionFiles, ncol=length(strahlerNumbers)); colnames(crossingRiversRan) = strahlerNumbers
						write.csv(crossingRiversObs, paste0(analyses[i],"_obs",suffix,".csv"), row.names=F, quote=F)
						write.csv(crossingRiversSim, paste0(analyses[i],"_sim",suffix,".csv"), row.names=F, quote=F)
						write.csv(crossingRiversRan, paste0(analyses[i],"_ran",suffix,".csv"), row.names=F, quote=F)
					}	else		{
						crossingRiversObs = read.csv(paste0(analyses[i],"_obs",suffix,".csv"), header=T)
						crossingRiversSim = read.csv(paste0(analyses[i],"_sim",suffix,".csv"), header=T)
						crossingRiversRan = read.csv(paste0(analyses[i],"_ran",suffix,".csv"), header=T)
					}				
			}
		for (j in 1:nberOfExtractionFiles)
			{
				if (sum(is.na(crossingRiversObs[j,])) == length(strahlerNumbers))
					{
						crossingRiversObsLine = matrix(0, nrow=1, ncol=length(strahlerNumbers)); colnames(crossingRiversObsLine) = strahlerNumbers
						crossingRiversSimLine = matrix(0, nrow=1, ncol=length(strahlerNumbers)); colnames(crossingRiversSimLine) = strahlerNumbers
						crossingRiversRanLine = matrix(0, nrow=1, ncol=length(strahlerNumbers)); colnames(crossingRiversRanLine) = strahlerNumbers
						obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",j,".csv"), header=T); branches_obs = list()
						sim = read.csv(paste0(localTreesDirectory,"/TreeSimulations_",j,".csv"), header=T); branches_sim = list()
						ran = read.csv(paste0(localTreesDirectory,"/TreeRandomisation_",j,".csv"), header=T); branches_ran = list()
						if (onInternalBranches == TRUE)
							{
								obs = obs[which(obs[,"node2"]%in%obs[,"node1"]),]
								sim = sim[which(sim[,"node2"]%in%sim[,"node1"]),]
								ran = ran[which(ran[,"node2"]%in%ran[,"node1"]),]
							}
						for (k in 1:dim(obs)[1])
							{
								x = c(obs[k,"startLon"],obs[k,"endLon"]); y = c(obs[k,"startLat"],obs[k,"endLat"])
								branches_obs[[k]] = Lines(Line(cbind(x,y)), ID=k)
							}
						for (k in 1:dim(sim)[1])
							{	
								x = c(sim[k,"startLon"],sim[k,"endLon"]); y = c(sim[k,"startLat"],sim[k,"endLat"])
								branches_sim[[k]] = Lines(Line(cbind(x,y)), ID=k)
							}
						for (k in 1:dim(ran)[1])
							{
								x = c(ran[k,"startLon"],ran[k,"endLon"]); y = c(ran[k,"startLat"],ran[k,"endLat"])
								branches_ran[[k]] = Lines(Line(cbind(x,y)), ID=k)
							}
						branches_obs = SpatialLines(branches_obs); crs(branches_obs) = rivers@proj4string
						branches_sim = SpatialLines(branches_sim); crs(branches_sim) = rivers@proj4string
						branches_ran = SpatialLines(branches_ran); crs(branches_ran) = rivers@proj4string
						# plot(branches_obs, col="gray30"); plot(branches_sim, add=T, col="red"); plot(branches_ran, add=T, col="green3")
						# k=4; plot(rivers_list[[k]], col="blue", add=T)
						for (k in 1:length(strahlerNumbers))
							{
								crossingRiversObs[j,k] = length(gIntersection(branches_obs, rivers_list[[k]]))
								crossingRiversSim[j,k] = length(gIntersection(branches_sim, rivers_list[[k]]))
								crossingRiversRan[j,k] = length(gIntersection(branches_ran, rivers_list[[k]]))
							}
						print(rbind(crossingRiversObs[j,],crossingRiversSim[j,],crossingRiversRan[j,]))
						if (i <= 4)
							{
								write.csv(crossingRiversObs, paste0(gsub("_genII","_obs",analyses[i]),suffix,".csv"), row.names=F, quote=F)
								write.csv(crossingRiversSim, paste0(gsub("_genII","_sim",analyses[i]),suffix,".csv"), row.names=F, quote=F)
								write.csv(crossingRiversRan, paste0(gsub("_genII","_ran",analyses[i]),suffix,".csv"), row.names=F, quote=F)
							}
						if (i >= 5)
							{
								write.csv(crossingRiversObs, paste0(analyses[i],"_obs",suffix,".csv"), row.names=F, quote=F)
								write.csv(crossingRiversSim, paste0(analyses[i],"_sim",suffix,".csv"), row.names=F, quote=F)
								write.csv(crossingRiversRan, paste0(analyses[i],"_ran",suffix,".csv"), row.names=F, quote=F)
							}
					}
			}
	}
dS = list(); BFs = dS; c = 0; suffix = 1; onInternalBranches = FALSE; # onInternalBranches = TRUE; suffix = 2
strahlerNumbers = c(1:7); datasets = c("LASV2_MRU_segL","LASV2_MRU_segS","LASV2_NGA_segL","LASV2_NGA_segS")
river_BFs1 = matrix(nrow=length(strahlerNumbers), ncol=4); row.names(river_BFs1) = strahlerNumbers; colnames(river_BFs1) = datasets
river_BFs2 = matrix(nrow=length(strahlerNumbers), ncol=4); row.names(river_BFs2) = strahlerNumbers; colnames(river_BFs2) = datasets
for (i in 1:length(datasets))
	{
		crossingRiversObs = read.csv(paste0(datasets[i],"_obs",suffix,".csv"), header=T)
		crossingRiversSim = read.csv(paste0(datasets[i],"_sim",suffix,".csv"), header=T)
		crossingRiversRan = read.csv(paste0(datasets[i],"_ran",suffix,".csv"), header=T)				
		# s = 7; plot(density(crossingRiversObs[,s]), xlim=c(1,1000), frame=F, col="gray30")
		# lines(density(crossingRiversSim[,s]), col="red")
		for (j in 1:dim(crossingRiversObs)[2])
			{
				c1 = 0; c2 = 0; n = 0
				for (k in 1:dim(crossingRiversObs)[1])
					{
						if ((!is.na(sum(crossingRiversObs[k,])))&(!is.na(sum(crossingRiversSim[k,])))&(!is.na(sum(crossingRiversRan[k,]))))
							{
								n = n+1
								if (sum(crossingRiversObs[k,j]) < sum(crossingRiversSim[k,j])) c1 = c1+1
								if (sum(crossingRiversObs[k,j]) < sum(crossingRiversRan[k,j])) c2 = c2+1
							}
					}
				p1 = c1/n; river_BFs1[j,i] = p1/(1-p1); p2 = c2/n; river_BFs2[j,i] = p2/(1-p2)
				if (j > 2)
					{
						c = c + 1; dS[[c]] = crossingRiversSim[,j]-crossingRiversObs[,j]; BFs[[c]] = p1/(1-p1)
					}
			}
	}
riverCol = "#AFBED5"; minX = 9999; maxX = -9999
for (i in 1:length(dS))
	{
		if (minX > min(dS[[i]])) minX = min(dS[[i]])
		if (maxX < max(dS[[i]])) maxX = max(dS[[i]])	
	}
dev.new(width=8, height=2.5); # pdf("Rivers_impact_NEW.pdf", width=10, height=4)
par(mfrow=c(1,2), oma=c(0,2,0,0), mar=c(2.5,3.0,0.5,2.7), mgp=c(0,0.1,0), bty="n")
labels = rev(c("MRU, S > 2","MRU, S > 3","MRU, S > 4","Nigeria, S > 2","Nigeria, S > 3","Nigeria, S > 4","Nigeria, S > 5","Nigeria, S > 6"))
plot(0:1, 0:1, type="n", xlim=c(-60,90), ylim=c(0,9), axes=F, ann=F)
vioplot(dS[[8]],dS[[7]],dS[[6]],dS[[5]],dS[[4]],dS[[3]],dS[[2]],dS[[1]],
		horizontal=T, drawRect=F, col=riverCol, add=T, wex=1, outline=F, box=F, lwd=0.2)
axis(side=2, at=1:8, labels=labels, las=1, cex.axis=0.6, lwd=0, lwd.ticks=0, col.axis="gray30", mgp=c(0,0,0))
axis(side=4, at=1:8, labels=c(rev(round(unlist(BFs[1:8]),1))), las=1, cex.axis=0.6, lwd=0, lwd.ticks=0, col.axis="red", mgp=c(0,0,0))
axis(side=1, cex.axis=0.6, lwd=0.2, lwd.ticks=0.17, col="gray30", col.axis="gray30", pos=0, tck=-0.025, mgp=c(0,0.0,0))
mtext("D (segment L)", side=1, col="gray30", cex=0.7, line=0.5)
mtext("BF", at=100, side=1, col="red", cex=0.7, line=-0.8)
segments(0, 0, 0, 14.5, col="gray30", lwd=0.2, lty=2)
plot(0:1, 0:1, type="n", xlim=c(-60,120), ylim=c(0,9), axes=F, ann=F)
vioplot(dS[[16]],dS[[15]],dS[[14]],dS[[13]],dS[[12]],dS[[11]],dS[[10]],dS[[9]],
		horizontal=T, drawRect=F, col=riverCol, add=T, wex=1, outline=F, box=F, lwd=0.2)
axis(side=2, at=1:8, labels=labels, las=1, cex.axis=0.6, lwd=0, lwd.ticks=0, col.axis="gray30", mgp=c(0,0,0))
axis(side=4, at=1:8, labels=c(rev(round(unlist(BFs[9:16]),1))), line=0, las=1, cex.axis=0.6, lwd=0, lwd.ticks=0, col.axis="red", mgp=c(0,0,0))
axis(side=1, cex.axis=0.6, lwd=0.2, lwd.ticks=0.17, col="gray30", col.axis="gray30", pos=0, tck=-0.025, mgp=c(0,0.0,0))
mtext("D (segment S)", side=1, col="gray30", cex=0.7, line=0.5)
mtext("BF", at=133, side=1, col="red", cex=0.7, line=-0.8)
segments(0, 0, 0, 14.5, col="gray30", lwd=0.2, lty=2); # dev.off()

	# 7.3. Analysing the impact of environmental factors on the dispersal velocity of lineages

		# 7.3.1. Analyses based on continuous environmental rasters

nberOfExtractionFiles = 1000; nberOfRandomisations = 0; randomProcedure = 3; fourCells = FALSE
showingPlots = FALSE; nberOfCores = 10; nberOfCores_CS = 1; OS = "Unix"
envVariableNames = c("Cover_land_croplands","Cover_land_forests","Cover_land_grasslands","Cover_land_savannas",
					  "Elevation","Annual_mean_temperature","Annual_precipitation")
for (i in c(1,3,5,6))
	{
		if (i <= 4)
			{
				analysis = gsub("_genII","",analyses[i])
				localTreesDirectory = gsub("genII","allG",localTreesDirectories[i])
				raster_names = paste0(envVariableNames,"_NGA_0.04")
			}
		if (i >= 5)
			{
				analysis = analyses[i]
				localTreesDirectory = localTreesDirectories[i]
				raster_names = paste0(envVariableNames,"_MRU_0.04")
			}
		envVariables = list(); resistances = list(); avgResistances = list(); c = 0
		for (k in c(10,100,1000))
			{
				for (j in 1:length(raster_names))
					{
						c = c+1
						rast = raster(paste("Environmental_files/LASV_rasters/",raster_names[j],".asc",sep=""))
						rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(raster_names[j], "_k", k, sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(raster_names[j],"_k",k,sep="")
						resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
					}
				for (j in 1:length(raster_names))
					{
						c = c+1
						rast = raster(paste("Environmental_files/LASV_rasters/",raster_names[j],".asc",sep=""))
						rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(raster_names[j], "_k", k, sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(raster_names[j],"_k",k,sep="")
						resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
					}
			}
		pathModel = 2; simulations = FALSE; outputName = paste0(analysis,"_seraphim_LC_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
		pathModel = 2; simulations = TRUE; outputName = paste0(analysis,"_seraphim_LC_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
		pathModel = 3; simulations = FALSE; outputName = paste0(analysis,"_seraphim_CS_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
		pathModel = 3; simulations = TRUE; outputName = paste0(analysis,"_seraphim_CS_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
	}

pathModels = c("Least-cost path model","Circuitscape path model")
envVariableNames = c("Cover_land_croplands","Cover_land_forests","Cover_land_grasslands","Cover_land_savannas",
					  "Elevation","Annual_mean_temperature","Annual_precipitation")
for (a in c(1,3,5,6))
	{
		if (a <= 4)
			{
				analysis = gsub("_genII","",analyses[a]); raster_names = paste0(envVariableNames,"_NGA_0.04")
			}
		if (a >= 5)
			{
				analysis = analyses[a]; raster_names = paste0(envVariableNames,"_MRU_0.04")
			}
		extractions = list(); simulations = list(); randomisations = list()
		extractions[[1]] = read.table(paste0(analysis,"_seraphim_LC_extractions_LR_results.txt"), header=T)
		extractions[[2]] = read.table(paste0(analysis,"_seraphim_CS_extractions_LR_results.txt"), header=T)
		simulations[[1]] = read.table(paste0(analysis,"_seraphim_LC_simulations_LR_results.txt"), header=T)
		simulations[[2]] = read.table(paste0(analysis,"_seraphim_CS_simulations_LR_results.txt"), header=T)
		allResults = matrix(nrow=length(envVariableNames)*2*2*3, ncol=7); kS = c(10,100,1000); CR = c("C","R"); L = 0
		colnames(allResults) = c("Path model","Environmental factor","k","Regression coefficient","Q statistic","p(Q) > 0","BF")
		for (i in 1:length(pathModels))
			{
				for (j in 1:length(envVariableNames))
					{
						for (k in 1:length(CR))
							{
								for (l in 1:length(kS))
									{
										L = L+1; c1 = 0; c2 = 0; allResults[L,1] = pathModels[i]; allResults[L,2] = paste0(envVariableNames[j]," (",CR[k],")"); allResults[L,3] = kS[l]
index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(envVariableNames[j],colnames(simulations[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
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
		cat("\n"); cat(analyses[a]); cat("\n"); print(allResults, quote=F)
	}
	
		# 7.3.2. Analyses based on rivers network rasters

regions = c("NGA","MRU"); extents = list(e_NGA,e_MRU); maxS = c(7,5)
for (i in 1:length(regions))
	{
		rivers = crop(readOGR(dsn="Environmental_files/LASV_rivers/", layer="FAO_GeoNetwork"), extents[[i]])
		template = raster(paste0("Environmental_files/LASV_rasters/Elevation_",regions[i],"_0.008.asc"))
		for (s in 3:maxS[i])
			{
				if (!file.exists(paste0("Environmental_files/LASV_rasters/Main_rivers_S",s,"_",regions[i],"_0.008.asc")))
					{
						selection = subset(rivers, rivers@data[,"Strahler"]>=s)
						rast = rasterize(selection, template)
						rast[!is.na(rast[])] = 1; rast[(is.na(rast[]))&(!is.na(template[]))] = 0
						writeRaster(rast, paste0("Environmental_files/LASV_rasters/Main_rivers_S",s,"_",regions[i],"_0.008.asc"))
					}
			}
	}

nberOfExtractionFiles = 1000; nberOfRandomisations = 0; randomProcedure = 3; fourCells = FALSE
showingPlots = FALSE; nberOfCores = 10; nberOfCores_CS = 1; OS = "Unix"
for (i in c(1,3,5,6))
	{
		if (i <= 4)
			{
				localTreesDirectory = gsub("genII","allG",localTreesDirectories[i]); analysis = gsub("_genII","",analyses[i])
				raster_names = c("Main_rivers_S3_NGA_0.008","Main_rivers_S4_NGA_0.008","Main_rivers_S5_NGA_0.008",
								 "Main_rivers_S6_NGA_0.008","Main_rivers_S7_NGA_0.008")
			}
		if (i >= 5)
			{
				localTreesDirectory = localTreesDirectories[i]; analysis = analyses[i]
				raster_names = c("Main_rivers_S3_MRU_0.008","Main_rivers_S4_MRU_0.008","Main_rivers_S5_MRU_0.008")
			}
		envVariables = list(); resistances = list(); avgResistances = list(); c = 0
		for (k in c(10,100,1000,10000))
			{
				for (j in 1:length(raster_names))
					{
						c = c+1
						rast = raster(paste("Environmental_files/LASV_rasters/",raster_names[j],".asc",sep=""))
						rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(raster_names[j], "_k", k, sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(raster_names[j],"_k",k,sep="")
						resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
					}
			}
		pathModel = 2; simulations = FALSE; outputName = paste0(analysis,"_rivers_LC_extractions")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
		pathModel = 2; simulations = TRUE; outputName = paste0(analysis,"_rivers_LC_simulations")
		spreadFactors(localTreesDirectory,nberOfExtractionFiles,envVariables,pathModel,resistances,avgResistances,fourCells,
					  nberOfRandomisations,randomProcedure,outputName,showingPlots,nberOfCores,OS,simulations)
	}

pathModels = c("Least-cost path model")
envVariableNames = c("Main_rivers_S3","Main_rivers_S4","Main_rivers_S5","Main_rivers_S6","Main_rivers_S7")
for (a in c(1,3,5,6))
	{
		if (a <= 4)
			{
				analysis = gsub("_genII","",analyses[a]); raster_names = paste0(envVariableNames,"_NGA_0.04")
			}
		if (a >= 5)
			{
				analysis = analyses[a]; raster_names = paste0(envVariableNames,"_MRU_0.04")
			}
		extractions = list(); simulations = list(); randomisations = list()
		extractions[[1]] = read.table(paste0(analysis,"_rivers_LC_extractions_LR_results.txt"), header=T)
		# extractions[[2]] = read.table(paste0(analysis,"_rivers_CS_extractions_LR_results.txt"), header=T)
		simulations[[1]] = read.table(paste0(analysis,"_rivers_LC_simulations_LR_results.txt"), header=T)
		# simulations[[2]] = read.table(paste0(analysis,"_rivers_CS_simulations_LR_results.txt"), header=T)
		allResults = matrix(nrow=length(envVariableNames)*3, ncol=7); kS = c(10,100,1000); CR = c("R"); L = 0
		colnames(allResults) = c("Path model","Environmental factor","k","Regression coefficient","Q statistic","p(Q) > 0","BF")
		for (i in 1:length(pathModels))
			{
				for (j in 1:length(envVariableNames))
					{
						if (sum(grepl(envVariableNames[j],colnames(extractions[[i]]))) > 1)
							{
								for (k in 1:length(CR))
									{
										for (l in 1:length(kS))
											{
												L = L+1; c1 = 0; c2 = 0; allResults[L,1] = pathModels[i]; allResults[L,2] = paste0(envVariableNames[j]," (",CR[k],")"); allResults[L,3] = kS[l]
index1 = which(grepl("LR_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
index2 = which(grepl("delta_R2",colnames(extractions[[i]]))&grepl(envVariableNames[j],colnames(extractions[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(extractions[[i]])))
index3 = which(grepl("delta_R2",colnames(simulations[[i]]))&grepl(envVariableNames[j],colnames(simulations[[i]]))&grepl(paste0("k",kS[l],"_",CR[k]),colnames(simulations[[i]])))
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
			}
		cat("\n"); cat(analyses[a]); cat("\n"); print(allResults, quote=F)
	}

# 8. Prepration of skygrid-GLM covariates

setwd(paste(wd,wd8,sep="/"))

	# 8.1. Preparation of climatic covariates
	
		# Preparation of annual temperature and precipitation estimates based on data from the Climate Research Unit (CRU)
		# Source: http://www.cru.uea.ac.uk (https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/cruts.1811131722.v4.02)

xMin = -16; xMax = 16; yMin = 3; yMax = 15
daysPerNormalYear = 365
daysInALeapYear = 366
cumulativeNberOfDays = c(0)
nberOfDays = 0; n = 0
for (i in 1:117)
	{
		n = n+1
		if (n == 4)
			{
				n = 0
				daysInTheYear = daysInALeapYear
			}	else 	{
				daysInTheYear = daysPerNormalYear
			}
		nberOfDays = nberOfDays+daysInTheYear
		cumulativeNberOfDays = c(cumulativeNberOfDays, nberOfDays)
	}
variables = c("temperature","precipitation"); suffixes = c("Temp","Prec")
yLabs = c("annual mean temperature","annual precipitation")
for (h in 1:length(variables))
	{
		if (!file.exists(paste0("CRU_",variables[h],".csv")))
			{
				nc = nc_open(paste0("CRU_",variables[h],".nc"))
				time = ncvar_get(nc,"time") # days since 1900-01-01 00:00:00.0
				lon = ncvar_get(nc,"lon"); lat = ncvar_get(nc,"lat")
				if (h == 1) var1 = ncvar_get(nc,"tmp")
				if (h == 2) var1 = ncvar_get(nc,"pre")
				if (h == 2) stns = ncvar_get(nc,"stn")
				var2 = c()
				indices1 = which((lon>xMin)&(lon<xMax))
				indices2 = which((lat>yMin)&(lat<yMax))
				lon = lon[indices1]; lat = lat[indices2]
				var1 = var1[indices1,indices2,]
				years = time; n = 0
				for (i in 1:length(years))
					{
						n = n+1
						if (n == 4)
							{
								n = 0
								daysInTheYear = daysInALeapYear
							}	else 	{
								daysInTheYear = daysPerNormalYear
							}
						index = which(cumulativeNberOfDays==max(cumulativeNberOfDays[cumulativeNberOfDays<=years[i]]))
						years[i] = 1899+index+((years[i]-cumulativeNberOfDays[index])/daysInTheYear)
					}
				n = 0; tot = dim(var1)[1]*dim(var1)[2]
				for (i in 1:dim(var1)[1])
					{
						for (j in 1:dim(var1)[2])
							{
								n = n+1; # print(paste0(n,"/",tot))
								row = cbind(lon[i],lat[j],t(var1[i,j,]))
								var2 = rbind(var2, row)
							}
					}
				colnames(var2) = c("longitude","latitude",years)
				write.csv(var2, paste0("CRU_",variables[h],".csv"), row.names=F, quote=F)
			}	else		{
				var2 = read.csv(paste0("CRU_",variables[h],".csv"), header=T)
			}
		curves1 = list(); curves2 = list(); years = c(1901:2017)
		analyses = c("LASV2_MRU_segL","LASV2_MRU_segS","LASV2_NIG_segL","LASV2_NIG_segS")
		for (i in 1:length(analyses))
			{
				values1 = matrix(nrow=length(years), ncol=2)
				for (y in 1:length(years))
					{
						buffer = var2[,3:dim(var2)[2]]; colNames = gsub("X","",colnames(buffer))
						selectedMeasures = cbind(var2[,1:2],buffer[,which(floor(as.numeric(colNames))==years[y])])
						hull = readOGR(dsn="Minimum_convex_hulls", layer=paste0(analyses[i],"_",years[y]))
						inHull = point.in.polygon(selectedMeasures[,"longitude"],selectedMeasures[,"latitude"],
								 		   		  hull@polygons[[1]]@Polygons[[1]]@coords[,1],
								 		   		  hull@polygons[[1]]@Polygons[[1]]@coords[,2])
						indices = which(inHull==1); values1[y,1] = years[y]
						values1[y,2] = mean(as.matrix(selectedMeasures[indices,3:dim(selectedMeasures)[2]]), na.rm=T)
					}
				curves1[[i]] = values1
				values2 = matrix(nrow=length(years), ncol=2)
				for (j in 1:length(years))
					{
						indices = which((values1[,1]>=(years[j]-2))&(values1[,1]<=(years[j]+2)))
						values2[j,1] = values1[j,1]
						values2[j,2] = mean(values1[indices,2])
					}
				values2[,1] = round(values2[,1],3); curves2[[i]] = values2
				colnames(values2) = c("time",variables[h])
				write.csv(values2, paste0(gsub("LASV2",suffixes[h],analyses[i]),".csv"), row.names=F)
			}
		if (h == 1)
			{
				dev.new(width=9, height=3); par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(1.0,1.0,0.5,0.5), mar=c(2.0,2.0,1,1))
				plot(curves1[[1]][,1], curves1[[1]][,2], lwd=0.5, type="l", axes=F, ann=F)
				lines(curves2[[1]][,1], curves2[[1]][,2], lwd=1.0, col="red")
				axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
				axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(24.4,26.6,0.2))
				title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
				title(ylab=yLabs[h], cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
				title(main="M-R-U", cex.main=0.8, col.main="gray30", line=-0.2)
				legend(1972, 24.95, legend=c("sliding window of 5 years"), lwd=1, cex=0.6, col=c("red"), text.col="gray30", border=NA, x.intersp=0.5, bty="n")
				plot(curves1[[3]][,1], curves1[[3]][,2], lwd=0.5, type="l", axes=F, ann=F)
				lines(curves2[[3]][,1], curves2[[3]][,2], lwd=1.0, col="red")
				axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
				axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(25.6,27.2,0.2))
				title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
				title(ylab=yLabs[h], cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
				title(main="Nigeria", cex.main=0.8, col.main="gray30", line=-0.2)
				dev.copy2pdf(file=paste0("CRU_",variables[h],"_NEW.pdf"))
			}
		if (h == 2)
			{
				dev.new(width=9, height=3); par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(1.0,1.0,0.5,0.5), mar=c(2.0,2.0,1,1))
				plot(curves1[[1]][,1], curves1[[1]][,2], lwd=0.5, type="l", axes=F, ann=F)
				lines(curves2[[1]][,1], curves2[[1]][,2], lwd=1.0, col="red")
				axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
				axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(135,240,10))
				title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
				title(ylab=yLabs[h], cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
				title(main="M-R-U", cex.main=0.8, col.main="gray30", line=-0.2)
				legend(1920, 140, legend=c("sliding window of 5 years"), lwd=1, cex=0.6, col=c("red"), text.col="gray30", border=NA, x.intersp=0.5, bty="n")
				plot(curves1[[3]][,1], curves1[[3]][,2], lwd=0.5, type="l", axes=F, ann=F)
				lines(curves2[[3]][,1], curves2[[3]][,2], lwd=1.0, col="red")
				axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
				axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(90,150,10))
				title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
				title(ylab=yLabs[h], cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
				title(main="Nigeria", cex.main=0.8, col.main="gray30", line=-0.2)
				dev.copy2pdf(file=paste0("CRU_",variables[h],"_NEW.pdf"))
			}
	}

	# 8.2. Preparation of land cover covariates
	
		# Preparation of past land use estimates based on data from the Harmonized Global Land Use for Years 1500-2100 (V1)
		# Source: https://daac.ornl.gov/VEGETATION/guides/Land_Use_Harmonization_V1.html
			# "We obtained scenarios of 21st century land-use from the Land-Use Harmonization project (i.e. LUHa_u2t1.v1; http://luh.umd.edu; 
			# Hurtt et al. (2011) Harmonization of land-use scenarios for the period 1500–2100: 600 years of global gridded annual land-use 
			# transitions, wood harvest, and resulting secondary land. Climatic Change 109: 117–161)"

selected_years = c(1901:2005); curves_list = list()
xMin = -16; xMax = 16; yMin = 3; yMax = 15
variables = c("gcrop","gothr","gsecd","gpast","gurbn")
suffixes = c("Crop","Prim","Secd","Past","Urbn")
variable_names = c("cropland","primary land","secondary land","pasture","urban land")
for (h in 1:length(variables))
	{
		nc = nc_open(paste0("Land_use_for_1500-2100/LUHa_u2t1.v1_",variables[h],".nc4"))
		time = ncvar_get(nc,"time")
		lon = ncvar_get(nc,"lon")
		lat = ncvar_get(nc,"lat")
		if (h == 1) var1 = ncvar_get(nc,"prop_crop")
		if (h == 2) var1 = ncvar_get(nc,"prop_primary")
		if (h == 3) var1 = ncvar_get(nc,"prop_secd")
		if (h == 4) var1 = ncvar_get(nc,"prop_pasture")
		if (h == 5) var1 = ncvar_get(nc,"prop_urbn")
		indices1 = which((lon>xMin)&(lon<xMax))
		indices2 = which((lat>yMin)&(lat<yMax))
		lon = lon[indices1]; lat = lat[indices2]
		var1 = var1[indices1,indices2,]
		years = time+1700; n = 0; var2 = c()
		tot = dim(var1)[1]*dim(var1)[2]
		for (i in 1:dim(var1)[1])
			{
				for (j in 1:dim(var1)[2])
					{
						n = n+1
						row = cbind(lon[i],lat[j],t(var1[i,j,]))
						var2 = rbind(var2, row)
					}
			}
		colnames(var2) = c("longitude","latitude",years); curves = list()
		analyses = c("LASV2_MRU_segL","LASV2_MRU_segS","LASV2_NIG_segL","LASV2_NIG_segS")
		for (i in 1:length(analyses))
			{
				values = matrix(nrow=length(selected_years), ncol=2)
				for (y in 1:length(selected_years))
					{
						buffer = var2[,3:dim(var2)[2]]; colNames = gsub("X","",colnames(buffer))
						selectedMeasures = cbind(var2[,1:2],buffer[,which(floor(as.numeric(colNames))==selected_years[y])])
						hull = readOGR(dsn="Minimum_convex_hulls", layer=paste0(analyses[i],"_",selected_years[y]))
						inHull = point.in.polygon(selectedMeasures[,"longitude"],selectedMeasures[,"latitude"],
								 		   		  hull@polygons[[1]]@Polygons[[1]]@coords[,1],
								 		   		  hull@polygons[[1]]@Polygons[[1]]@coords[,2])
						indices = which(inHull==1); values[y,1] = selected_years[y]
						values[y,2] = mean(as.matrix(selectedMeasures[indices,3:dim(selectedMeasures)[2]]), na.rm=T)
					}
				curves[[i]] = values; colnames(values) = c("time",variables[h])
				write.csv(values, paste0(gsub("LASV2",suffixes[h],analyses[i]),".csv"), row.names=F)
			}
		curves_list[[h]] = curves
	}
analyses = c("LASV2_MRU_segL","LASV2_MRU_segS","LASV2_NIG_segL","LASV2_NIG_segS")
for (h in 1:length(variables))
	{
		curves = list()
		for (i in 1:length(analyses))
			{
				curves[[i]] = read.csv(paste0(gsub("LASV2",suffixes[h],analyses[i]),".csv"), header=T)
			}
		curves_list[[h]] = curves
	}
entropies = list()
for (i in 1:length(analyses))
	{
		proportions = matrix(nrow=dim(curves_list[[1]][[i]])[1], ncol=length(variables)+1)
		proportions[,1] = curves_list[[1]][[i]][,1]
		for (h in 1:length(variables))
			{
				proportions[,h+1] = curves_list[[h]][[i]][,2]
			}
		entropy = matrix(nrow=dim(proportions)[1], ncol=2); entropy[,1] = proportions[,1]
		for (j in 1:dim(proportions)[1])
			{
				p = 0
				for (k in 2:dim(proportions)[2])
					{
						if (proportions[j,k] != 0)
							{
								p_i = proportions[j,k]
								p = p - ((p_i)*(log(p_i)))
							}
					}
				entropy[j,2] = p
			}
		entropies[[i]] = entropy; colnames(values) = c("time","entropy")
		write.csv(values, paste0(gsub("LASV2","Enty",analyses[i]),".csv"), row.names=F)
	}
cols = c("olivedrab3","mediumseagreen","steelblue3","orange1","tan3")
cols = c("orange1","olivedrab3","mediumseagreen","gray60","steelblue3")
dev.new(width=9, height=3); par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(1.0,1.0,0.5,0.5), mar=c(2.0,2.0,1,1))
plot(curves_list[[1]][[1]][,1], curves_list[[1]][[1]][,2], col=cols[1], lwd=1.5, type="l", axes=F, ann=F, ylim=c(0,0.9))
for (i in 2:5) lines(curves_list[[i]][[1]][,1], curves_list[[i]][[1]][,2], col=cols[i], lwd=1.5)
axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(0,0.9,0.1))
title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
title(ylab="land use proportion", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
title(main="M-R-U", cex.main=0.8, col.main="gray30", line=-0.2)
plot(curves_list[[1]][[3]][,1], curves_list[[1]][[3]][,2], col=cols[1], lwd=1.5, type="l", axes=F, ann=F, ylim=c(0,0.9))
for (i in 2:5) lines(curves_list[[i]][[3]][,1], curves_list[[i]][[3]][,2], col=cols[i], lwd=1.5)
axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(0,0.9,0.1))
title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
title(ylab="land use proportion", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
title(main="Nigeria", cex.main=0.8, col.main="gray30", line=-0.2)
legend(1975, 0.97, legend=variable_names, lwd=2, cex=0.6, col=cols, text.col="gray30", border=NA, x.intersp=0.5, bty="n")
dev.copy2pdf(file=paste0("Land_use_1901-2005.pdf"))
dev.new(width=9, height=3); par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(1.0,1.0,0.5,0.5), mar=c(2.0,2.0,1,1))
plot(entropies[[1]][,1], entropies[[1]][,2], col="red", lwd=1.5, type="l", axes=F, ann=F)
axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(0.75,1.00,0.05))
title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
title(ylab="land use entropy", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
title(main="M-R-U", cex.main=0.8, col.main="gray30", line=-0.2)
plot(entropies[[3]][,1], entropies[[3]][,2], col="red", lwd=1.5, type="l", axes=F, ann=F)
axis(side=1, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,-0.05,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1900,2020,10))
axis(side=2, lwd.tick=0.2, cex.axis=0.5, mgp=c(0,0.25,0), lwd=0.2, tck=-0.020, col.axis="gray30", at=seq(1.03,1.17,0.01))
title(xlab="time (year)", cex.lab=0.6, mgp=c(0.7,0,0), col.lab="gray30")
title(ylab="land use entropy", cex.lab=0.6, mgp=c(1.1,0,0), col.lab="gray30")
title(main="Nigeria", cex.main=0.8, col.main="gray30", line=-0.2)
legend(1975, 0.97, legend=variable_names, lwd=2, cex=0.6, col=cols, text.col="gray30", border=NA, x.intersp=0.5, bty="n")
dev.copy2pdf(file=paste0("Entropy_1901-2005.pdf"))

	# 8.3. Preparation of the BEAST skygrid-GLM analyses

		# Models used: GTR+G substitution model, relaxed lognormal molecular clock model, skygrid coalescent model (grid points defined to have one point until 1901)

analyses = c("MRU_segL","MRU_segS","NIG_segL","NIG_segS")
mostRecentSamplingDates = c(2018.1671232876713, 2018.1671232876713, 2019.1780821917807, 2019.1780821917807)
for (a in 1:length(analyses))
	{
		mostRecentSamplingDatum = mostRecentSamplingDates[a]
		covariates = list(); dates = list(); ylims = list()
		covariates[[1]] = read.csv(paste0("Temp_",analyses[a],".csv"), header=T)
		covariates[[2]] = read.csv(paste0("Prec_",analyses[a],".csv"), header=T)
		# covariates[[3]] = read.csv(paste0("Enty_",analyses[a],".csv"), header=T)
		covariates[[3]] = read.csv(paste0("Crop_",analyses[a],".csv"), header=T)
		covariates[[4]] = read.csv(paste0("Prim_",analyses[a],".csv"), header=T)
		covariates[[5]] = read.csv(paste0("Secd_",analyses[a],".csv"), header=T)
		covariates[[6]] = read.csv(paste0("Past_",analyses[a],".csv"), header=T)
		covariate_names = c(paste0("Temp_",analyses[a]), paste0("Prec_",analyses[a]), # paste0("Enty_",analyses[a]),
							paste0("Crop_",analyses[a]), paste0("Prim_",analyses[a]),
							paste0("Secd_",analyses[a]), paste0("Past_",analyses[a]))
		gridPoints = round(mostRecentSamplingDatum-rev(1902:floor(mostRecentSamplingDatum)),4)
		gridPoints_txt = gridPoints[1]
		for (i in 2:length(gridPoints))
			{
				gridPoints_txt = paste(gridPoints_txt,gridPoints[i],sep=" ")
			}
		timePoints = matrix(nrow=length(gridPoints)+1,ncol=1+length(covariates))
		colnames(timePoints) = c("time",covariate_names)
		timePoints[,1] = c(1901,rev(mostRecentSamplingDatum-gridPoints))
		for (i in 1:dim(timePoints)[1])
			{
				timePoint = as.integer(timePoints[i,1])
				for (j in 1:length(covariates))
					{
						if (sum(timePoint==as.numeric(covariates[[j]][,1])) == 1)
							{
								timePoints[i,covariate_names[j]] = covariates[[j]][which(covariates[[j]][,1]==timePoint),2]
							}
					}
			}
		nberOfIntervals = dim(timePoints)[1]
		txt = c()
		txt = c(txt, paste0("\t\t<populationSizes>"))
		txt = c(txt, paste0("\t\t\t<parameter id=\"skygrid.logPopSize\" dimension=\"",nberOfIntervals,"\"  value=\"1.0\"/>"))	
		txt = c(txt, paste0("\t\t</populationSizes>"))
		txt = c(txt, "")
		txt = c(txt, paste0("\t\t<lastObservedIndex>"))
		txt = c(txt, paste0("\t\t\t<parameter id=\"skygrid.lastObservedIndex.",1,"\" value=\"",length(gridPoints)+1,"\"/>"))
		if (length(covariates) > 1)
			{
				for (i in 2:length(covariates))
					{
						txt = c(txt, paste0("\t\t\t<parameter id=\"skygrid.lastObservedIndex.",i,"\" value=\"",length(gridPoints)+1,"\"/>"))
					}
			}
		txt = c(txt, paste0("\t\t</lastObservedIndex>"))
		txt = c(txt, paste0("\t\t<covariatePrecision>"))
		for (i in 1:length(covariates))
			{
				maxDiff = 0
				for (j in 2:dim(timePoints)[1])
					{
						if ((!is.na(as.numeric(timePoints[j,covariate_names[i]])))&(!is.na(as.numeric(timePoints[j-1,covariate_names[i]]))))
							{
								if (maxDiff < abs(as.numeric(timePoints[j,covariate_names[i]])-as.numeric(timePoints[j-1,covariate_names[i]])))
									{
										maxDiff = abs(as.numeric(timePoints[j,covariate_names[i]])-as.numeric(timePoints[j-1,covariate_names[i]]))
									}
							}
					}
				values = as.numeric(timePoints[,covariate_names[i]]); values = values[!is.na(values)]
				covPrec = 1/((sd(values))^2)
				txt = c(txt, paste0("\t\t\t<parameter id=\"skygrid.covPrec.",i,"\" value=\"",covPrec,"\"/>"))
			}
		txt = c(txt, paste0("\t\t</covariatePrecision>"))
		txt = c(txt, "")
		txt = c(txt, paste0("\t\t<gridPoints>"))
		txt = c(txt, paste0("\t\t\t<parameter id=\"skygrid.gridPoints\" value=\"",gridPoints_txt,"\"/>"))
		txt = c(txt, paste0("\t\t</gridPoints>"))
		txt = c(txt, paste0("\t\t<precisionParameter>"))
		txt = c(txt, paste0("\t\t\t<parameter id=\"skygrid.precision\" value=\"0.1\" lower=\"0.0\"/>"))
		txt = c(txt, paste0("\t\t</precisionParameter>"))
		txt = c(txt, "")
		txt = c(txt, paste0("\t\t<betaParameter>"))
		for (i in 1:length(covariates))
			{
				txt = c(txt, paste0("\t\t\t<parameter id=\"skygrid.beta.",i,"\" value=\"0\" lower=\"-1000\" upper=\"1000\"/>"))
			}
		txt = c(txt, paste0("\t\t</betaParameter>"))
		txt = c(txt, "")
		txt = c(txt, paste0("\t\t<covariates>"))
		for (i in 1:length(covariates))
			{
				c = 0; lastValue = "NA"; meanValue = mean(as.numeric(timePoints[,covariate_names[i]]), na.rm=T)
				txt = c(txt, paste0("\t\t\t<matrixParameter id=\"",covariate_names[i],"\">"))
				for (j in dim(timePoints)[1]:1)
					{
						c = c+1
						if (!is.na(timePoints[j,covariate_names[i]]))
							{
								txt = c(txt, paste0("\t\t\t\t<parameter id=\"cov",i,".",c,"\" value=\"",timePoints[j,covariate_names[i]],"\"/>"))
								lastValue = timePoints[j,covariate_names[i]]
							}	else	{
								startValue = timePoints[max(which(!is.na(timePoints[,covariate_names[i]]))),covariate_names[i]]
								txt = c(txt, paste0("\t\t\t\t<parameter id=\"cov",i,".",c,"\" value=\"",startValue,"\" lower=\"0\"/>"))
								# txt = c(txt, paste0("\t\t\t\t<parameter id=\"cov",i,".",c,"\" value=NA","\" lower=\"0\"/>"))
							}
					}
				txt = c(txt,"\t\t\t</matrixParameter>")
			}
		txt = c(txt,"\t\t</covariates>")
		txt = c(txt,"","XXXXXXXXXXXXXXXXXX","")
		for (i in 1:length(covariates))
			{
				txt = c(txt, paste0("\t\t<randomWalkOperator windowSize=\"1.0\" weight=\"30\">"))
				txt = c(txt, paste0("\t\t\t<parameter idref=\"skygrid.beta.",i,"\"/>"))
				txt = c(txt, paste0("\t\t</randomWalkOperator>"))
			}
		txt = c(txt,"","XXXXXXXXXXXXXXXXXX","")	
		for (i in 1:length(covariates))
			{
				txt = c(txt, paste0("\t\t\t<normalPrior mean=\"0.0\" stdev=\"1\" offset=\"0.0\">"))
				txt = c(txt, paste0("\t\t\t\t<parameter idref=\"skygrid.beta.",i,"\"/>"))
				txt = c(txt, paste0("\t\t\t</normalPrior>"))
			}
		txt = c(txt,"","XXXXXXXXXXXXXXXXXX","")
		for (i in 1:length(covariates)) txt = c(txt, paste0("\t\t\t<parameter idref=\"skygrid.beta.",i,"\"/>"))
		txt = c(txt, paste0("\t\t\t<parameter idref=\"skygrid.gridPoints\"/>"))
		write(txt, paste0("GLM_",analyses[a],".xml"))
	}

# 9. Analysis of skygrid-GLM analyses

XXXX

# 10. Species distribution modelling

setwd(paste(wd,wd10,sep="/"))

	# 10.1. Preparation of environmental rasters

envVariableNames = read.csv("Environmental_rasters/LC_vars.csv")[,2]
# Source: http://luh.umd.edu/data.shtml; see also Lawrence et al. (2016, Geosci. Model Dev.)
fileName = "Environmental_rasters/Historical/GFDL-ESM2M/landcover_historical_annual_1986_2005_timmean.nc4"
land_cover = nc_open(fileName); varNames = names(land_cover$var); land_covers = list(); cols = list()
for (i in 2:13) land_covers[[i-1]] = brick(fileName, varname=varNames[i])
for (i in 1:length(land_covers))
	{
		names(land_covers[[i]]) = varNames[i+1]
		cols[[i]] = colorRampPalette(brewer.pal(9,"YlGn"))(120)[11:110]
		if (i == 1)
			{
				r_global = land_covers[[1]]
			}	else		{
				r_global[] = r_global[]+land_covers[[i]][]	
			}
	}
if (!file.exists("Land_cover_rasters.pdf"))	
	{
		pdf("Land_cover_rasters.pdf", width=7.5, height=5.0); par(mfrow=c(4,3), oma=c(1.5,2.0,1,0.5), mar=c(0,0,0,0), mgp=c(1,0.2,0), lwd=0.2)
		for (i in 1:12)
			{
				plot(land_covers[[i]], bty="n", box=F, axes=F, legend=F, col=c("gray90",cols[[i]]), colNA="white")
				plot(land_covers[[i]], legend.only=T, add=T, col=cols[[i]], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.09,0.105,0.14,0.49), adj=3,
					 axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1, side=3)
				mtext(envVariableNames[i], side=1, adj=0.5, line=-1.8, at=40, cex=0.55, font=1, col="gray30")
			}
		dev.off()
	}
mask = raster("Environmental_rasters/Mask.nc4")
nc_open("Environmental_rasters/Historical/GFDL-ESM2M/population_histsoc_0p5deg_annual_1986_2005_timmean.nc4")
nc_open("Environmental_rasters/Historical/GFDL-ESM2M/tas_day_GFDL-ESM2M_historical_r1i1p1_EWEMBI_1986_2005_timmean.nc4")
nc_open("Environmental_rasters/Historical/GFDL-ESM2M/pr_day_GFDL-ESM2M_historical_r1i1p1_EWEMBI_1986_2005_timmean.nc4")
population = raster("Environmental_rasters/Historical/GFDL-ESM2M/population_histsoc_0p5deg_annual_1986_2005_timmean.nc4")
temperature = raster("Environmental_rasters/Historical/GFDL-ESM2M/tas_day_GFDL-ESM2M_historical_r1i1p1_EWEMBI_1986_2005_timmean.nc4")
precipitation = raster("Environmental_rasters/Historical/GFDL-ESM2M/pr_day_GFDL-ESM2M_historical_r1i1p1_EWEMBI_1986_2005_timmean.nc4")
population[is.na(mask[])] = NA; # unit: number of people # plot(population)
temperature[is.na(mask[])] = NA; # unit: K (Kelvin, °C+273.15) # plot(temperature)
precipitation[is.na(mask[])] = NA; # unit: kg/m2/second # plot(precipitation)

