library(lubridate)

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
dev.off(); # dev.copy2pdf(file=paste("Mismatch_distributions.pdf",sep=""))

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
