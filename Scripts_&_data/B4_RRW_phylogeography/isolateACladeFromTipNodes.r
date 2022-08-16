isolateACladeFromTipNodes = function(tab, selectedTipNodes)
	{
		nodesInClade = matrix(nrow=dim(tab)[1], ncol=2)
		nodesInClade[which(tab[,"node2"]%in%selectedTipNodes),2] = 1
		tipNodes = tab[which(!tab[,"node2"]%in%tab[,"node1"]),"node2"]
		nodesInClade[which((tab[,"node2"]%in%tipNodes)&(!tab[,"node2"]%in%selectedTipNodes)),2] = 0
		branchesToScan1 = which(tab[,"node2"]%in%tipNodes)
		while (sum(is.na(nodesInClade)) != 0)
			{
				for (i in 1:dim(tab)[1])
					{
						if ((is.na(nodesInClade[i,1]))&(!is.na(nodesInClade[i,2])))
							{
								endNodesIndices = which(tab[,"node1"]==tab[i,"node1"])
								startNodesIndex = which(tab[,"node2"]==tab[i,"node1"])
								values = nodesInClade[endNodesIndices,2]
								if (sum(is.na(values)) == 0)
									{
										if (sum(values) == 0)
											{
												nodesInClade[endNodesIndices,1] = 0; nodesInClade[startNodesIndex,2] = 0
											}
										if (sum(values) == 1)
											{
												nodesInClade[endNodesIndices,1] = 0; nodesInClade[startNodesIndex,2] = 0
											}
										if (sum(values) == 2)
											{
												nodesInClade[endNodesIndices,1] = 1; nodesInClade[startNodesIndex,2] = 1
											}
									}
							}
					}
			}
		selectedBranches = nodesInClade[,1]+nodesInClade[,2]
		selectedBranches = which(selectedBranches == 2)
		return(tab[selectedBranches,])
	}
