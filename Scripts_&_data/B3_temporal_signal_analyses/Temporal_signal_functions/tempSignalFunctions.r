############################################################################################################
### Script with various functions for tests of temporal signal                                                                                      							  ###
### (c) 2015 Gemma GR Murray, John J Welch                                                                      							                                      ###
############################################################################################################

############################################################################################################
## Return a tree rooted at a given node
prepare.phy <- function(phy,node,nt)
{
if(node<=nt){
	#root with external branch
	rooted_tree<-root(phy,node,resolve.root=T)	
	b<-sum(rooted_tree$edge.length[rooted_tree$edge==(nt+1)])
	rooted_tree$edge.length[rooted_tree$edge==(nt+1)]<-c(b/2,b/2)
} else {
	#root with internal branch
	rooted_tree<-root(phy,node=node,resolve.root=T)	
	b<-sum(rooted_tree$edge.length[rooted_tree$edge==(nt+1)])
	rooted_tree$edge.length[rooted_tree$edge==(nt+1)]<-c(b/2,b/2)
}
return(rooted_tree)
}
############################################################################################################
## Return all possible rootings of a phylogeny
get.all.roots <- function(phy)
{
	nt <- length(phy$tip.label) 
	n <- nt+phy$Nnode
	phy.all.roots <- vector('list',n)
	for(i in c(1:nt,(nt+2):n))
	{
		phy.all.roots[[i]] <- prepare.phy(phy,i,nt)
	}
	phy.all.roots[-(nt+1)]	
}	
############################################################################################################
get.altered.tree <- function(phy,p,root.node)
{
	w <- which(phy$edge[,1]==root.node)
	s <- phy$edge.length[w[1]] + phy$edge.length[w[2]]
	phy$edge.length[w[1]] <- p*s
	phy$edge.length[w[2]] <- (1-p)*s
	phy
}
############################################################################################################
## Return a summary stat of a given regression of root-to-tip distance on tip sampling date.
myfun <- function(p,phy,dates,stat='nve-rms',nt)
{
	phy <- get.altered.tree(phy,p,nt+1)
	div <- dist.nodes(phy)[nt+1,1:nt] # Root-to-tip distances
	
	mylm <- lm(div~dates)
	if(stat=='r') {
		s <- summary(mylm)
		return(sqrt(s$r.squared)*sign(s$coefficients[2,1]))
	} else if(stat=='r2')
		return(summary(mylm)$r.squared)
	else if(stat=='nve-rms') {
		#return(-1e5*var(residuals(mylm))*(nt-1)/(nt-2))
		return(-1e5*(sum((mean(residuals(mylm))-residuals(mylm))^2)/mylm$df.residual))
	}
	else
		stop('unrecognised stat')	
}
############################################################################################################
optfun <- function(phy,dates,nt,stat='nve-rms')
{
	ans <- optim(runif(1),fn=myfun,method='L-BFGS',lower=0,upper=1,control=list(fnscale=-1),phy=phy,dates=dates,stat=stat,nt=nt)
	return(c(ans$par,ans$value))
}
############################################################################################################
## Find the "best" root for a tree
get.best.root <- function(phy.all.roots,dates,stat='nve-rms')
{
	nt <- length(phy.all.roots[[1]]$tip.label)
	p <- lapply(phy.all.roots,FUN=optfun,dates=dates,stat=stat,nt=nt)	
	p <- matrix(unlist(p),length(phy.all.roots),2,byrow=T)
	w   <- which.max(p[,2])
	phy <- get.altered.tree(phy.all.roots[[w]],p[w,1],nt+1)
	phy
}
############################################################################################################
## Return nreps *unique* random permutations (each column is a permutation). 
## If luc is not NULL, it permutes over clusters of sequences, each sharing the same date
perm.rand <- function(nreps, mm, luc=NULL, clusters=NULL) { # Obtain nreps unique permutations of mm
    # Function to obtain a new permutation.
    newperm <- function() {
        count <- 0                # Protects against infinite loops
        repeat {
            # Generate a permutation and check against previous ones.
    		if(is.null(luc)) {
	            s <- sample(c(1:length(mm[1,])))
	            p <-mm[2,s]
    		} else {	
    			s <- sample.int(luc)[clusters]
    			p <- c()
    			for(i in 1:length(mm[1,])){  			
    				p[i]<-sample(as.character(mm[2,which(mm[1,]==s[i])]),1)
            }}
            
            hash.s <- paste(s, collapse="")
            if (is.null(cache[[hash.s]])) break

            # Prepare to try again.
            count <- count+1
            if (count > 1000) {   # 1000 is arbitrary; adjust to taste
                p <- NA           # NA indicates a new permutation wasn't found
                hash.s <- ""
                break
            }
        }
        cache[[hash.s]] <<- TRUE  # Update the list of permutations found
        p                         # Return this (new) permutation
    }
    # Obtain m unique permutations.
    cache <- list()
    cache[[ paste(clusters, collapse="")]] <- TRUE
    replicate(nreps, newperm())  
}
############################################################################################################
### Returns a grouping of strains in clades
get.clades <- function(phy)
{
	phy <- reorder(phy)
	ue <- unique(phy$edge[,1])	
	fm <- vector("list",phy$Nnode)
	root.lab <- length(phy$tip.label)+1 
	for(i in 1:length(phy$tip.label))
	{
		z <- i
		while(z != root.lab) {
			new.z <- phy$edge[which(phy$edge[,2]==z),1]
			m <- match(new.z,ue)
			fm[[m]] <- c(fm[[m]],phy$tip.label[i])
			z <- new.z
		}
	}
	return(fm)
}
############################################################################################################
## Return the number of unique sampling dates associated with each clade
checkndates <- function(clade,dates,noms)
{
	length(unique(dates[match(clade,noms)]))
}
############################################################################################################
## Plot the clusters
plot.clusters <- function(phy,cl)
{
	cols <-rainbow(max(cl), start=0, end=(1-1/(max(cl))), v=rep(c(1, 0.9, 0.8),length(cl)),s=rep(c(1,0.7),length(cl)))[sample(max(cl))]
	
	phy$tip.label <- rep('.',length(phy$tip.label))
	plot(phy,tip.color=cols[cl],cex=5,root.edge=T)
	title(main=paste('No clusters =',max(cl)))
}
############################################################################################################
## Get clusters of sequences that share the same date, and produce random permutations of these clusters
get.clusters <- function(phy, dates, auto.round.dates=F, rounded.dates=NULL, do.plot=F,stat='nve-rms', reroot=T)
{
	if(!is.rooted(phy)| reroot==T)
	{	
		phy.allroots <- get.all.roots(phy)	
		phy <- get.best.root(phy.allroots,dates,stat)
	}

	if(auto.round.dates==T) 
	{
		dates<-round(dates+0.0000000000000000000000001)
	}
	
	cl <- 1:length(dates)
		
	gc <- get.clades(phy)
	if(!is.null(rounded.dates)) da <- rounded.dates
	if(is.null(rounded.dates)) da <- dates
	
	single.date.clades <- which(unlist(lapply(gc,FUN=checkndates,dates=da,noms=phy$tip.label))==1)
	if(length(single.date.clades)>0)
	{
		for(i in 1:length(single.date.clades))
		{
			m <- match(gc[[single.date.clades[i]]],phy$tip.label)		
			cl[m] <- cl[m[1]]
		}
	}
	cl <- match(cl,unique(cl))
	if(length(unique(cl))==1) stop("Only one cluster in data. Dates need to be finer grained.")
	if(do.plot)
		plot.clusters(phy,cl)
	return(cl)
}
############################################################################################################