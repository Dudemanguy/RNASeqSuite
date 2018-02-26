exactTest <- function(object, pair=1:2, dispersion="auto", rejection.region="doubletail", big.count=900, prior.count=0.125, adjust.method="BH", sort.by="FDR")
#	Calculates exact p-values for the differential expression levels of tags in the two groups being compared.
#	Davis McCarthy, Gordon Smyth.
#	Created September 2009. Last modified 8 July 2012.
{
#	Check input
	if(!is(object,"DataList")) stop("Currently only supports DataList objects as the object argument.")
	if(length(pair)!=2) stop("Pair must be of length 2.")
	rejection.region <- match.arg(rejection.region,c("doubletail","deviance","smallp"))

#	Get group names
	group <- as.factor(object$samples$group)
	levs.group <- levels(group)
	if(is.numeric(pair))
		pair <- levs.group[pair]
	else
		pair <- as.character(pair)	
	if(!all(pair %in% levs.group)) stop("At least one element of given pair is not a group.\n Groups are: ", paste(levs.group, collapse=" "))

#	Get dispersion vector
	if(is.null(dispersion)) dispersion <- "auto"
	if(is.character(dispersion)) {
		dispersion <- match.arg(dispersion,c("auto","common","trended","tagwise"))
		dispersion <- switch(dispersion,
			"common"=object$common.dispersion,
			"trended"=object$trended.dispersion,
			"tagwise"=object$tagwise.dispersion,
			"auto"=getDispersion(object)
		)
		if(is.null(dispersion)) stop("specified dispersion not found in object")
		if(is.na(dispersion[1])) stop("dispersion is NA")
	}
	ldisp <- length(dispersion)
	ntags <- nrow(object$counts)
	if(ldisp!=1 && ldisp!=ntags) stop("Dispersion provided by user must have length either 1 or the number of tags in the DataList object.")
	if(ldisp==1) dispersion <- rep(dispersion,ntags)

#	Reduce to two groups
	group <- as.character(group)
	j <- group %in% pair
	y <- object$counts[,j,drop=FALSE]
	lib.size <- object$samples$lib.size[j]
	norm.factors <- object$samples$norm.factors[j]
	group <- group[j]
	if(is.null(rownames(y))) rownames(y) <- paste("tag",1:ntags,sep=".")

#	Normalized library sizes
	lib.size <- lib.size * norm.factors
	offset <- log(lib.size)
	lib.size.average <- exp(mean(offset))

#	logFC
	prior.count <- prior.count*lib.size/mean(lib.size)
	offset.aug <- log(lib.size+2*prior.count)
	j1 <- group==pair[1]
	n1 <- sum(j1)
	if(n1==0) stop("No libraries for",pair[1])
	y1 <- y[,j1,drop=FALSE]
	abundance1 <- mglmOneGroup(y1+matrix(prior.count[j1],ntags,n1,byrow=TRUE),offset=offset.aug[j1],dispersion=dispersion)
	j2 <- group==pair[2]
	n2 <- sum(j2)
	if(n1==0) stop("No libraries for",pair[2])
	y2 <- y[,j2,drop=FALSE]
	abundance2 <- mglmOneGroup(y2+matrix(prior.count[j2],ntags,n2,byrow=TRUE),offset=offset.aug[j2],dispersion=dispersion)
	logFC <- (abundance2-abundance1)/log(2)

#	Equalize library sizes
	abundance <- mglmOneGroup(y,dispersion=dispersion,offset=offset)
	e <- exp(abundance)
	input.mean <- matrix(e,ntags,n1)
	output.mean <- input.mean*lib.size.average
	input.mean <- t(t(input.mean)*lib.size[j1])
	y1 <- q2qnbinom(y1,input.mean=input.mean,output.mean=output.mean,dispersion=dispersion)
	input.mean <- matrix(e,ntags,n2)
	output.mean <- input.mean*lib.size.average
	input.mean <- t(t(input.mean)*lib.size[j2])
	y2 <- q2qnbinom(y2,input.mean=input.mean,output.mean=output.mean,dispersion=dispersion)

	exact.pvals <- switch(rejection.region,
		doubletail=exactTestDoubleTail(y1,y2,dispersion=dispersion,big.count=big.count),
		deviance=exactTestByDeviance(y1,y2,dispersion=dispersion),
		smallp=exactTestBySmallP(y1,y2,dispersion=dispersion)
	)

#	grab raw counts and average by group
	width1 <- table(object$samples$group)[[1]]
	a <- object$counts[,1:width1]
	b <- object$counts[,(width1+1):ncol(object$counts)]
	c <- rowMeans(a)
	d <- rowMeans(b)

#	adjust.pvalues
	FWER.methods <- c("holm", "hochberg", "hommel", "bonferroni")
	FDR.methods <- c("BH", "BY", "fdr")
	adjust.method <- match.arg(adjust.method,c(FWER.methods,FDR.methods,"none"))
	adj.p.val <- p.adjust(exact.pvals, method=adjust.method)

#	add to DataList
	AveLogCPM <- object$AveLogCPM
	if(is.null(AveLogCPM)) AveLogCPM <- aveLogCPM(object)
	de.out <- data.frame('Avg Ct A'=c, 'Avg Ct B'=d, logFC=logFC, logCPM=AveLogCPM, PValue=exact.pvals, FDR=adj.p.val)
	rn <- rownames(object$counts)
	if(!is.null(rn)) rownames(de.out) <- make.unique(rn)
	o <- switch(sort.by,
		"logFC" = order(de.out$logFC, decreasing=TRUE),
		"logCPM" = order(de.out$logCPM, decreasing=TRUE),
		"PValue" = order(de.out$PValue, decreasing=FALSE),
		"FDR" = order(de.out$FDR, decreasing=FALSE),
		"none" = 1:nrow(de.out)
	)
	de.out <- de.out[o,]
	object[["comparison"]] <- pair
	object[["adjust.method"]] <- adjust.method
	object[["et_results"]] <- de.out
	object
}

exactTestBetaApprox <- function(y1,y2,dispersion=0)
#	Test for differences in means between two negative binomial
#	or Poisson random variables, or between two groups of variables,
#	using a beta distribution approximation.
#	Test is naturally conditional on total sum.
#	Left and right rejection regions have equal probability.

#	Gordon Smyth
#	28 Sep 2019.  Last modified 28 Sep 2011.
{
#	Convert matrices to vectors
	ntags <- NROW(y1)
	n1 <- NCOL(y1)
	n2 <- NCOL(y2)
	if(n1>1) y1 <- rowSums(y1)
	if(n2>1) y2 <- rowSums(y2)
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)

#	Null fitted values
	y <- y1+y2
	mu <- y/(n1+n2)

#	Compute p-values
	pvals <- rep(1,ntags)
	all.zero <- y<=0
	alpha1 <- n1*mu/(1+dispersion*mu)
	alpha2 <- n2/n1*alpha1
	med <- rep(0,ntags)
	med[!all.zero] <- qbeta(0.5,alpha1[!all.zero],alpha2[!all.zero])
	left <- (y1+0.5)/y<med & !all.zero
	if(any(left)) {
		pvals[left] <- 2*pbeta((y1[left]+0.5)/y[left],alpha1[left],alpha2[left])
	}
	right <- (y1-0.5)/y>med & !all.zero
	if(any(right)) {
		pvals[right] <- 2*pbeta((y1[right]-0.5)/y[right],alpha1[right],alpha2[right],lower.tail=FALSE)
	}
	names(pvals) <- names(y1)
	pvals
}

exactTestBySmallP <- function(y1,y2,dispersion=0) 
#	Test for differences in means between two groups of
#	negative binomial or Poisson random variables,
#	using exact enumeration conditional on total sum.

#	Rejection region is by method of small probability, i.e.,
#	all values with probability equal or less than that observed.

#	Mark Robinson, Davis McCarthy, Gordon Smyth.
#	Created 17 June 2009.  Last modified 9 Dec 2013.
{
	y1 <- as.matrix(y1)
	y2 <- as.matrix(y2)
	ntags <- nrow(y1)
	if(ntags!=nrow(y2)) stop("Number of rows of y1 not equal to number of rows of y2")
	if(any(is.na(y1)) || any(is.na(y2))) stop("NAs not allowed")
	n1 <- ncol(y1)
	n2 <- ncol(y2)

	if(n1==n2) return(exactTestDoubleTail(y1=y1,y2=y2,dispersion=dispersion))

	sum1 <- round(rowSums(y1))
	sum2 <- round(rowSums(y2))
	N <- sum1+sum2
	mu <- N/(n1+n2)
	if(all(dispersion==0)) return(binomTest(sum1,sum2,p=n1/(n1+n2)))
	if(any(dispersion==0)) stop("dispersion must be either all zero or all positive")
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)
	r <- 1/dispersion
	all.zeros <- N==0

	pvals <- rep(1,ntags)
	if(ntags==0) return(pvals)
	if(any(all.zeros)) {
		pvals[!all.zeros] <- Recall(y1=y1[!all.zeros,,drop=FALSE],y2=y2[!all.zeros,,drop=FALSE],dispersion=dispersion[!all.zeros])
		return(pvals)
	}
	for (i in 1:ntags) {
		ind <- 0:N[i]
		p.top <- dnbinom(ind,size=n1*r[i],mu=n1*mu[i])*dnbinom(N[i]-ind,size=n2*r[i],mu=n2*mu[i])
		p.obs <- dnbinom(sum1[i],size=n1*r[i],mu=n1*mu[i]) * dnbinom(sum2[i],size=n2*r[i],mu=n2*mu[i])
		keep <-  p.top<=p.obs
		p.bot <- dnbinom(N[i],size=(n1+n2)*r[i],mu=(n1+n2)*mu[i])
		pvals[i] <- sum(p.top[keep]/p.bot)
	}
	min(pvals,1)
}

exactTestDoubleTail <- function(y1,y2,dispersion=0,big.count=900)
#	Test for differences in means between two groups of
#	negative binomial or Poisson random variables,
#	using exact enumeration conditional on total sum.

#	Smaller tail probability is doubled to get p-value.
#	QUESTION: should we use sign(logFC) to choose which tail to evaluate
#	instead of trying to find smaller of tail probabilities?

#	Gordon Smyth
#	28 Sep 2019.  Last modified 10 Jan 2012.
{
#	Convert matrices to vectors
	ntags <- NROW(y1)
	n1 <- NCOL(y1)
	n2 <- NCOL(y2)
	if(n1>1) s1 <- round(rowSums(y1)) else s1 <- round(y1)
	if(n2>1) s2 <- round(rowSums(y2)) else s2 <- round(y2)
	if(length(dispersion)==1) dispersion <- rep(dispersion,ntags)

#	Null fitted values
	s <- s1+s2
	mu <- s/(n1+n2)
	mu1 <- n1*mu
	mu2 <- n2*mu

	pvals <- rep(1,ntags)
	names(pvals) <- names(y1)

#	Poisson case
	pois <- dispersion<=0
#	BINOMTEST DOESN'T USE EQUAL TAILED REJECTION REGION
	if(any(pois)) pvals[pois] <- binomTest(s1[pois],s2[pois],p=n1/(n1+n2))

#	Use beta approximation for large counts
	big <- s1>big.count & s2>big.count
	if(any(big)) {
		y1 <- as.matrix(y1)
		y2 <- as.matrix(y2)
		pvals[big] <- exactTestBetaApprox(y1[big,,drop=FALSE],y2[big,,drop=FALSE],dispersion[big])
	}

	p.bot <- size1 <- size2 <- rep(0,ntags)
	left <- s1<mu1 & !pois & !big
	if(any(left)) {
		p.bot[left] <- dnbinom(s[left],size=(n1+n2)/dispersion[left],mu=s[left])
		size1[left] <- n1/dispersion[left]
		size2[left] <- n2/dispersion[left]
		for (g in which(left)) {
			x <- 0:s1[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(s[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[left] <- pvals[left]/p.bot[left]
	}
	right <- s1>mu1 & !pois & !big
	if(any(right)) {
		p.bot[right] <- dnbinom(s[right],size=(n1+n2)/dispersion[right],mu=s[right])
		size1[right] <- n1/dispersion[right]
		size2[right] <- n2/dispersion[right]
		for (g in which(right)) {
			x <- s1[g]:s[g]
			p.top <- dnbinom(x,size=size1[g],mu=mu1[g]) * dnbinom(s[g]-x,size=size2[g],mu=mu2[g])
			pvals[g] <- 2*sum(p.top)
		}
		pvals[right] <- pvals[right]/p.bot[right]
	}
	pmin(pvals,1)
}
