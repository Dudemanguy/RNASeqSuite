#  SUBSET DATA SETS

assign("[.DataList",
function(object, i, j, keep.lib.sizes=TRUE) {
#  Subsetting for DataList objects
#  24 September 2009.  Last modified 8 Feb 2015.
	if(nargs() < 3) {
		stop("Two subscripts required",call.=FALSE)
	}

#	Recognized components
	IJ <- c("counts","pseudo.counts","offset","weights")
	IX <- c("et_results","genes")
	JX <- c("samples")
	I  <- c("AveLogCPM","trended.dispersion","tagwise.dispersion","prior.n","prior.df")
#	Obsolete <- c("conc","infos","all.zeros")

	out <- subsetListOfArrays(object,i,j,IJ=IJ,IX=IX,I=I,JX=JX)
	if (!(missing(i) || keep.lib.sizes)) {
		out$samples$lib.size <- colSums(out$counts)
	}
	out
})

