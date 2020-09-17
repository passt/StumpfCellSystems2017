# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Neuronal time course - Fluidigm arrays - Load and Normalise Data
# 2016-02-09
# Patrick Stumpf
# University of Southampton
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# load data

t0 <- read.table("./data/2015-01-08 0H sample Ct values.txt", header=T, sep="\t")
rownames(t0) <- t0[,1]
t0 <- t0[,c(-1,-98)]

t24 <- read.table("./data/2015-01-08 24H sample Ct values.txt", header=T, sep="\t")
rownames(t24) <- t24[,1]
t24 <- t24[,-1]

t48 <- read.table("./data/2015-01-08 48H sample Ct values.txt", header=T, sep="\t")
rownames(t48) <- t48[,1]
t48 <- t48[,-1]

t72 <- read.table("./data/2015-01-08 72H sample Ct values.txt", header=T, sep="\t")
rownames(t72) <- t72[,1]
t72 <- t72[,-1]

t96 <- read.table("./data/2015-01-08 96H sample Ct values.txt", header=T, sep="\t")
rownames(t96) <- t96[,1]
t96 <- t96[-97,c(-1,-98)]

t120 <- read.table("./data/2015-01-08 120H sample Ct values.txt", header=T, sep="\t")
rownames(t120) <- t120[,1]
t120 <- t120[,-1]

t168 <- read.table("./data/2015-01-08 168H sample Ct values.txt", header=T, sep="\t")
rownames(t168) <- t168[,1]
t168 <- t168[-97,c(-1,-98)]

# combine individual timepoints
timecourse <- cbind(t0, t24, t48, t72, t96, t120, t168)

# remove redundant variables
rm(t0, t120, t168, t24, t48, t72, t96)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Sample annotation and colour
	sample.factor <- as.factor(matrix(rep(c(0,24,48,72,96,120,168),96), ncol=7, byrow=T))
	line.factor   <- factor(rep("E14", length(sample.factor)), levels=c("E14", "R1"))
	line.factor[grep("R1.", colnames(timecourse))] <- "R1"				   

# Gene Annotation
	gene.info <- read.table("./data/GeneClass.txt", sep="\t", row.names=1, header=F,
							col.names=c("Gene","Category", "TaqManProbeID"))							
	# set colours for gene categories
	gene.info[,3] <- factor(gene.info[,1], labels=brewer.pal(12, "Paired")[c(5,2,3,9,11,7,4,6,8,10,1,12)])
	colnames(gene.info) <- c(colnames(gene.info[1:2]),"Colour")
      
# colourschemes
	# samples:
	col.time <- brewer.pal(9, "YlGnBu")[2:8]	
	# lines:	
	col.line <- c("midnightblue", "darkmagenta")


# Heatmap
	heatmap.col <- maPalette(low="#FFFFFF", high="#000000", k=24)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Quality control (1) for dataset:
# Remove cells with missing readings for housekeeping genes
	no.actb 		<- which(timecourse[1,] >15)
	no.gapdh 	<- which(timecourse[26,]>15)
# union == 27 samples
	no.hk <- union(no.actb, no.gapdh)
# samples without housekeeping genes are missing 86% of readings
	sum(timecourse[, no.hk]==999) / prod(dim(timecourse[, no.hk]))
# other samples are missing only 47% of readings (a lot)
	sum(timecourse[, -no.hk]==999) / prod(dim(timecourse[, -no.hk]))

# remove samples with missing HK expression from datset and other factors
	timecourse 		<- timecourse[, -no.hk]
	sample.factor 	<- sample.factor[-no.hk]
	line.factor		<- line.factor[-no.hk]

# visualise removal
	barplot(rbind(table(sample.factor)/96, 1-table(sample.factor)/96),
			ylim=c(0,1), col=c("grey40", "firebrick"), border=F, las=1,
			ylab="Fraction", main="QC (1) - Fraction of cells removed")
  legend('bottomright', legend=c('removed', 'retained'), fill=c("firebrick","grey40"))
	
# remove redundant variables
	rm(no.actb, no.gapdh, no.hk)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# normalisation per array (equivalent to previous PDC arrays)
# expression cut-off: Ct 28
# create new data.frame for normalisation
	timecourse.log2 <- timecourse
# set values na beyond expression cut-off
	is.na(timecourse.log2) <- timecourse > 28
# calculate median per array for mean of hk genes Actb: #1 and Gapdh: #26
	median.hk <-  tapply( t( colSums(timecourse.log2[c(1,26),])/2 ),   sample.factor, median, na.rm=T)
# create vector with median.hk values of same length as columns in data
	med.fac <- sample.factor
	levels(med.fac) <- median.hk
	med.vec <- as.numeric(as.vector(med.fac))
# substract Ct from each median
	timecourse.log2 <- t(apply(timecourse.log2, 1, function(x)(med.vec-x)))
# find minimum and substract from each datum to create baseline at 0
	timecourse.log2 <- timecourse.log2 - min(timecourse.log2, na.rm=T)
# set na values to new baseline
	timecourse.log2[is.na(timecourse.log2)] <- 0

# remove redundant variables
rm(timecourse, med.fac, med.vec, median.hk)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Quality control (2) for dataset:
# Remove cells with cumulative expression (E)  (E < Q1 - 2 * IQR) | (E > Q3 + 2 * IQR)
E <- colSums(timecourse.log2)

# Summary Stats
E.IQR <- tapply(E, sample.factor, IQR)
E.quantiles <- matrix(unlist(tapply(E, sample.factor, quantile, probs=c(.25,.75))), ncol=2, byrow=T)
rownames(E.quantiles) <- levels(sample.factor)
colnames(E.quantiles) <- c("Q1","Q3")

# lower and upper limits (factor)
E.ll <- E.ul <- sample.factor
levels(E.ll) <- E.quantiles[,1] - 2*E.IQR
levels(E.ul) <- E.quantiles[,2] + 2*E.IQR

# Identify outliers with (E < Q1 - 2 * IQR)
E.low  <- E < as.numeric(as.vector(E.ll))
# Identify outliers with (E > Q3 + 2 * IQR)
E.high <- E > as.numeric(as.vector(E.ul))

# subset dataset and factors
timecourse.log2 <- timecourse.log2[ , !(E.low | E.high)]
sample.factor 	<- sample.factor[!(E.low | E.high)]
line.factor 	<- line.factor[!(E.low | E.high)]

# set sample colours
sample.colours <- factor(sample.factor, labels=col.time)

# visualise removal
	barplot(rbind(table(sample.factor)/96, 1-table(sample.factor)/96),
			ylim=c(0,1), col=c("grey40", "firebrick"), border=F, las=1,
			ylab="Fraction", main="QC (2) - Fraction of cells removed")
	legend('bottomright', legend=c('removed', 'retained'), fill=c("firebrick","grey40"))
	
	
cat("Cumulative number of cells excluded in QC (1 & 2): \n\n")
print(48-matrix(unlist(tapply(sample.factor, line.factor, table)), nrow=2, byrow=T, 
                dimnames=list(Line=levels(line.factor),Time=levels(sample.factor))))

# remove redundant variables
rm(list=ls(pattern="E"))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# to make code downstream more legible
E14 <- grepl("E14", line.factor)
R1  <- grepl("R1",  line.factor)