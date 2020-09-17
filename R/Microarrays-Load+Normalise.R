# Download data from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5861/

## Load and Normalise Microarrays
  # load time course data
  timecourse.microarrays <- lumiR("./data/GS_Report.txt")
  # normalise microarrays using RSN method
  tcn <- lumiN(timecourse.microarrays, method='rsn')

  # load annotation files
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("lumiMouseIDMapping")
  getChipInfo(timecourse.microarrays, species="Mouse") # check annotation 
  tcourse <- addNuID2lumi(timecourse.microarrays, 'MouseWG6_V2_0_R3_11278593_A') 	# convert to nuID
  # obtain affymetrix mouse gene 1.0 st v1 probe ids in exchange for illumina mouse whole genome 6 v2 ids
  ILMN <- nuID2IlluminaID(featureNames(tcourse), species="Mouse", chipVersion='MouseWG6_V2_0_R3_11278593_A')

## Find differentially expressed genes (relative to t0):
  # log scale data
	tcn.log <- log(exprs(tcn))
	# linear scale (rel. changes centered around zero)
	tcn.rel <- (((tcn.log[,1:6] - tcn.log[,1]) + (tcn.log[,7:12] - tcn.log[,7])) /2)
	
	# goi (combined) qualifies if sum over relative changes is larger than Q3 + 3*IQR
	goi.inc		<- rowSums(tcn.rel) > (summary(rowSums(tcn.rel))[5] + 3*IQR(rowSums(tcn.rel)))
	goi.dec		<- rowSums(tcn.rel) < (summary(rowSums(tcn.rel))[2] - 3*IQR(rowSums(tcn.rel)))
	percent.goi <- table(goi.inc)
	percent.inc <- paste(round(100 * percent.goi[2] / sum(percent.goi),2), '%', sep='')
	percent.goi <- table(goi.dec)
	percent.dec <- paste(round(100 * percent.goi[2] / sum(percent.goi),2), '%', sep='')

	# Visualise selection for combined data
	hist(rowSums(tcn.rel), breaks=100, border=F, col='black', xlab="Cumulative Relative Change", main="", las=1, axes=F, xlim=c(-15,15))
	abline(v=(summary(rowSums(tcn.rel))[5] + 3*IQR(rowSums(tcn.rel))), col="red", lwd=2)
	abline(v=(summary(rowSums(tcn.rel))[2] - 3*IQR(rowSums(tcn.rel))), col="red", lwd=2)	
	text( 6, 2500, label=percent.inc, pos=4, col="red")
	text(-6, 2500, label=percent.dec, pos=2, col="red")	
	axis(side=1, las=1, at=c(-10,0,10))
	axis(side=2, at=seq(0, 10000, 5000), las=2)


  # Find differentially expressed genes for individual cell lines
	tcn.rel <- (cbind((tcn.log[,1:6] - tcn.log[,1]), (tcn.log[,7:12] - tcn.log[,7])))
	# goi (E14) qualifies if sum over relative changes is larger than Q3 + 1.5*IQR
	goi.e14.inc <- rowSums(tcn.rel[,1:6])  > (summary(rowSums(tcn.rel[,1:6]))[5] 	+ 3*IQR(rowSums(tcn.rel[,1:6])))
	goi.e14.dec <- rowSums(tcn.rel[,1:6])  < (summary(rowSums(tcn.rel[,1:6]))[2] 	- 3*IQR(rowSums(tcn.rel[,1:6])))	
	# goi (R1) qualifies if sum over relative changes is larger than Q3 + 1.5*IQR
	goi.r1.inc	<- rowSums(tcn.rel[,7:12]) > (summary(rowSums(tcn.rel[,7:12]))[5] + 3*IQR(rowSums(tcn.rel[,7:12])))
	goi.r1.dec	<- rowSums(tcn.rel[,7:12]) < (summary(rowSums(tcn.rel[,7:12]))[2] - 3*IQR(rowSums(tcn.rel[,7:12])))
		
  # Visualise overlap
  venn(list(Combined=unique(fData(tcn)[goi.inc,4]), E14tg2a=unique(fData(tcn)[goi.e14.inc,4]),R1=unique(fData(tcn)[goi.r1.inc,4])))
  venn(list(Combined=unique(fData(tcn)[goi.dec,4]), E14tg2a=unique(fData(tcn)[goi.e14.dec,4]),R1=unique(fData(tcn)[goi.r1.dec,4])))

  # Save GOI list to files
  # combined
  write.table(ILMN[goi.inc, 5], file="./Output/MicroarrayGOI-any-inc-comb.txt", sep="", row.names=F, col.names=F, quote=F)
  write.table(ILMN[goi.dec, 5], file="./Output/MicroarrayGOI-any-dec-comb.txt", sep="", row.names=F, col.names=F, quote=F)
  # E14tg2a
  write.table(ILMN[goi.e14.inc, 5], file="./Output/MicroarrayGOI-any-inc-e14.txt", sep="", row.names=F, col.names=F, quote=F)
  write.table(ILMN[goi.e14.dec, 5], file="./Output/MicroarrayGOI-any-dec-e14.txt", sep="", row.names=F, col.names=F, quote=F)
  # R1
  write.table(ILMN[goi.r1.inc, 5], file="./Output/MicroarrayGOI-any-inc-r1.txt", sep="", row.names=F, col.names=F, quote=F)
  write.table(ILMN[goi.r1.dec, 5], file="./Output/MicroarrayGOI-any-dec-r1.txt", sep="", row.names=F, col.names=F, quote=F)


