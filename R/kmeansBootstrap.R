





par(mar=c(1,1,2,1), pty="s")
celllines <- list(E14=E14, R1=R1)
clusters <- list()
clusters.bs <- list()
bootstraps=5

for (line in c('E14','R1')) {
  
  # Principal Component Analysis on binary expr data (ET > 2)
  pca.tc <- princomp(t(timecourse.log2[,celllines[[line]]]>2))
  
  # Extract scores
  xs <- pca.tc$scores[,1]
  ys <- pca.tc$scores[,2]
  xsign <- -sign(xs[1])
  ysign <- -sign(ys[1])
  
  # flip axis to have 0h in bottom left corner (uses first cell as reference)
  xs <-  xsign* xs
  ys <-  ysign* ys
  
  # Percent variance explained:
  varexpl <- (pca.tc$sdev^2/sum(pca.tc$sdev^2))[1:10]
  
  # Plot Scores
  plot(x=xs, y=ys, las=1, bty='n', main=line, xlab='', ylab='',
       col=as.vector(sample.colours[celllines[[line]]]), pch=19, cex=1.5, axes=F)
  mtext(paste0("PC1 (", signif(varexpl[1]*100,3),"%)"), side=1)
  mtext(paste0("PC2 (", signif(varexpl[2]*100,3),"%)"), side=2)
  lines(tapply(xs, sample.factor[celllines[[line]]], mean), 
        tapply(ys, sample.factor[celllines[[line]]], mean), col='#555555', lwd=4)
  
  # Colour by kmeans
  clu <- kmeans(cbind(xs,ys), centers = 3, nstart = 10)$cluster
  clu <- factor(clu, levels=order(tapply(xs, INDEX = clu, FUN = mean)), labels=1:3)
  clusters[[line]] <- as.numeric(as.vector(clu))
  
  # Viz PCA with cluster annotation
  plot(x=xs, y=ys, las=1, bty='n', main=line,
       col=clusters[[line]]+5, pch=19, cex=1.5, axes=F)
  mtext(paste0("PC1 (", signif(varexpl[1]*100,3),"%)"), side=1)
  mtext(paste0("PC2 (", signif(varexpl[2]*100,3),"%)"), side=2)
  
  
  # bootstrapping
  n = sum(celllines[[line]]);
  clusters.bs[[line]] <- array(dim = c(3, 7, bootstraps), dimnames=list(state=c('ESC', 'EPI','NPC'),
                                                                        time=levels(sample.factor),
                                                                        bootstrap=1:bootstraps))
  for (strap in 1:bootstraps){
      ixs <- sample(1:n, size = n, replace = T)
      clu.bs <- kmeans(cbind(xs,ys)[ixs,], centers = 3, nstart = 10)$cluster
      clu.bs <- factor(clu.bs, levels=order(tapply(xs[ixs], INDEX = clu.bs, FUN = mean)), labels=1:3)
      clu.mat <- table(clu.bs, sample.factor[celllines[[line]]][ixs])
      clusters.bs[[line]][, , strap] <- (clu.mat / (matrix(clu.mat %>% colSums, nrow=3, ncol=7, byrow=T)))
    }

}
