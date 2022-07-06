#prcomp(t(input.data), scale. = T)$x
norm.mtx <- norm.peak[,c(4:ncol(norm.peak))]
norm.mtx %>% dim()
input.data <- norm.mtx
pca = as.data.frame(prcomp(t(input.data), scale. = T)$x)
pca %>% str()
pca
PCA <- prcomp(t(input.data),scale.=T)
PCA <- (PCA$sdev/sum(PCA$sdev))[1:2]
PCA_percentage <- paste(round(PCA,4)*10^2,'%',sep = '')

ggplot(pca, aes(PC1,PC2,label=rownames(pca), 
                color= substr(rownames(pca),1,2))) + geom_point(size=4) +
  xlab(paste0('PC1: ',PCA_percentage[1])) +
  ylab(paste0('PC2: ',PCA_percentage[2])) + theme(legend.position = 'none')

ggplot(pca, aes(PC1,PC2,label=rownames(pca), 
                color= substr(rownames(pca),1,2))) + geom_point(size=4, alpha=0.5) +
  xlab(paste0('PC1: ',PCA_percentage[1])) +
  ylab(paste0('PC2: ',PCA_percentage[2])) +geom_text()+ theme(legend.position = 'none')
