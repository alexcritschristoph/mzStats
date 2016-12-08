library(vegan)
library(ggplot2)
#Read data file
args = commandArgs(trailingOnly=TRUE)
d = read.csv(args[1])
d$X = NULL
data = d[,2:length(d)]
#Read mapping file
map = read.csv(args[2])

#Create two distance matrices: jaccard and bray-curtis
distances = vegdist(log(1+t(data)), method="jaccard")
distances2 = vegdist(log(1+t(data)), method="bray")

#Run PCoA on log-normalized distances 
pcoa = cmdscale(distances)
pcoa2 = cmdscale(distances2)

#Output results
png("jaccard.png")
qplot(pcoa[,1], pcoa[,2], colour=days)
dev.off()

png("bray.png")
qplot(pcoa2[,1], pcoa2[,2], colour=days)
dev.off()

#Run clustering
h = hclust(distances)
h2 = hclust(distances2)

#Output results
png("jaccard_clust.png")
plot(h)
dev.off()

png("bray_clust.png")
plot(h2)
dev.off()

group = c("60","60","60","60","90","90","90","90")

#Run ANOSIM for categorical variables
anosim(distances, group)

#Run differential abundance/presence testing

pvalues = c()
for (i in 1:nrow(d2)){
	w = wilcox.test(dnew[,i]~group)
	pvalues = c(pvalues, w$p.value)
}

print(d2[,which(pvalues<0.05)])


#Create MZxml of differentially abundant compounds

