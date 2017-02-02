library(vegan)
library(ggplot2)
library(optparse)

#Read arguments
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="compounds file name", metavar="character"),
  make_option(c("-m", "--map"), type="character", default=NULL, 
              help="mapping file name", metavar="character"),
	make_option(c("-g", "--grouping"), type="character", default=NULL, 
              help="grouping column from mapping file to color PCoA plot by.", metavar="character"),
	make_option(c("-c", "--category"), type="character", default=NULL, 
              help="category column from mapping file to run ANOSIM on. This column should only have two values.", metavar="character")
	)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# #Get data; clean up
d = read.csv(opt$file)
rownames(d) = paste0(d$Compound, d$Mass)
d$Compound = NULL
d$Mass = NULL
d$X = NULL
data = d

#Read mapping file; order correctly
map = read.csv(opt$map, sep="\t")
map = map[order(map$file),]
map$file = as.character(map$file)
map$file = gsub('/', '.', map$file)
data = as.matrix(data[ , order(names(data))])

# #Create two distance matrices: jaccard and bray-curtis
distances = vegdist(log(1+t(data)), method="jaccard")
distances2 = vegdist(log(1+t(data)), method="bray")

# #Run PCoA on log-normalized distances 
pcoa = cmdscale(distances)
pcoa2 = cmdscale(distances2)

#Output results
png("pcoa_jaccard.png", width=600, height=400)
qplot(pcoa[,1], pcoa[,2], colour=map[,opt$grouping]) + geom_line()
dev.off()

png("pcoa_bray.png", width=600, height=400)
qplot(pcoa2[,1], pcoa2[,2], colour=map[,opt$grouping]) + geom_line()
dev.off()

#Run ANOSIM for categorical variables
print("Results of ANOSIM test on category:")
anosim(distances, map[,opt$category])