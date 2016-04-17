install.packages('biom',repo='http://cran.wustl.edu')

library('biom')
library('data.table')

gg.otus.biom <- fread('data/QIIME/v35_psn_otu.genus.fixed.txt')
gg.otus <- as.matrix(biom_data(gg.otus.biom))
meta = fread('data/metadata/ppAll_V35_map.txt')

# subset for stool
stool.indx = meta$HMPBodySubsit == 'Stool'
meta.stool = meta[stool.indx,]

# grab the names of stool samples
otus.stool.indx = rownames(gg.otus) %in% meta.stool$SampleID
gg.otus.stool = gg.otus[otus.stool.indx,]
