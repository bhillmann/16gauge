library('data.table')

gg.otus.biom <- read_biom('results/ninja_closed-4-5-2016/ninja_otutable.biom')
gg.otus <- as.matrix(biom_data(gg.otus.biom))
meta = fread('data/metadata/ppAll_V35_map.txt')

# subset for stool
stool.indx = meta$HMPBodySubsit == 'Stool'
meta.stool = meta[stool.indx,]

# Choose 20 stool samples for testing methods
set.seed(1)
meta.stool.subset = meta.stool[sapply(meta.stool$SampleID, function(x) {substr(x, 1, 1) != '7'}),]
meta.stool.subset = meta.stool.subset[sample(1:nrow(meta.stool.subset), 20, replace=FALSE),]

# Save the subset files
write.csv(meta.stool.subset, 'data/metadata/meta.stool.subset.csv')
