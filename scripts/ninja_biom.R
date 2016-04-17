library('data.table')
library('biom')

gg.otus.biom <- read_biom('data/ninja_otus_16s/ninja_otutable.biom')
gg.otus <- as.matrix(biom_data(gg.otus.biom))
meta = fread('data/metadata/ppAll_V35_map.txt')
