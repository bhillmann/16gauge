library(ggplot2)
library(reshape2)
library(cowplot)
library(vegan)

sufs <- c('10','100','1k', '10k','100k','1m', '10m', 'full')
nobs <- c('10', '100', '1000', '10000', '100000', '1000000', '10000000', 'full')

costs_shotgun <- function(x) {51.29 + (10.77/1000000)*x}
costs_16s <- function(x) {23.20 + (100/1000000)*x}

USD = sapply(c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 40000000), costs_shotgun)

USD.labels = apply(rbind(sufs, USD), 2, function(x) {print(x);sprintf('%s ($%.02f)', x['sufs'], as.numeric(x['USD']))})

levels <- c('kingdom','phylum','class', 'order','family','genus', 'species')
prefs <- c('results/simulated-4-16-16/', 'results/simulated_metaphlan2-4-16-16/', 'results/simulated_shogun-4-16-16/')
methods <- c('truth', 'mp2', 'shogun')

predictions <- list()

for(i in 1:length(levels)) {
  level <- levels[i]
  level.rownames <- c()
  preds <- list()

  for (j in 1:length(prefs)) {
    pref <- prefs[j]
    filename <- file.path(pref, sprintf('%s.csv', level))
    pred <- t((read.csv(filename, head=T, row=1, stringsAsFactors = F)))
    preds[[methods[j]]] <- pred
    level.rownames <- c(level.rownames, sapply(rownames(pred), function(x) paste(methods[j], x, sep='-')))
  }
  allnames <- sort(unique(c(unlist(sapply(preds, colnames)))))

  # matrix of observations, each row is a different depth
  predx <- matrix(0, nrow=length(level.rownames), ncol=length(allnames)); rownames(predx) <- level.rownames; colnames(predx) <- allnames


  for(i in 1:length(level.rownames)){
    rowname = level.rownames[i]
    names = strsplit(rowname, "-")[[1]]
    predx[i, colnames(preds[[names[1]]])] <- preds[[names[1]]][names[2],]
  }
  predx[is.na(predx)] <- 0
  predx <- t(apply(predx, 1, function(x) x/sum(x)))
  predx[is.na(predx)] <- 0
  laplace.smooth = 1/length(allnames)
  predx <- as.data.frame(log10(predx + laplace.smooth))
  predx <- as.data.frame(predx)

  #predx$bray <- 1 - unlist(apply(predx, 1, function(x) vegdist(rbind(predx[1,], x))[1]))
  predx$corr.spear <- unlist(apply(predx, 1, function(x) cor(unlist(predx[1,]), x, method='spear')))
  #predx$corr.pear <- unlist(apply(predx, 1, function(x) cor(unlist(predx[1,]), x, method='pear')))

  predx$group <- sapply(rownames(predx), function(x) strsplit(x, '-')[[1]][1])
  predx$depth <- sapply(rownames(predx), function(x) {
    depth = as.numeric(strsplit(x, '[.]')[[1]][2])
    USD.labels[ceiling(log10(depth))]
    })
  predx$depth = factor(predx$depth, levels=USD.labels)
  predictions[[level]] <- predx
}

for (i in 2:length(levels)) {
  level <- levels[i]
  plot <- ggplot(data=predictions[[level]], aes(x=depth, y=corr.spear, group=group)) + geom_smooth(aes(colour=group), method='loess', se=F) + geom_point(aes(shape=group), size=5) + xlab('Sequencing Depth (USD)') + ylab('Spearman Correlation with true metagenome') + theme_light() + ggtitle(paste('Simulated rarefaction accuracy at the', level, 'level'))
  plot
  ggsave(paste('docs/shogun-mp2_', level, '.pdf', sep=''))
}

species = t(predictions[['species']])
species = data.frame(species[1:(dim(species)[1]-3),], stringsAsFactors = F)
species[1:length(colnames(species))] = lapply(species[1:length(colnames(species))], as.numeric)
cor(species$truth.X.40000000, species$mp2.X.40000000)
cor(species$truth.X.40000000, species$shogun.X.40000000, method='spear')

ggplot(data=species, aes(x=shogun.X.40000000, y=truth.X.40000000)) + geom_point() + xlab('Shogun log10 relative abundance') + ylab('Truth log10 relative abundance') + ggtitle('Shogun Species Relative Abundance') + geom_abline(linetype='dotted')
ggsave('docs/shogun-relative_abundance_plot_shogun.pdf')
ggplot(data=species, aes(x=mp2.X.40000000, y=truth.X.40000000)) + geom_point() + xlab('MP2 log10 relative abundance') + ylab('Truth log10 relative abundance') + ggtitle('Metaphlan2 Species Relative Abundance') + geom_abline(linetype='dotted')
ggsave('docs/shogun-relative_abundance_plot_mp2.pdf')

# "There is all this noise", sole basis.

genus = t(predictions$genus)
genus = data.frame(genus[1:(dim(genus)[1]-3),], stringsAsFactors = F)
genus[1:length(colnames(genus))] = lapply(genus[1:length(colnames(genus))], as.numeric)
cor(genus$truth.X.40000000, genus$mp2.X.40000000)
cor(genus$truth.X.40000000, genus$shogun.X.40000000, method='spear')

ggplot(data=genus, aes(x=shogun.X.40000000, y=truth.X.40000000)) + geom_point() + xlab('Shogun log10 relative abundance') + ylab('Truth log10 relative abundance') + ggtitle('Shogun Genus Relative Abundance') + geom_abline(linetype='dotted')
ggsave('docs/shogun-relative_abundance_genus_plot_shogun.pdf')
ggplot(data=genus, aes(x=mp2.X.40000000, y=truth.X.40000000)) + geom_point() + xlab('MP2 log10 relative abundance') + ylab('Truth log10 relative abundance') + ggtitle('Metaphlan2 Genus Relative Abundance') + geom_abline(linetype='dotted')
ggsave('docs/shogun-relative_abundance_genus_plot_mp2.pdf')

mp2.indx = genus$mp2.X.40000000 <= -2.2
write.csv(row.names(genus)[mp2.indx], 'results/mp_genus.csv')

truth.indx = genus$truth.X.40000000 <= -2.2
write.csv(row.names(genus)[truth.indx], 'results/truth_genus.csv')

family = t(predictions$family)
family = data.frame(family[1:(dim(family)[1]-3),], stringsAsFactors = F)
family[1:length(colnames(family))] = lapply(family[1:length(colnames(family))], as.numeric)
cor(family$truth.X.40000000, family$mp2.X.40000000)
cor(family$truth.X.40000000, family$shogun.X.40000000, method='spear')

ggplot(data=family, aes(x=shogun.X.40000000, y=truth.X.40000000)) + geom_point() + xlab('Shogun log10 relative abundance') + ylab('Truth log10 relative abundance') + ggtitle('Shogun Family Relative Abundance') + geom_abline(linetype='dotted')
ggsave('docs/shogun-relative_abundance_family_plot_shogun.pdf')
ggplot(data=family, aes(x=mp2.X.40000000, y=truth.X.40000000)) + geom_point() + xlab('MP2 log10 relative abundance') + ylab('Truth log10 relative abundance') + ggtitle('Metaphlan2 Family Relative Abundance') + geom_abline(linetype='dotted')
ggsave('docs/shogun-relative_abundance_family_plot_mp2.pdf')

genus.lm = lm(mp2.X.40000000~truth.X.40000000, genus)
genus.res = sort(abs(resid(genus.lm)), decreasing = T)[1:10]
write.csv(row.names(genus)[truth.indx], 'results/top_10_resids_mp2_genus.csv')

