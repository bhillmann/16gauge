library(ggplot2)
library(reshape2)
library(cowplot)

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
    pred <- t((read.csv(filename, head=T, row=1)))
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

  predx$corr.spear <- unlist(apply(predx, 1, function(x) cor(unlist(predx[1,]), x, method='spear')))
  predx$corr.pear <- unlist(apply(predx, 1, function(x) cor(unlist(predx[1,]), x, method='pear')))
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
