library(ggplot2)
library(reshape2)

sufs <- c('10','100','1k', '10k','100k','1m', '10m', 'full')
nobs <- c('10', '100', '1000', '10000', '100000', '1000000', '10000000', 'full')

costs_shotgun <- function(x) {51.29 + (10.77/1000000)*x}
costs_16s <- function(x) {23.20 + (100/1000000)*x}

USD = sapply(c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 40000000), costs_shotgun)

USD.labels = apply(rbind(sufs, USD), 2, function(x) {print(x);sprintf('%s ($%.02f)', x['sufs'], as.numeric(x['USD']))})

depths <- nobs / .06

preds <- list()
obss <- list()
pre = 'results/predictions-4-16-16'

full.obs = as.matrix(read.table(file.path(pre, 'SRS011271.consensus.full.csv'), sep=',', head=T, row=1))

for(i in 1:length(sufs)) {
    nob <- nobs[i]
    suf <- sufs[i]
    filename <- file.path(pre, sprintf('SRS011271.intersection.%s.csv', nob))
    pred <- as.matrix(read.table(filename, sep=',', head=T, row=1)); tmp <- colnames(pred); pred <- as.numeric(pred); names(pred) <- tmp
    filename <- file.path(pre, sprintf('SRS011271.consensus.%s.csv',nob))
    obs <- as.matrix(read.table(filename, sep=',', head=T, row=1)); tmp <- colnames(obs); obs <- as.numeric(obs); names(obs) <- tmp
    allnames <- union(names(full.obs), union(names(obs), names(pred)))

    preds[[suf]] <- pred
    obss[[suf]] <- obs
}

# obss[['full']] <- full.obs

allnames <- sort(unique(c(unlist(sapply(preds, names)), unlist(sapply(obss, names)), names(full.obs))))

# matrix of observations, each row is a different depth
predx <- matrix(0,nrow=length(preds), ncol=length(allnames)); rownames(predx) <- sufs; colnames(predx) <- allnames
obsx <- matrix(0,nrow=length(obss), ncol=length(allnames)); rownames(obsx) <- sufs; colnames(obsx) <- allnames

truth <- matrix(0, nrow=1, ncol=length(allnames)); rownames(truth) <- 'full'; colnames(truth) <- allnames

for(i in 1:length(preds)){
    predx[i, names(preds[[i]])] <- preds[[i]]
    obsx[i, names(obss[[i]])] <- obss[[i]]
}

truth[1, allnames %in% colnames(full.obs)] <- full.obs[1,]

laplace.smooth = 1/length(allnames)

# add eps to obsx and predx
predx.log <- predx
obsx.log <- obsx
truth.log <- truth
for(i in 1:nrow(predx.log)) predx.log[i,] <- predx.log[i,] + laplace.smooth
for(i in 1:nrow(obsx.log)) obsx.log[i,] <- obsx.log[i,] + laplace.smooth
truth.log[1, ] <- truth.log[1, ] + laplace.smooth
predx.log <- log10(predx.log)
obsx.log <- log10(obsx.log)
truth.log <- log10(truth.log)

cc.pred <- apply(predx, 1, function(xx) cor(xx, truth[1,], method='spear'))
cc.obs <- apply(obsx,1,function(xx) cor(xx, truth[1,], method='spear'))
cc.pred.log <- apply(predx.log, 1, function(xx) cor(xx, truth.log[1,], method='spear'))

cc.pred.pearson <- apply(predx,1,function(xx) cor(xx, truth[1,],method='pear'))
cc.obs.pearson <- apply(obsx,1,function(xx) cor(xx, truth[1,],method='pear'))

# ignore unobserved values at each depth
ignore.ix <- truth[1,] == 0
cc.pred.ignore.no.unobs <- apply(predx,1,function(xx) cor(xx[!ignore.ix],truth[1, !ignore.ix],method='spear'))
cc.pred.ignore.no.unobs.log <- apply(predx.log,1,function(xx) cor(xx[!ignore.ix],truth.log[1, !ignore.ix],method='pear'))


cc.pred.pearson.no.unobs <- apply(predx,1,function(xx) cor(xx[!ignore.ix],truth[1, !ignore.ix],method='pear'))
cc.obs.pearson.no.unobs <- apply(predx.log,1,function(xx) cor(xx[!ignore.ix],truth.log[1, !ignore.ix],method='pear'))


obsx <- obsx[,order(obsx[5,])]
predx <- predx[,colnames(obsx)]
obsx.log <- obsx.log[,order(obsx.log[5,])]
predx.log <- predx.log[,colnames(obsx.log)]

fit.res <- resid(lm(predx[5,] ~ obsx[5,]))
# make bins of 100
sd.bins <- c(seq(1, 5701,100),length(obsx[5,]))
# calculate in each bin
sd.bin <- numeric(length(sd.bins) - 1)
for(i in 1:(length(sd.bins)-1)) sd.bin[i] <- sd(fit.res[sd.bins[i]:sd.bins[i+1]])

# calculate mean obsx value in each bin
obs.bin <- numeric(length(sd.bins) - 1)
for(i in 1:(length(sd.bins)-1)) obs.bin[i] <- sd(obsx[5,][sd.bins[i]:sd.bins[i+1]])

dat = data.frame(Shogun=cc.pred.ignore.no.unobs, Observed=cc.obs, depths=USD.labels, PICRUSt=rep(.82, length(USD.labels)))
dat = melt(dat, id.vars=c('depths'))
dat$depths = factor(dat$depths, levels=USD.labels)

ggplot(data=dat, aes(x=depths, y=value, group=variable)) + geom_smooth(aes(colour=variable), method='loess', se=F) + geom_point(aes(shape=variable), size=5) + xlab('Sequencing Depth (USD)') + ylab('Spearman Correlation with true metagenome') + ggtitle("Subsampling Predicted Potential Functional Repertoire of Metagenome")
ggsave("docs/shogun-functional_predicted.pdf")
dim(obsx)
obsx.indx = obsx[dim(obsx)[1],] > 0
predx.genes = predx[,obsx.indx]

predx.counts = apply(predx.genes, 1, function(x) sum(x > 0)/sum(obsx.indx))

predx.indx = predx[dim(predx)[1],] > 0
obsx.genes = obsx[,predx.indx]

obsx.counts = apply(obsx.genes, 1, function(x) sum(x > 0)/sum(predx.indx))

counts = data.frame(Observed=obsx.counts, Shogun=predx.counts, labels=USD.labels, depth=c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 40000000))
counts$labels = factor(counts$labels, levels=USD.labels)
counts = counts[-(dim(counts)-1),]

line = lm(counts$Observed~log10(counts$depth))

ggplot(data=counts, aes(x=labels, y=Observed)) + geom_point(shape=2, size=5) + ylab("Fraction of Predicted Metagenome Observed") + xlab("Sequencing Depth (USD)") + ggtitle('Fraction Predicted Genome Missed')
#ggplot(data=counts, aes(x=labels, y=Observed)) + geom_point(shape=2, size=5) + ylab("Fraction of predicted metagenome observed") + xlab("Sequencing Depth (USD)") + geom_smooth(aes(x=log10(depth), y=Observed), method='glm', se=F)
ggsave('docs/shogun-functional_rarefaction.pdf')
