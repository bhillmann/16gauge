library('data.table')

sam = fread('results/bt2_sam-4-11-16/bt2_sam_files.txt')
sam.indx <- sapply(sam$picrust, function(x) {substr(x, 1, 1) == 'S'})
shotgun.samples <- sam$picrust[sam.indx]
shotgun.samples = sapply(shotgun.samples, function(x) {substr(x, 1, 9)})

meta = fread('data/metadata/ppAll_V35_map.txt')

# subset for stool
stool.indx = meta$HMPBodySubsit == 'Stool'
meta.stool = meta[stool.indx,]

# Number of Samples in Both
sum(meta.stool$SRS_SampleID %in% shotgun.samples)

# Save the metadata file
shotgun.samples.indx = meta.stool$SRS_SampleID %in% shotgun.samples
meta.stool = meta.stool[shotgun.samples.indx]
meta.stool
