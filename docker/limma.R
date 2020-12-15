# edgeR also has limma
suppressMessages(suppressWarnings(library("edgeR", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
BASE_CONDITION_SAMPLES <- args[2]
EXPERIMENTAL_CONDITION_SAMPLES <- args[3]
CONDITION_A<-args[4]
CONDITION_B<-args[5]
CUTOFF<-as.integer(args[6])

# the "prefix" for the output file
OUTPUT_FILE_BASE <- 'limma_voom_results'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# convert names in case they are not "proper" for R
CONDITION_A <- make.names(CONDITION_A)
CONDITION_B <- make.names(CONDITION_B)

# create a string to identify the contrast:
contrast_str <- paste0(CONDITION_B, '_vs_', CONDITION_A)

# the sample names are given as a comma-delimited string. Split them
base_samples <- make.names(strsplit(BASE_CONDITION_SAMPLES, ',')[[1]])
exp_samples <- make.names(strsplit(EXPERIMENTAL_CONDITION_SAMPLES, ',')[[1]])
all_samples <- c(base_samples, exp_samples)

condition_a_list <- rep(CONDITION_A, length(base_samples))
condition_b_list <- rep(CONDITION_B, length(exp_samples))
all_conditions <- c(condition_a_list, condition_b_list)

# full annotation matrix:
annotations <- as.data.frame(cbind(all_samples, all_conditions), stringsAsFactors = F)
colnames(annotations) <- c('sample', 'condition')

# read the raw count matrix, genes as row names:
count_data <- read.table(RAW_COUNT_MATRIX, sep='\t', header = T, row.names = 1, stringsAsFactors = F)

# subset to keep only the samples in the count table.  This is important if the annotation
# file has more samples than the count matrix:
count_mtx_cols <- colnames(count_data)
annotations <- annotations[annotations$sample %in% count_mtx_cols,]

# subset to only keep samples corresponding to the current groups in the count_data dataframe
count_data <- count_data[,annotations$sample]

if (dim(count_data)[2] == 0){
    message('After subsetting the matrix for the samples of interest, the matrix was empty. Please check the input samples and matrix')
    quit(status=1)
}

# create a DGEList object for working with the counts, annotations, etc.
d <- DGEList(count_data)
d <- calcNormFactors(d)

# drop lowly expressed genes
low_genes <- which(apply(cpm(d), 1, max) < CUTOFF)
d <- d[-low_genes,] 

# if there is nothing left, exit...
if (dim(d)[1] == 0){
    message('After subsetting the matrix to remove genes with maximum counts less than the cutoff, the matrix was empty. Please check the expression cutoff and your matrix')
    quit(status=1)
}

# Need to set the condition as a factor since it's going to be used as a design matrix
conditions <- as.factor(annotations$condition)
conditions <- relevel(conditions, ref=CONDITION_A)

# create the model matrix. This model will produce a coefficient that directly gives the difference between the experimental and base condition
# since we set CONDITION_A to be the reference level of the factor above
design <- model.matrix(~conditions)

# run voom to perform the normalization and weights estimation.
v <- voom(d, design)

# now fit this and run the shrinkage:
fitted_data <- lmFit(v, design)
efit <- eBayes(fitted_data)

# Get the top hits, sorted by raw p-value
t<-topTreat(efit, coef=paste0('conditions',CONDITION_B), sort.by='P', n=Inf)

# Produce the normalized counts-per-million values
nc <- cpm(d, normalized.lib.sizes = TRUE, log = FALSE)

# merge the two matrices, which makes it easier to work with on the front-end
m <- merge(t, nc, by="row.names")
rownames(m) <- m[,'Row.names']
drop_cols <- c('Row.names')
m <- m[, !(names(m) %in% drop_cols)]

# rename columns so we don't have to adjust interfaces on the front-end
cols <- colnames(m)
cols[which(names(m) == 'P.Value')] = 'pvalue'
cols[which(names(m) == 'adj.P.Val')] = 'padj'
cols[which(names(m) == 'logFC')] = 'log2FoldChange'
colnames(m) <- cols

contrast_str <- paste0(CONDITION_B, '_vs_', CONDITION_A)
f <- paste0(OUTPUT_FILE_BASE, '.', contrast_str, '.tsv')
write.table(m, f, sep='\t', quote=F)

json_str <- paste0('{"dge_results":"', f, '"}')
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)


