# edgeR also has limma
suppressMessages(suppressWarnings(library("edgeR", character.only=T, warn.conflicts = F, quietly = T)))

# args from command line:
args<-commandArgs(TRUE)
RAW_COUNT_MATRIX<-args[1]
ANN_FILE<-args[2]
ORIG_ANN_COL<-args[3]
CONDITION_A<-args[4]
CONDITION_B<-args[5]
CUTOFF<-as.integer(args[6])

# the "prefix" for the output file
OUTPUT_FILE_BASE <- 'limma_voom_results'
OUTPUT_NORMALIZED_COUNTS_BASE <- 'voom_cpm_normalized_counts'

# change the working directory to co-locate with the counts file:
working_dir <- dirname(RAW_COUNT_MATRIX)
setwd(working_dir)

# convert names in case they are not "proper" for R
ANN_COL <- make.names(ORIG_ANN_COL)

# create a string to identify the contrast:
contrast_str = paste0(ANN_COL, '_', CONDITION_B, '_vs_', CONDITION_A)

# read the annotation file:
annotations <- read.table(ANN_FILE, sep='\t', header = T, row.names = 1)

if (!(ANN_COL %in% colnames(annotations))) {
    message(sprintf('The column "%s" was not found in your annotation file. Please check your inputs. Note that this input is case-sensitive and must match exactly.', ORIG_ANN_COL))
    quit(status=1)
}

# filter to only those samples which match the conditions
base_samples_filter <- annotations[, ANN_COL] == CONDITION_A
exp_samples_filter <- annotations[, ANN_COL] == CONDITION_B
orig_base_samples <- rownames(annotations)[base_samples_filter]
orig_exp_samples <- rownames(annotations)[exp_samples_filter]
base_samples <- make.names(orig_base_samples)
exp_samples <- make.names(orig_exp_samples)

if(
    (length(base_samples) == 0)
    ||
    (length(exp_samples) == 0)
){
    message(sprintf('One or both of your sample sets contained zero samples. Check that both %s and %s are valid values in column "%s".', CONDITION_A, CONDITION_B, ORIG_ANN_COL))
    quit(status=1)
}

if(
    (length(base_samples) < 2)
    ||
    (length(exp_samples) < 2)
){
    message(sprintf('One or both of your sample sets contained fewer than 2 samples. To perform differential expression analysis, replicates are required. Check that both conditions %s and %s have 2 or more samples each.', CONDITION_A, CONDITION_B))
    quit(status=1)
}

intersection_list = intersect(base_samples, exp_samples)

if (length(intersection_list) > 0){
    sample_list = paste0(intersection_list, collapse=',')
    message(paste(
       'The following samples were in both contrast groups. Fix this and try again: ',
       sample_list
    ))
    quit(status=1)
}

# create a full list. Note that we have to track the original sample names
# so we can map back at the end
all_samples <- c(base_samples, exp_samples)
original_sample_names <- c(orig_base_samples, orig_exp_samples)
colname_mapping = data.frame(
    orig_names = original_sample_names,
    row.names=all_samples,
    stringsAsFactors=F)

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

# Check that after subsetting samples, we actually have more than one 'condition' represented.
# This also covers the edge case where only a single sample name was valid. In that case, the 
# column subset below produces an array and `dim(count_data)[2]` throws an error.
fl <- length(levels(as.factor(annotations$condition)))
if(fl < 2){
    message(sprintf('After subsetting the matrix for the samples of interest (%d found), only one cohort of samples was present. Please double-check your inputs or sample names.',  dim(annotations)[1]))
    quit(status=1)
}

# subset to only keep samples corresponding to the current groups in the count_data dataframe
count_data <- count_data[,annotations$sample]

if (dim(count_data)[2] == 0){
    message('After subsetting the matrix for the samples of interest, the matrix was empty. Please check the input samples and matrix')
    quit(status=1)
}

lowerthan <- rowSums(count_data < CUTOFF) == dim(count_data)[2]
count_data <- count_data[!lowerthan,]

# if there is nothing left, exit...
if (dim(count_data)[1] == 0){
    message('After subsetting the matrix to remove lowly expressed genes, the matrix was empty. Please check the expression cutoff and your matrix')
    quit(status=1)
}

# create a DGEList object for working with the counts, annotations, etc.
d <- DGEList(count_data)
d <- calcNormFactors(d)

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
t<-topTable(efit, coef=paste0('conditions',CONDITION_B), sort.by='P', n=Inf)

# Produce the normalized counts-per-million values
nc <- cpm(d, normalized.lib.sizes = TRUE, log = FALSE)
nc_cols = colnames(nc)
remapped_cols = colname_mapping[nc_cols, 'orig_names']
colnames(nc) = remapped_cols

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
cols[which(names(m) == 'AveExpr')] = 'overall_mean'
cols[which(names(m) == 't')] = 'statistic'
colnames(m) <- cols

# the merge ends up changing the sort order, so resort:
m = m[order(m$pvalue),]

# write the concatenated table to file
contrast_str <- paste0(CONDITION_B, '_vs_', CONDITION_A)
f1 <- paste0(OUTPUT_FILE_BASE, '.', contrast_str, '.tsv')
f1 <- paste(working_dir, f1, sep='/')
write.table(m, f1, sep='\t', quote=F)

f2 <- paste(OUTPUT_NORMALIZED_COUNTS_BASE, contrast_str, 'tsv', sep='.')
f2 <- paste(working_dir, f2, sep='/')
write.table(nc, f2, sep='\t', quote=F, row.names=T)

# create the expected outputs file:
json_str = paste0(
       '{"dge_results":"', f1, '",',
       '"normalized_counts":"', f2, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)


