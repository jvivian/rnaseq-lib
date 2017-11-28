import os
import shutil
import textwrap
from multiprocessing import cpu_count
from subprocess import call, Popen, PIPE

import numpy as np
import pandas as pd

from rnaseq_lib.docker import fix_directory_ownership, base_docker_call
from rnaseq_lib.tissues import get_tumor_samples, get_gtex_samples, get_normal_samples, map_genes
from rnaseq_lib.utils import mkdir_p


def log2fc(a, b):
    """
    Calculate the log2 Fold Change between two arrays, floats, or integers
    a and b cannot be, nor contain, values less than 0

    :param (int/float/np.array) a: Value or array
    :param (int/float/np.array) b: Value or array
    :return: L2FC array or value
    :rtype: (int/float/np.array)
    """
    return np.log2(a + 1) - np.log2(b + 1)


def run_deseq2(df_path, tissue, output_dir, gtex=True, cores=None):
    """
    Runs DESeq2 on a specific tissue

    :param str df_path: Path to samples by genes dataframe
    :param str tissue: Tissue to run
    :param str output_dir: Full path to output directory
    :param bool gtex: If True uses GTEx as normal tissue. Otherwise uses TCGA Normal
    :param int cores: Number of cores to use. Defaults to # of cores on machine.
    """
    # Make workspace directories
    work_dir = os.path.join(output_dir, 'work_dir')
    mkdir_p(work_dir)

    # Get samples for tissue
    tumor = get_tumor_samples(tissue)
    normal = get_gtex_samples(tissue) if gtex else get_normal_samples(tissue)

    # Write out vectors
    tissue_vector = os.path.join(work_dir, 'tissue.vector')
    with open(tissue_vector, 'w') as f:
        f.write('\n'.join(tumor + normal))

    disease_vector = os.path.join(work_dir, 'disease.vector')
    with open(disease_vector, 'w') as f:
        f.write('\n'.join(['T' if x in tumor else 'N' for x in tumor + normal]))

    # Write out script
    cores = cores if cores else int(cpu_count())
    script_path = os.path.join(work_dir, 'deseq2.R')
    with open(script_path, 'w') as f:
        f.write(
            textwrap.dedent("""
            library('DESeq2'); library('data.table'); library('BiocParallel')
            register(MulticoreParam({cores}))

            # Argument parsing
            args <- commandArgs(trailingOnly = TRUE)
            df_path <- args[1]
            tissue_path <- args[2]
            disease_path <- args[3]
            tissue <- '{tissue}'
            output_dir <- '/data/'

            # Read in vectors
            tissue_vector <- read.table(tissue_path)$V1
            disease_vector <- read.table(disease_path)$V1

            # Read in table and process
            n <- read.table(df_path, sep='\\t', header=1, row.names=1, check.names=FALSE)
            sub <- n[, colnames(n)%in%tissue_vector]
            setcolorder(sub, as.character(tissue_vector))

            # Preprocessing
            countData <- round(sub)
            colData <- data.frame(disease=disease_vector, row.names=colnames(countData))
            y <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ disease)

            # Run DESeq2
            y <- DESeq(y, parallel=TRUE)
            res <- results(y, parallel=TRUE)
            summary(res)

            # Write out table
            resOrdered <- res[order(res$padj),]
            res_name <- paste(tissue, '.tsv', sep='')
            res_path <- paste(output_dir, res_name, sep='/')
            write.table(as.data.frame(resOrdered), file=res_path, col.names=NA, sep='\\t',  quote=FALSE)

            # MA Plot
            ma_name <- paste(tissue, '-MA.pdf', sep='')
            ma_path <- paste(output_dir, ma_name, sep='/')
            pdf(ma_path, width=7, height=7)
            plotMA(res, main='DESeq2')
            dev.off()

            # Dispersion Plot
            disp_name <- paste(tissue, '-dispersion.pdf', sep='')
            disp_path <- paste(output_dir, disp_name, sep='/')
            pdf(disp_path, width=7, height=7)
            plotDispEsts( y, ylim = c(1e-6, 1e1) )
            dev.off()

            # PVal Hist
            hist_name <- paste(tissue, '-pval-hist.pdf', sep='')
            hist_path <- paste(output_dir, hist_name, sep='/')
            pdf(hist_path, width=7, height=7)
            hist( res$pvalue, breaks=20, col="grey" )
            dev.off()

            # Ratios plots
            qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
            bins <- cut( res$baseMean, qs )
            levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
            ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
            ratio_name <- paste(tissue, '-ratios.pdf', sep='')
            ratio_path <- paste(output_dir, ratio_name, sep='/')
            pdf(ratio_path, width=7, height=7)
            barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
            dev.off()                                           
            """.format(cores=cores, tissue=tissue)))

    # Call DESeq2
    docker_parameters = ['docker', 'run',
                         '-v', '{}:/data'.format(output_dir),
                         '-v', '{}:/df'.format(os.path.dirname(df_path)),
                         'jvivian/deseq2']

    parameters = ['/data/work_dir/deseq2.R',
                  '/df/{}'.format(os.path.basename(df_path)),
                  '/data/{}'.format(os.path.join('work_dir', 'tissue.vector')),
                  '/data/{}'.format(os.path.join('work_dir', 'disease.vector'))]

    print '\nCalling: {}\n'.format(' '.join(docker_parameters + parameters))
    p = Popen(docker_parameters + parameters, stderr=PIPE, stdout=PIPE)
    out, err = p.communicate()
    if out or err:
        print out
        print err

    # Fix output of files
    fix_directory_ownership(output_dir=output_dir, tool='jvivian/deseq2')

    # Add gene names to output
    output_tsv = os.path.join(output_dir, '{}.tsv'.format(tissue))
    df = pd.read_csv(output_tsv, index_col=0, sep='\t')
    df.index = map_genes(df.index)
    df.to_csv(output_tsv, sep='\t')

    # Clean up
    shutil.rmtree(work_dir)


def deseq2_normalize(df_path,
                     output_dir='.',
                     map_gene_names=True,
                     clean_workdir=True,
                     normalize_fn=None,
                     suffix='deseq2-normalized.tsv'):
    """
    Accepts a gene by sample expression matrix normalized values with DESeq2
    Output filename: <INPUT>.deseq2-normalized.tsv

    :param str df_path: Path to input expression gene by sample dataframe
    :param str output_dir: Output directory
    :param bool map_gene_names: If True, maps gene IDs to gene names
    :param bool clean_workdir: If True, deletes temporary work directory
    :param fn normalize_fn: Pass a function to apply to the dataframe before normalization. e.g. lambda x: 2**x + 1
    :param str suffix: Suffix added after df_paths basename
    :return: Path to normalized dataframe
    :rtype: str
    """
    # Make workspace directory
    output_dir = os.path.abspath(output_dir)
    work_dir = os.path.join(output_dir, 'work_dir')
    mkdir_p(work_dir)

    # Write out pseudo-vector
    # Tested that results are the same regardless of vector grouping
    samples = [x.strip() for x in open(df_path, 'r').readline().split()[1:]]
    tissue_vector = os.path.join(work_dir, 'tissue.vector')
    with open(tissue_vector, 'w') as f:
        f.write('\n'.join(samples))

    if normalize_fn:
        df = pd.read_csv(df_path, sep='\t', index_col=0)
        df = df.apply(normalize_fn)
        df_path = os.path.join(os.path.dirname(df_path), 'processed.' + os.path.basename(df_path))
        df.to_csv(df_path, sep='\t')

    # Write normalization script
    output_path = os.path.join(output_dir, os.path.basename(df_path).split('.')[0] + suffix)
    script_path = os.path.join(work_dir, 'deseq2.R')
    with open(script_path, 'w') as f:
        f.write(
            textwrap.dedent("""
            suppressMessages(library('DESeq2'))

            # Argument parsing
            output_dir <- '/data/'

            # Read in pseudo-vector
            tissue_vector <- read.table('/data/work_dir/tissue.vector')$V1

            # Read in table and process
            print("Reading in dataframe")
            n <- read.table('/df/{df_path}', sep='\\t', header=1, row.names=1, check.names=FALSE)

            # Preprocessing
            print("Rounding data to integers")
            countData <- round(n)

            print("Creating DESeq2 Dataset Object")
            colData <- data.frame(tissue=tissue_vector, row.names=colnames(countData))
            dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ tissue)

            # Estimate size factors
            print("Estimating Size Factors")
            dds <- estimateSizeFactors(dds)

            # Extract Normalized Counts and Write
            print("Extracting normalized counts")
            norm <- counts(dds, normalized = TRUE)

            print("Writing output: {output_name}")
            write.table(norm, file="/data/{output_name}", sep='\\t', quote=F, dec='.', col.names=NA) 
            """.format(df_path=os.path.basename(df_path), output_name=os.path.basename(output_path))))

    # Call Docker
    base_params = base_docker_call(os.path.dirname(output_path))
    parameters = base_params + ['-v', '{}:/df'.format(os.path.abspath(os.path.dirname(df_path))),
                                'jvivian/deseq2',
                                '/data/work_dir/deseq2.R']

    print '\nCalling: {}\n'.format(' '.join(parameters))
    call(parameters)

    # Fix output of files
    fix_directory_ownership(output_dir=output_dir, tool='jvivian/deseq2')

    # Map gene IDs to gene names
    if map_gene_names:
        df = pd.read_csv(output_path, index_col=0, sep='\t')
        df.index = map_genes(df.index)
        df.to_csv(output_path, sep='\t')

    # Clean up
    if clean_workdir:
        shutil.rmtree(work_dir)

    return output_path
