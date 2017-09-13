import os
import textwrap
from multiprocessing import cpu_count
from subprocess import call

import shutil

from rnaseq_lib.tissues import get_tumor_samples, get_gtex_samples, get_normal_samples
from rnaseq_lib.utils import mkdir_p


def run_deseq2(df_path, tissue, output_dir, gtex=True):
    # Make workspace directories
    work_dir = os.path.join(output_dir, 'deseq2-results', 'work_dir')
    mkdir_p(work_dir)

    # Get samples for tissue
    tumor = [x.replace('-', '.') for x in get_tumor_samples(tissue)]
    normal = get_gtex_samples(tissue) if gtex else get_normal_samples(tissue)
    normal = [x.replace('-', '.') for x in normal]

    # Write out vectors
    tissue_vector = os.path.join(work_dir, 'tissue.vector')
    with open(tissue_vector, 'w') as f:
        f.write('\n'.join(tumor + normal))

    disease_vector = os.path.join(work_dir, 'disease.vector')
    with open(disease_vector, 'w') as f:
        f.write('\n'.join(['T' if x in tumor else 'N' for x in tumor + normal]))

    # Write out script
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
            tissue <- {tissue}
            output_dir <- '/data/deseq2-results'
            
            # Read in table and process
            n <- read.table(df_path, sep='\\t', header=1, row.names=1)
            sub <- n[, colnames(n)%in%vector]
            setcolorder(sub, as.character(vector))
            
            # Read in vectors
            vector <- read.table(vector_path)$V1
            disease_vector <- read.table(disease_path)$V1
            
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
            """.format(cores=cpu_count(), tissue=tissue)))

    # Call DESeq2
    docker_parameters = ['docker', 'run',
                         '-v', '{}:/data'.format(output_dir),
                         '-v', '{}:/df'.format(os.path.dirname(df_path)),
                         'jvivian/deseq2']

    parameters = ['/df/{}'.format(os.path.basename(df_path)),
                  '/data/{}'.format(os.path.join('deseq2-results', 'work_dir', 'tissue.vector')),
                  '/data/{}'.format(os.path.join('deseq2-results', 'work_dir', 'disease.vector'))]

    print '\nCalling: {}\n'.format(docker_parameters + parameters)
    call(docker_parameters + parameters)
    # shutil.rmtree(work_dir)