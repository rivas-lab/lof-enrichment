import argparse
import cPickle
import os
import sys

import pandas as pd

def describe_overlap(genes, piemm_results, afr_ss, amr_ss, lof_gene):
    # Read summary stats
    afr_ss = pd.read_table(afr_ss, index_col=0)
    amr_ss = pd.read_table(amr_ss, index_col=0)

    # Read per-gene LoF data.
    lof_gene = pd.read_table(lof_gene, index_col=0)

    # Read PIEMM results.
    piemm_res = cPickle.load(open(piemm_results))

    afr_gamma = piemm_res['AFR'][2]
    amr_gamma = piemm_res['AMR'][2]
    afr_enriched = set(afr_gamma.idxmax(axis=1)[afr_gamma.idxmax(axis=1) ==
                                                1].index)
    amr_enriched = set(amr_gamma.idxmax(axis=1)[amr_gamma.idxmax(axis=1) ==
                                                1].index)
    
    # Read Gencode gene info.
    gtable = pd.read_table('gencode.v19.annotation.table.tsv', index_col=0)
    gtable['gencode_id'] = gtable.index
    gtable.index = [x.split('.')[0] for x in gtable.index]
        
    # Read gene list
    gene_list = pd.read_table(genes, header=None, squeeze=True)
    total_genes = gene_list.shape[0]

    # If the following is true, then the list is probably gene symbols.
    if (len(set(gene_list) & set(gtable.index)) < 
        len(set(gene_list) & set(gtable.gene_name))):
        ind = []
        for i in gene_list.index:
            if len(gtable[gtable.gene_name == gene_list[i]].index) > 0:
                ind.append(gtable[gtable.gene_name == gene_list[i]].index[0])
            else:
                ind.append('NaN')
        gene_list.index = ind
        gene_list = gene_list.drop('NaN')
    else:
        gene_list.index = gene_list.values
    gene_list.drop_duplicates()
    
    sys.stdout.write('{:,} of {:,} provided genes recognized.\n'.format(
        gene_list.shape[0], total_genes))
    sys.stdout.write('{:,} of {:,} provided genes tested for enrichment.\n'.format(
        len(set(gene_list.index) & set(afr_gamma.index)), total_genes))
    sys.stdout.write('{:,} genes enriched for LoF variants in AFR.\n'.format(
        len(afr_enriched & set(gene_list.index))))
    sys.stdout.write('{:,} genes enriched for LoF variants in AMR.\n'.format(
        len(amr_enriched & set(gene_list.index))))
    sys.stdout.write('{:,} shared enriched genes for AFR and AMR.\n'.format(
        len(afr_enriched & amr_enriched &
            set(gene_list.index))))
    
    out = pd.DataFrame(index=(afr_enriched | amr_enriched) & set(gene_list.index),
                                  columns=['symbol', 'AFR_beta', 'AMR_beta'])
    out['symbol'] = gtable.ix[out.index, 'gene_name']
    out.ix[afr_enriched & set(out.index), 'AFR_beta'] = afr_ss.ix[afr_enriched & set(out.index), 'beta']
    out.ix[amr_enriched & set(out.index), 'AMR_beta'] = amr_ss.ix[amr_enriched & set(out.index), 'beta']
    out = out.sort_values(by='symbol')
    out['AF_AFR'] = lof_gene.AF_AFR.ix[out.index]
    out['AF_AMR'] = lof_gene.AF_AMR.ix[out.index]
    out['AF_NFE'] = lof_gene.AF_NFE.ix[out.index]
    out['AC_AFR'] = lof_gene.AC_AFR.ix[out.index]
    out['AC_AMR'] = lof_gene.AC_AMR.ix[out.index]
    out['AC_NFE'] = lof_gene.AC_NFE.ix[out.index]
    out_fn = os.path.splitext(os.path.split(genes)[1])[0]
    out.to_csv(out_fn + '_lof_enriched.tsv', sep='\t')

def main():
    piemm_results = 'piemm_res_01_05_05.pickle'
    afr_ss = 'AFR_ss.tsv'
    amr_ss = 'AMR_ss.tsv'
    lof_gene = 'lof_gene.tsv'
    parser = argparse.ArgumentParser(
        description=('Query a gene list of Ensembl IDs or gene symbols to '
                     'see which genes are enriched for LoF variants in AFR '
                     'or AMR populations. Outputs a TSV file with overlap '
                     'information.'))
    parser.add_argument('genes', help=('Gene list with one gene per line and '
                                       'no header. Must be Ensembl gene IDs '
                                       'or gene symbols.'))
    parser.add_argument('--piemm_results', 
                        help=('Pickle file containing PIEMM results. Default: '
                              '{}.'.format(piemm_results)), 
                        default=piemm_results)
    parser.add_argument('--afr_ss', 
                        help=('Summary stats for AFR. Default: '
                              '{}.'.format(afr_ss)), 
                        default=afr_ss)
    parser.add_argument('--amr_ss', 
                        help=('Summary stats for AMR. Default: '
                              '{}.'.format(amr_ss)), 
                        default=amr_ss)
    parser.add_argument('--lof_gene', 
                        help=('File with LoF counts per gene. Default: '
                              '{}.'.format(lof_gene)), 
                        default=lof_gene)
    args = parser.parse_args()

    genes = args.genes
    piemm_results = args.piemm_results

    describe_overlap(genes, piemm_results, afr_ss, amr_ss, lof_gene)

if __name__ == '__main__':
    main()
