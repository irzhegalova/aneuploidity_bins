# %%
# general packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import re
import os

# bio packages
import bioframe
from pybedtools import BedTool
import gseapy as gp
from gseapy import barplot, dotplot
from pybiomart import Dataset

# add custom function
def pValue_featureOccurenceInLoop(
    file_loops="results/loops.bed",
    N_shuffle=1000,
    file_features="results/promoters.bed",
    genome_file="data/genome/hg38.chrom.sizes",
    pic_path="results"):

    '''
    Plot feature occurence with permutation test p-value for occurences of features in Hi-C feature

    Parameters
    ----------
    file_loops : str
        Path to the bed file of loops
    N_shuffle : int
        How many times to shuffle the loops
    file_features : str
        File with features to count
    genome_file : str
        Path to the chrom.sizes file used in shuffling the features
    pic_path : str
        Where to save the resulting figures
    '''
    # load loops
    df_loops = BedTool(file_loops)
    # load features
    df_feature = BedTool(file_features)
    # create list to store shuffle counts
    inter_shuffle_list = []
    for i in range(N_shuffle):
        # create shuffles of loops
        loops_shuffled = df_loops.shuffle(g=genome_file, chrom=True, seed=i)
        # intersect features with shuffled loops and count the number of interactions
        inter_shuffle = df_feature.intersect(loops_shuffled, wa=True).to_dataframe().drop_duplicates()
        # append the number to the shuffle list
        inter_shuffle_list.append(inter_shuffle.shape[0])
    # create df with shuffled counts
    inter_shuffle_df = pd.DataFrame({"Shuffle": inter_shuffle_list})
    # intersect real loops with features and count the number of interactions
    inter_Loops = df_feature.intersect(df_loops, wa=True).to_dataframe().drop_duplicates()
    # plot shuffles as a histogram
    fig = sns.histplot(data=inter_shuffle_df, x="Shuffle", kde=True, stat="percent", binwidth=1)
    # add real count as an axline
    fig.axvline(inter_Loops.shape[0], color="red", lw=3)
    # compute shuffle test p-value
    p_value = np.round(
        np.min([len(inter_shuffle_df[inter_shuffle_df['Shuffle'] > inter_Loops.shape[0]]) / N_shuffle,
                len(inter_shuffle_df[inter_shuffle_df['Shuffle'] < inter_Loops.shape[0]]) / N_shuffle]), 3)
    # add p-value to the plot
    if p_value == 0:
        p_value = 1 / N_shuffle
        plt.annotate("p-value < " + str(p_value),
                     xy=(0.95, 0.95), xycoords='axes fraction',
                     bbox=dict(facecolor='pink',
                               alpha=0.5),
                     horizontalalignment='right',
                     fontsize=12)
    else:
        plt.annotate("p-value = " + str(p_value),
                     xy=(0.95, 0.95), xycoords='axes fraction',
                     bbox=dict(facecolor='pink',
                               alpha=0.5),
                     horizontalalignment='right',
                     fontsize=12)
    # save plot
    plt.savefig("%s/%s.pdf"
                % (pic_path, name),
                bbox_inches='tight')
    plt.show()

# %%

os.chdir('~/projects/human/aneuploidity_LADs')
genes_gff = bioframe.read_table('./data/genome/gencode.v43.basic.annotation.gtf', schema='gtf').query('feature == "gene"').sort_values(['chrom', 'start', 'end'], ascending=True).astype({'start': 'int64', 'end': 'int64'})

# add column with gene names
ids = genes_gff['attributes'].str.split('";', expand=True)[0].str.split('.', expand=True)[0]
genes_gff['Gene ID'] = [re.sub('gene_id "', "", x) for x in ids.tolist()]
# add column with gene types
gene_types = genes_gff['attributes'].str.split('";', expand=True)[1]
genes_gff['gene_type'] = [re.sub(' gene_type "', "", x) for x in gene_types.tolist()]

# select only protein coding genes
genes_gff_annGenes = genes_gff.query('gene_type == "protein_coding"')
genes_gff_annGenes.loc[:,['chrom', 'start', 'end']].to_csv('./data/genome/gencode.v43.basic.annotation.genes.bed', index=False, header=False, sep='\t')

# %% HK genes
HK_exons = pd.read_excel('./data/HK_exons.xlsx', sheet_name='HK_exons.both')
HK_genes = genes_gff_ann.loc[genes_gff_ann['Gene name'].isin(HK_exons['Gene Name'].tolist()),:]

# %% load out data
allsums_lads = pd.read_table('./data/allsums_1Mb_new_50kb.tsv')
# sHF18/sPFCH6 (трисомия по 18 хромосоме, деление на норму)
allsums_lads['tri_18'] = allsums_lads.sHF18 / allsums_lads.sPFCH6
# s3471/s3494 (трисомия по 16), 
allsums_lads['tri_16'] = allsums_lads.s3471 / allsums_lads.s3494
# s3475/s3492 (трисомия по 13), 
allsums_lads['tri_13'] = allsums_lads.s3475 / allsums_lads.s3492

allsums_lads['norm1'] = allsums_lads.s3492 / allsums_lads.s3525
allsums_lads['norm2'] = allsums_lads.s3492 / allsums_lads.s3524
allsums_lads['norm3'] = allsums_lads.s3494 / allsums_lads.s3518

for col_name in ['tri_13', 'tri_18', 'tri_16', 'norm1', 'norm2',  'norm3']: #  
    # thres = [0.78, 1.23]
    thres = np.nanquantile(allsums_lads[col_name].tolist(), [0.05, 0.85])
    
    for i in range(0,2):
        # DOWN
        if i == 0:
            df_bins = allsums_lads[allsums_lads[col_name] < thres[i]]
        # UP
        else:
            df_bins = allsums_lads[allsums_lads[col_name] > thres[i]]
        df_bins = df_bins.loc[:,['chrom', 'start', 'end', col_name]]
        # cluster selected bins
        df_binsClustered = bioframe.cluster(df_bins, min_dist=50000)
        # calculate mean value for all the clusters
        meanValue = df_binsClustered.groupby('cluster')[col_name].mean()
        # assign mean value to cluster
        df_bins = df_binsClustered[['chrom', 'cluster_start', 'cluster_end', 'cluster']].set_index('cluster').join(meanValue, on='cluster').drop_duplicates()
        df_bins.columns = ['chrom', 'start', 'end', 'value']
        # discard clusters consisting of 1 bin
        df_bins = df_bins.query('end - start > 50000')

        # calculate coverage by HK genes
        df_HK = bioframe.read_table('data/HK_genes.GeneCode.hg38.bed', schema='bed3')
        df_coverage = bioframe.coverage(df_bins, df_HK)
        df_coverage['HK_cov_fraction_file_%s' % (state_dic[col_name][j])] = df_coverage['coverage'] / (df_coverage['end'] - df_coverage['start'])
        
        # create dic with data to iterate through
        state_dic = {'norm1': ['s3492', 's3524'], 'norm2': ['s3492', 's3525'], 'norm3': ['s3494', 's3518'], 'tri_18': ['sHF18', 'sPFCH6'], 'tri_16': ['s3471', 's3494'], 'tri_13': ['s3475', 's3492']}

        for j in range(2):
            for comp in range(2):
                # compartemtns are encoded in column 4 as binary 
                # A - 1; B - 0
                df_compartments = bioframe.read_table('data/compartments/%s_100000_compartments.cis.bed' % state_dic[col_name][j], schema='bed6').query('name == @comp')
                # calculate coverage of clusters by compartments
                tmp_coverage = bioframe.coverage(df_bins, df_compartments)
                df_coverage['Comp%scov_fraction_file_%s' % (comp, state_dic[col_name][j])] = tmp_coverage['coverage'] / (tmp_coverage['end'] - tmp_coverage['start'])
        # save to file
        df_coverage.drop(['value', 'coverage'], axis=1).to_csv("./results/our_data_%s.%s.Clustered.withFracCovByCompart.tsv" % (col_name, i), sep="\t", index=False, header=True)


# %% public data
allsums_lads = pd.read_table('./data/1Mb_to_whole_cis_UPD.tsv')
#3492 vs 3524, 3492 vs 3525, 3494 vs 3518, IsoE vs IsoT, NPC_IsoE vs NPC_IsoT
allsums_lads['norm1'] = allsums_lads.s3492 / allsums_lads.s3524
allsums_lads['norm2'] = allsums_lads.s3492 / allsums_lads.s3525
allsums_lads['norm3'] = allsums_lads.s3494 / allsums_lads.s3518
allsums_lads['Itri_21'] = allsums_lads.IsoE / allsums_lads.IsoT
allsums_lads['tri_21'] = allsums_lads.NPC_IsoE / allsums_lads.NPC_IsoT


for col_name in ['norm1', 'norm2', 'norm3', 'Itri_21', 'tri_21']: 
    #thres = np.nanquantile(allsums_lads[col_name].tolist(), [0.05, 0.85])
    thres = [0.78, 1.23]
    
    for i in range(0,2):
        # DOWN
        if i == 0:
            df_bins = allsums_lads[allsums_lads[col_name] < thres[i]]
        # UP
        else:
            df_bins = allsums_lads[allsums_lads[col_name] > thres[i]]
        df_bins = df_bins.loc[:,['chrom', 'start', 'end', col_name]]
        # cluster selected bins
        df_binsClustered = bioframe.cluster(df_bins, min_dist=50000)
        # calculate mean value for all the clusters
        meanValue = df_binsClustered.groupby('cluster')[col_name].mean()
        # assign mean value to cluster
        df_bins = df_binsClustered[['chrom', 'cluster_start', 'cluster_end', 'cluster']].set_index('cluster').join(meanValue, on='cluster').drop_duplicates()
        df_bins.columns = ['chrom', 'start', 'end', 'value']
        # discard clusters consisting of 1 bin
        df_bins = df_bins.query('end - start > 50000')

        # HK genes occurence
        print(pValue_featureOccurenceInLoop(
                file_loops="./data/tmp1.bed",
                N_shuffle=1000,
                file_features='data/HK_genes.GeneCode.hg38.bed', 
                name='HKatClusters.%s.%s' % (col_name, i),
                genome_file='data/genome/hg38.chrom.sizes',
                pic_path='results/'))

        #protein-coding genes occurence
        print(pValue_featureOccurenceInLoop(
                file_loops="./data/tmp1.bed",
                N_shuffle=1000,
                file_features='./data/genome/gencode.v43.basic.annotation.genes.bed', 
                name='CodingGenesAtClusters.%s.%s.%s' % (col_name, i),
                genome_file='data/genome/hg38.chrom.sizes',
                pic_path='results/')) 
        