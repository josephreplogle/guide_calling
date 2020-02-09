import pandas as pd
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns; sns.set()
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from pomegranate import *

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
            
def load_gbc_reads(experiment):
    output_dir = os.path.join(experiment, 'outs')

    # output of cellranger count and cellranger aggr have slightly different directory structures
    if os.path.isdir(os.path.join(experiment, os.path.normpath("outs/filtered_gene_bc_matrices"))):
         matrix_dir = os.path.join(experiment, os.path.normpath("outs/filtered_gene_bc_matrices"))
    else:
         matrix_dir = os.path.join(experiment, os.path.normpath("outs/filtered_gene_bc_matrices_mex"))

    print "Loading observed cell barcodes..."

    barcodes_path = find('barcodes.tsv', matrix_dir)
    cell_barcodes = pd.read_csv(barcodes_path, sep='\t', header=None, names=['cell_barcode'])
    num_cells = len(cell_barcodes)

    print "Loading guides..."
    gbc_file = os.path.join(experiment, 'outs/guide_barcode_reads.txt.gz')

    gbc_reads = pd.read_table(gbc_file, header=None, names=('guide_identity', 'cell_barcode', 'UMI'))
    
    return gbc_reads, cell_barcodes

def capture_reads(gbc_reads, cell_barcodes):
    print "Filtering by cell barcode..."
    # we will be conservative and merge reads not just based on cell BC and UMI, but also based on mapping 
    # (as bowtie sometimes fails to assign a read to a particular guide, resulting in a *)
    gbc_umis = gbc_reads.drop_duplicates().groupby(['cell_barcode', 'guide_identity']).count().rename(columns={'UMI': 'UMI_count'}).reset_index()
    cell_barcode_umi_counts = gbc_umis[['cell_barcode', 'UMI_count']].groupby('cell_barcode').sum().sort_values('UMI_count', ascending=False)

    # gbc_reads = gbc_reads.query('guide_identity != "*"').copy()

    # count reads and UMIs for cell barcodes that appear in RNA-seq experiment
    # first reads, where we simply count the number of things with same cell_barcode and guide_identity call
    captured_gbc_reads = gbc_reads[gbc_reads['cell_barcode'].isin(cell_barcodes['cell_barcode'])]
    captured_gbc_read_counts = captured_gbc_reads.groupby(['cell_barcode', 'guide_identity']).count().rename(columns={'UMI': 'read_count'})

    # next UMIs... we will be conservative and merge reads not just based on cell BC and UMI, but also based on mapping 
    # (as bowtie sometimes fails to assign a read to a particular guide, resulting in a *)
    captured_gbc_UMI_counts = captured_gbc_reads.drop_duplicates().groupby(['cell_barcode', 'guide_identity']).count().rename(columns={'UMI': 'UMI_count'})

    # merge table
    captured_gbc_table = pd.merge(captured_gbc_read_counts, captured_gbc_UMI_counts, left_index=True, right_index=True)
    # look at (rough) coverage per UMI as a metric of quality of the call (note this is assuming reads distribute evenly over UMIs)
    captured_gbc_table['coverage'] = captured_gbc_table['read_count']/captured_gbc_table['UMI_count']

    captured_gbc_table['gemgroup'] = pd.Series(captured_gbc_table.index.get_level_values(0).map(lambda x: x.split('-')[1]).astype(int),                                                     index=captured_gbc_table.index)
    return captured_gbc_table

def find_threshold(coverage_data, gemgroup):
    # find location of upper mode by fitting KDE and finding extrema
    model = gaussian_kde(coverage_data)
    x = np.linspace(0, coverage_data.quantile(0.99), 200)
    extrema = argrelextrema(model(x), np.greater)[0]
    upper_mode = x[extrema[-1]]
    # now find local minimum between this and the rest of the data
    x_low = x[0:extrema[-1] - 1]
    low_extrema = argrelextrema(model(x), np.less)[0][-1]
    cutoff = x[low_extrema]

    ax = sns.distplot(coverage_data, kde_kws={"bw":0.5})
    yrange = ax.get_ylim()
    plt.plot([cutoff, cutoff], [0, yrange[1]], 'r')
    plt.plot([upper_mode, upper_mode], [0, yrange[1]], 'g')
    sns.despine()
    plt.title('Gemgroup {0}'.format(gemgroup))
    
    print 'Gemgroup {0}: threshold is {1:0.2f}'.format(gemgroup, cutoff)
    print 'Gemgroup {0}: average coverage of high quality identities is {1:0.2f}'.format(gemgroup, coverage_data[coverage_data > cutoff].mean())
    
    return cutoff

def identify_cells(captured_gbc_table, coverage_thresholds, read_threshold=50, umi_threshold=3):
    # for each cell we will take the top identity by read count (UMI count is finicky because then a bunch of *'s will infiltrate at low coverage)
    cell_identities = captured_gbc_table.reset_index().sort_values('read_count', ascending=False).groupby('cell_barcode').head(1)

    captured_gbc_table['good_coverage'] = False

    for gemgroup in captured_gbc_table['gemgroup'].unique():
        entries = cell_identities['gemgroup'] == gemgroup
        best_coverage_mean = cell_identities[entries]['coverage'].mean()
        best_coverage_std = cell_identities[entries]['coverage'].std()

        entries = captured_gbc_table['gemgroup'] == gemgroup
        captured_gbc_table.loc[entries, 'good_coverage'] = (captured_gbc_table[entries]['coverage'] > coverage_thresholds[gemgroup]) &                                                (captured_gbc_table[entries]['read_count'] >= read_threshold) &                                                (captured_gbc_table[entries]['UMI_count'] >= umi_threshold)

    cell_identities = captured_gbc_table.reset_index().sort_values(['good_coverage', 'read_count'], ascending=[False, False])                                              .groupby('cell_barcode').head(1).set_index('cell_barcode')
    cell_identities['number_of_cells'] = captured_gbc_table.reset_index().groupby('cell_barcode').apply(lambda x: x['good_coverage'].sum())

    return cell_identities

def write_identities(experiment, cell_identities, cell_barcodes):
    output_dir = os.path.join(experiment, 'outs')
    
    cell_identities_summary = cell_identities[cell_identities['number_of_cells']==1].groupby('guide_identity').count()[['number_of_cells']].sort_values('number_of_cells', ascending=False)
    dummy = range(1, len(cell_identities_summary)) + [' ',]*5

    num_cells = len(cell_barcodes)
    num_identified = cell_identities_summary[~cell_identities_summary.index.isin({'*'})]['number_of_cells'].sum()
    num_multiplets = cell_identities[cell_identities['number_of_cells']>1].groupby('guide_identity').count()['number_of_cells'].sum()
    num_unidentified = num_cells - num_identified - num_multiplets

    cell_identities_summary_stats = pd.DataFrame([num_multiplets, num_identified, num_unidentified, num_cells],                                                  index=['Multiplets', 'Total uniquely identified', 'Total unidentifiable', 'Total number of cells'], columns=['number_of_cells'])
    cell_identities_summary = cell_identities_summary.append(cell_identities_summary_stats)
    guide_percentages = cell_identities_summary[['number_of_cells']].rename(columns={'number_of_cells': 'percentage'})/num_cells*100
    cell_identities_summary = pd.merge(cell_identities_summary,guide_percentages,left_index=True, right_index=True).round(2)
    cell_identities_summary = pd.merge(pd.DataFrame(dummy, columns=['rank']), cell_identities_summary.reset_index(), left_index=True, right_index=True).rename(columns={'index': 'guide_identity'}).set_index(['rank', 'guide_identity'])

    cell_identities_summary.to_csv(os.path.join(output_dir, 'cell_identities_summary.csv'))
    cell_identities_summary.to_html(os.path.join(output_dir, 'cell_identities_summary.html'))
    filename = os.path.join(output_dir, 'cell_identities.csv')
    cell_identities.to_csv(filename)
    
    print "Identities written to " + filename
    return cell_identities_summary

def plot_umi_distribution(table, gemgroup):
    ax = sns.distplot(table.query('good_coverage')['UMI_count'])
    yrange = ax.get_ylim()
    umi_mean = table.query('good_coverage')['UMI_count'].mean()
    plt.plot(np.array([1,1])*(umi_mean), [0, yrange[1]])
    plt.title('Gemgroup {0}'.format(gemgroup))
    print 'Gemgroup {0}: UMI mean is {1:0.2f}'.format(gemgroup, umi_mean)
    sns.despine()
    
def plot_umi_read_distribution(table, gemgroup):
    num_identities = table.groupby(level=0).count()['read_count']

    single_identities = num_identities[num_identities >= 1]
    double_identities = num_identities[num_identities >= 2]
    most_likely_identities = table.loc[single_identities.index.tolist()].sort_values('read_count', ascending=False).reset_index().groupby('cell_barcode').apply(lambda x: x.iloc[0])
    second_most_likely_identities = table.loc[double_identities.index.tolist()].sort_values('read_count', ascending=False).reset_index().groupby('cell_barcode').apply(lambda x: x.iloc[1])
    confident_most_likely_identities = most_likely_identities.query('good_coverage')
    nonconfident_most_likely_identities = most_likely_identities.query('~good_coverage')
    confident_second_most_likely_identities = second_most_likely_identities.query('good_coverage')
    
    plt.scatter(np.log2(most_likely_identities['read_count']),
                np.log2(most_likely_identities['UMI_count']), s=3, alpha=0.5, color='blue', label='Single identity')
    plt.scatter(np.log2(nonconfident_most_likely_identities['read_count']),
                np.log2(nonconfident_most_likely_identities['UMI_count']), s=3, color='black', label='Nonconfident single identity')
    plt.scatter(np.log2(second_most_likely_identities['read_count']),
                np.log2(second_most_likely_identities['UMI_count']), s=3, alpha=0.5, color='red', label='Rejected second identity')
    plt.scatter(np.log2(confident_second_most_likely_identities['read_count']), 
                np.log2(confident_second_most_likely_identities['UMI_count']), s=3, alpha=0.5, color='green', label='Multiplet')

    plt.xlabel('log_2 guide barcode read count')
    plt.ylabel('log_2 guide barcode UMI count')
    sns.despine()
    plt.legend(loc='upper left')

    
def MixedModelCall(guide,gbc_table,library,directory):
    data=np.array(np.log2(gbc_table.reset_index()[gbc_table.reset_index()['guide_identity']==guide]['UMI_count']))
    data=data.reshape(-1,1)
    
    ##this loop checks that the model converges 
    ##and that the Poisson distribution is the lower component 
    i=0
    gmm_x = np.linspace(-2,max(data)+2,1000)
    while i==0:
        model = GeneralMixtureModel.from_samples([PoissonDistribution,NormalDistribution],2,data)
        if numpy.isnan(model.probability(gmm_x)).any():
            i=0
        else:
            if model.distributions[0].name is 'PoissonDistribution':
                if model.distributions[0].parameters[0]<model.distributions[1].parameters[0]:
                    i=1    
                else:
                    i=0
            elif model.distributions[0].parameters[0]>model.distributions[1].parameters[0]:
                    i=0
    
    ##plot
    MixedModelPlot(guide,gbc_table,library,directory,model,data,gmm_x)

    ##append positive calls to table
    return gbc_table.reset_index()[gbc_table.reset_index()['guide_identity']==guide].loc[model.predict_proba(data)[:,1]>0.5]

def MixedModelPlot(guide,gbc_table,library,directory,model,data,gmm_x):
    ## Plot histograms and gaussian curves
    fig, ax = plt.subplots()
    sns.distplot(data,kde=False,norm_hist=True,bins=22,color='grey')
    ax.plot(gmm_x, model.distributions[0].probability(gmm_x)*np.exp(model.weights[0]),color='red',linestyle=':')
    ax.plot(gmm_x, model.distributions[1].probability(gmm_x)*np.exp(model.weights[1]),color='black',linestyle=':')
    ax.set(xlabel='log2 '+guide+' UMI per cell', ylabel='Frequency')
    sns.despine()
    plt.tight_layout()
    plt.savefig(directory+library+'_'+guide+'_MixedModel.pdf',bbox_inches='tight')
    plt.close()
