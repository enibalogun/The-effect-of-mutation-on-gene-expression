"""
<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Assumptions-for-ANOVA-Test---DOES-NOT-FULFILL!!" data-toc-modified-id="Assumptions-for-ANOVA-Test---DOES-NOT-FULFILL!!-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Assumptions for ANOVA Test - DOES NOT FULFILL!!</a></span></li></ul></div>
"""

# %matplotlib inline

"""
<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Assumptions-for-ANOVA-Test---DOES-NOT-FULFILL!!" data-toc-modified-id="Assumptions-for-ANOVA-Test---DOES-NOT-FULFILL!!-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Assumptions for ANOVA Test - DOES NOT FULFILL!!</a></span></li></ul></div>
"""

# %matplotlib inline

""
# %matplotlib inline
import matplotlib.pyplot as plt
from Bio import SeqIO
import random, glob
import json
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec
import seaborn as sns
import scipy
from scipy import stats
from scipy.stats import mannwhitneyu
stats.junk = lambda chisq, df: stats.chi2.sf(chisq, df)
import csv
import gffpandas.gffpandas as gffpd
import statsmodels.formula.api as smf
from statsmodels.formula.api import ols

CC2344 = ["CC2344-L1", "CC2344-L2", "CC2344-L3", "CC2344-L4", "CC2344-L5", "CC2344-L6", "CC2344-L7", "CC2344-L8", "CC2344-L9", "CC2344-L10", "CC2344-L11", "CC2344-L12", "CC2344-L13", "CC2344-L14", "CC2344-L15"]
CC2931 = ["CC2931-L1", "CC2931-L2", "CC2931-L3", "CC2931-L4", "CC2931-L5", "CC2931-L6", "CC2931-L7", "CC2931-L9", "CC2931-L10", "CC2931-L11", "CC2931-L13", "CC2931-L14", "CC2931-L15"]

#### Dataframe recording the generation time per sample ####
dic_gen = {}
generation = pd.read_csv('/research/projects/chlamydomonas/MAexpression/genome_info/mutation_info/Mutation_Fitness.txt', delimiter = '\t')
generation['Sample'] = generation['Sample'].str.replace('_', '-L')
generation = generation.loc[generation['Sample'].isin(CC2344 + CC2931)]

for i in generation.index.values:
    dic_gen[generation.at[i,'Sample']] = generation.at[i, 'Generation']
generations = pd.Series(dic_gen)

#### Dataframe recording the number of mutations per sample ####
mutations = pd.read_csv('/research/projects/chlamydomonas/MAexpression/genome_info/mutation_info/all_mutations.csv', delimiter = '\t')
dic_mut = {maline:mutations.loc[mutations['sample'] == maline].count()[0] for maline in mutations['sample'].values.tolist()}

##########################################################
# #### Assumptions for ANOVA Test - DOES NOT FULFILL!! ####
# #########################################################
# (1) Normal distribution
# (2) Homogeneity of variance

import re
from collections import OrderedDict


""
class GFF_line:
	"""This class basically parses a GFF line and allows you to interact with different components that I have deemed interesting
	Most components are simple strings or intgers.
	The attributes field which is a ;-separated list is returned as a dictionary """
	def __init__(self, l, info_delimiter=";", info_field_delimiter = '='):
		self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attribs = l.split('\t')
		self.attributes = self.attribute_dict(self.attribs, info_delimiter, info_field_delimiter)
		self.start = int(self.start)
		self.end = int(self.end)
		self.line = l
	def attribute_dict(self, attributes, info_delimiter=";", info_field_delimiter = '='): ###this is fragile
		d = OrderedDict()
		attributes = attributes.strip(info_delimiter)
		for i in [x.strip() for x in re.split(info_delimiter, attributes)]:
			if len(i)>0:
				field = re.split(info_field_delimiter, i)[0]
				d[field] = re.split(info_field_delimiter, i)[1].strip('"')
		return d
	def retrieve_sequence(self, ref_dict, reverse_complement=False):
		seq = str(ref_dict[self.seqid].seq[self.start-1: self.end])
		if reverse_complement:
			seq = reverse_complement(seq)
		return seq

def randomDEGs(genes, nDEGs):
    random_DEGs = random.choices(list(genes.keys()), k=nDEGs)
    return random_DEGs

def nearest_mutation(DEG, genes, gene_bounds, muts):
    """
    #set DEGs with mutations in them as 0
    #otherwise use the minimum absolute distance to the genes flank

    DEG= gene name
    genes = a look up table of where genes are
    gene_bounds = a look up table of the chromosome position where a gene is
    muts = a dictionary of mutations keyed by chromosome
    """
    if genes[DEG] == 'contig':
        return None
    else:
        c = genes[DEG]
        bounds = gene_bounds[c][DEG]
        min_distance = 1e9 #set to some huge number bigger than the whole genome
        if c not in muts:
            return None
        for p in muts[c]: # for this DEG we will check its distance to every mutation
            if bounds[0]< p <bounds[1]: #if its in the gene
                min_distance = 1 # set to one - you can use 0 but it messes with log scales
                return min_distance #you're done
            current_distance = min(abs(p-bounds[0]),abs(p-bounds[1])) #the nearest distance of this DEG to this mutation is the min of its flanks absolute distance 
            if current_distance < min_distance: # if this current min is closer than the running min then set it has the new running min
                min_distance = current_distance
    return min_distance


""
GFF_file  = "/research/references/chlamydomonas/6.0_chlamy/CC4532.v1_1.gene_exons.gff3"

gene_types = ['gene', 'transposable_element_gene', 'tRNA', 'group_II_intron', 'rRNA', 'group_I_intron','ncRNA']
genes ={}
gene_bounds = {}

for l in open(GFF_file):
    if l.startswith("#"):continue
    g = GFF_line(l)
    if 'contig' in g.seqid: 
        ID = g.attributes['ID']
        genes[ID] = 'contig'
    if g.type in gene_types:
        ID = g.attributes['ID']
        if g.seqid in gene_bounds: gene_bounds[g.seqid][ID] = [g.start, g.end]
        else:
            gene_bounds[g.seqid] = {ID:[g.start, g.end]}
        genes[ID] = g.seqid

""
import glob
DEG_folder = "/research/projects/chlamydomonas/MAexpression/analysis/genes_log2fold"

deg_files = glob.glob(DEG_folder+"/CC*-L*") 
obs_DEGs = {}
for f in deg_files:
    ma_line = f.split("/")[-1]
    obs_DEGs[ma_line] = []
    #print(ma_line)
    for l in open(f).readlines()[1:]:
        gene, baseMean,log2FoldChange,lfcSE,stat,pvalue,padj = l.strip().split("\t")
        try:
            padj = float(padj)
            if padj < 0.05:
                obs_DEGs[ma_line].append(gene)
                #could apply other filters here - ie strong down regulation
#                 if abs(float(log2FoldChange)) > 2.0:
#                     obs_DEGs[ma_line].append(gene)
        except:
            pass
            #print(l.strip())

            
for i in obs_DEGs:
    print(i, len(obs_DEGs[i]))

""
mutations_file  = "/research/projects/chlamydomonas/MAexpression/genome_info/mutation_info/v6_coordinates.bed"

obs_mutations = {}
for l in open(mutations_file):
    c,s,e,ma_line =l.strip().split()
    ma_line = "-L".join(ma_line.split("_"))
    if ma_line not in obs_mutations:obs_mutations[ma_line] = {c:[int(e)]}
    elif c not in obs_mutations[ma_line]:obs_mutations[ma_line][c] = [int(e)]
    else:obs_mutations[ma_line][c].append(int(e))
    
for ma_line in obs_mutations:
    print(ma_line, sum([len(obs_mutations[ma_line][c]) for c in obs_mutations[ma_line]]))

""
CC2931 = ["CC2931-L1", "CC2931-L2", "CC2931-L3", "CC2931-L4", "CC2931-L5", "CC2931-L6", "CC2931-L7", "CC2931-L9", "CC2931-L10", "CC2931-L11", "CC2931-L13", "CC2931-L14", "CC2931-L15"]
CC2344 = ["CC2344-L1", "CC2344-L2", "CC2344-L3", "CC2344-L4", "CC2344-L5", "CC2344-L6", "CC2344-L7", "CC2344-L8", "CC2344-L9", "CC2344-L10", "CC2344-L11", "CC2344-L12", "CC2344-L13", "CC2344-L14", "CC2344-L15"]

p_value = pd.DataFrame()

for f in [100]: 
    n_sims = 10000
    bins = [0] + [10**i for i in range(1,8)]
    flank_region = f 
    total_sim_ds = []
    total_obs_ds = []
    significant= 0
    for ma_line in CC2931 + CC2344:
        DEGs = obs_DEGs[ma_line]
        muts = obs_mutations[ma_line]
        obs_distances = [nearest_mutation(d, genes, gene_bounds, muts) for d in DEGs]
        fraction_in_flank = 0
        for i in obs_distances:
            if i != None:
                if i < flank_region:fraction_in_flank+=1
        n=0
        sim_distances = []
        closer_than_obs = 0
        while n < n_sims:
            simDEGs = randomDEGs(genes, len(DEGs))
            these_sim_distances = [nearest_mutation(d, genes, gene_bounds, muts) for d in simDEGs]
            sim_distances = sim_distances+these_sim_distances
            sim_fraction_in_flank = 0
            for i in these_sim_distances:
                if i != None:
                    if i < flank_region:
                        sim_fraction_in_flank+=1
            if sim_fraction_in_flank > fraction_in_flank:closer_than_obs+=1
            n+=1
        title_line = "{} DEGs: {}, muts: {}, p-value={}".format(ma_line, len(DEGs), len(muts), closer_than_obs/n)
        p_value.at[ma_line, 'pval'] = closer_than_obs/n

        if closer_than_obs/n <=0.05: 
            significant +=1
            print(significant)
#         plt.figure(figsize = (20,15))
#         sns.histplot(sim_distances, stat = 'density', log_scale=True, bins=25)
#         sns.histplot(obs_distances, stat = 'density', log_scale=True, color='red', bins=25, alpha=.4)
#         plt.xlabel('Distance from nearest mutation (bps)', fontsize = 35)
#         plt.ylabel('Density', fontsize = 35)
#         plt.xticks(fontsize = 25)
#         plt.yticks(fontsize = 25)
#         plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/DEGs/cis_mutations/CC2344-L14_distribution.eps', dpi = 600, format = 'eps', bbox_inches = 'tight')
#         plt.savefig('/research/projects/chlamydomonas/MAexpression/analysis/DEGs/cis_mutations/CC2344-L14_distribution.pdf', dpi = 600, format = 'pdf', bbox_inches = 'tight')
        
#         total_sim_ds +=sim_distances
#         total_obs_ds +=obs_distances
p_value.to_csv('/research/projects/chlamydomonas/MAexpression/analysis/DEGs/cis_mutations/p_values.csv', sep = '\t', index = True, header = True)
