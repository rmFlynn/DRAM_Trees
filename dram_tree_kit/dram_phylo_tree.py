#! /opt/Miniconda2/miniconda2/envs/DRAM2/bin/python
"""If you want to use trees in dram, here is your opertunity"""
import os
import click
import pandas as pd
import tempfile
import matplotlib.pyplot as plt
from Bio import Phylo as phy
from collections import namedtuple
from mag_annotator.pull_sequences import pull_sequences
Tree = namedtuple('Tree', ('name', 'tree', 'alignment_hmm', 'mapping'))

TREES = [
    Tree('nxr_nar', 
         phy.read('/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metabolism/nxr_nar/bestTree.nxr-nar_seqs_for_tree_aligned.faa_mode_low.renamed', "newick"), 
         './second_try/nxr-nar_seqs_for_tree_aligned.hmm',
         pd.read_csv('./second_try/nxr-nar-tree-mapping.tsv', sep='\t', index_col=0)
         )]


def color_tree(tree):
    pass


tree = TREES[0].tree
mapping = TREES[0].mapping
terminals = {i.name for i in tree.get_terminals()}
mapping = mapping.loc[list(terminals.intersection(set(mapping.index)))]
mapping.loc['gunnisonriver_2019_sw_WHONDRS-S19S_0062_A_bin.22_Ga0451722_0001801_3']
mapping['Call based on tree'].unique()
color_map = {'nxr/nar': 'purple', 'other': 'black', 'narG': 'green', 'nxr': 'red'}
mapping['color'] = mapping['Call based on tree'].apply(lambda x: color_map[x])

for cl in tree.get_terminals():
    cl.color = 'black' if cl.name not in mapping.index else mapping.loc[cl.name, 'color']

plt.clf()
fig = plt.subfigure(figsize=(50, 10))
fig = plt.figure(figsize=(50, 10))
axes = fig.add_axes([0.5, 1, 0.5, 1])

phy.draw(tree, rcParams["figure.figsize"]=(20,3))
phy.draw_ascii(tree)
plt.savefig('tree_nodes.pdf', dpi=600, )

import networkx as nx
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
treenx = phy.to_networkx(tree)
pos = graphviz_layout(treenx, prog="dot")
nx.draw(T, pos)
networkx.graphviz_layout = nx_agraph.graphviz_layout
phy.draw_ascii(tree, prog='dot')

tree.tree
tree.name
py

tree.get_path(tree.get_terminals()[1].name)[-1]

# tree = phy.read("./bipartitionsBranchLabels.dsr_genes_for_tree_aligned.faa_mode_low.renamed", "newick")
ls $NXR_NAR
TREE=$NXR_NAR
print(tree)
tree.get_path(tree.get_terminals()[1].name)[-1]
__version__ = '0.0.1'


def extract_enigmatic_genes(annotations:pd.DataFrame, gene_fasta) -> (pd.DataFrame, str):
    return annotations_enigma, gene_enigma

def extract_enigmatic_genes(gene_enigma):
    return genes_alined
mapping

def make_adjective_apbundace(gene_adj, adjective_abundance):
    adj_abu = gene_adj.groupby('adjective').apply(lambda x: pd.concat([x.drop('adjective', axis=1).sum(), pd.Series({"genes": len(x.index.unique())})]))
    adj_abu.to_csv(adjective_abundance, sep='\t')


@click.command()
@click.version_option(__version__)
@click.option('-a', 'dram_annotations', type=click.Path(exists=True), required=False,
                help="The DRAM annotations file, not necessary"
                "  if you use the dram_directory option")
@click.option('-g', 'gene_fasta', type=click.Path(exists=True), required=False,
                help="The gene fasta file, genes.faa file from dram output.")
@click.option('-d', 'dram_directory', type=click.Path(exists=True), required=False,
                help="The dram input file, with no names changed so it contains annotations.txt and genes.faa, genes.fna")
@click.option('-o', 'annotations_out', type=click.Path(exists=False), required=False,
                help="The gene_abundances with out adjectives,"
                " optional if you have the gene_adjectives file")
@click.option('-o', 'annotations_out', type=click.Path(exists=False), required=False,
def adjective_abundances(
    dram_annotations:str=None,
    gene_fasta:str=None,
    dram_directory:str=None,
    annotations_out:str=None
    ):
        
        
    if adjectives_strainer is not None and \
            input_gene_abundance is not None and \
            (gene_adjectives is None or not os.path.exists(gene_adjectives)):
        gene_adj_data = make_gene_adjectives(adjectives_strainer, input_gene_abundance, gene_adjectives)
    elif gene_adjectives is not None and os.path.exists(gene_adjectives):
        gene_adj_data = pd.read_csv(gene_adjectives, sep='\t', index_col=0) 
    else:
        ValueError("Something went wrong, check your arguments")
        

    if gene_adj_data is not None and \
            adjective_abundance is not None:
        make_adjective_apbundace(gene_adj_data, adjective_abundance)

if __name__ == '__main__':
    dram_placer()
