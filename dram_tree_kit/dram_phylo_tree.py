"""This program add profiles based on phylogentic trees into dram. The key to this process is pplacer which places leaves into pre-exitsting trees"""
#! /opt/Miniconda2/miniconda2/envs/DRAM2/bin/python
import os
from io import StringIO       
import click
import logging
import pandas as pd
import tempfile
# import matplotlib.pyplot as plt
from dram_tree_kit import __version__

from tempfile import TemporaryDirectory
from mag_annotator.utils import run_process
from mag_annotator.utils import setup_logger
from mag_annotator.summarize_genomes import get_ids_from_annotations_by_row
from Bio import Phylo as phy
from dram_tree_kit.pplacer import DramTree
from mag_annotator.pull_sequences import pull_sequences
from skbio import write as write_sq
from skbio import read as read_sq
import networkx as nx
from networkx.drawing.nx_agraph import write_dot
# Tree = namedtuple('Tree', ('name', 'tree', 'alignment_hmm', 'mapping'))

#     Tree('nxr_nar', 
#          phy.read('/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metabolism/nxr_nar/bestTree.nxr-nar_seqs_for_tree_aligned.faa_mode_low.renamed', "newick"), 
#          './second_try/nxr-nar_seqs_for_tree_aligned.hmm',
#          )]
"""
os.system('dram_tree_kit//dram_phylo_tree.py -a ./example_one/all_bins_combined_3217db_ACTIVE_GENES_annotations.txt -g ./example_one/all_bins_combined_3217db_genes.faa -c 30')
os.system('ls ')

"""

NXR_NAR_TREE = DramTree(pplacer_profile='./data/dram_trees/nxr_nar/nxr_nar.refpkg', 
                target_ids = ['K11180', 'dsrA', 'dsrB', 'K11181'],
                reference_seq=('./data/dram_trees/nxr_nar/'
                               'nxr-nar_seqs_for_tree_aligned.faa'),
                tree_mapping_path='data/dram_trees/nxr_nar/nxr-nar-tree-mapping.tsv')
TREES = [NXR_NAR_TREE]

def color_tree(tree):
    pass


def extract_enigmatic_genes(annotation_ids:pd.DataFrame, gene_fasta, work_dir, 
                            target_ids:set, logger) -> (pd.DataFrame, str):
    """
    :param annotation_ids: annotations from dram run
    :param gene_fasta: faa from dram run
    :param work_dir: Temp files here
    :param target_ids: ID set needing phylo info, used to filter genes
    :param logger: Standard DRAM logger
    :returns: The path to the ambiguous genes in a fasta file

    Takes in a fasta file of genes and a list of ids in the dram annotation, and returns a filtered fasta to match.
    """
    output_fasta = os.path.join(work_dir, 'trim.faa')
    ids_keep = ids[annotation_ids.apply(lambda x: len(x.intersection(target_ids))>0)]
    output_fasta_generator = (i for i in read_sq(gene_fasta, format='fasta') 
                              if i.metadata['id'] in ids_keep.index)
    write_sq(output_fasta_generator, format='fasta', into=output_fasta)
    return output_fasta

def test_extract_enigmatic_genes(temp_dir):
    temp_dir='../temp' 
    logger = logging.getLogger()
    annotations = pd.DataFrame({'ko_id': ['K00001', 'K00002']},
                               index=['one', 'two'])
    annotation_ids = get_ids_from_annotations_by_row(annotations, logger)
    start_fa = os.path.join(temp_dir, 'in.faa')
    end_fa = os.path.join(temp_dir, 'trim.faa')
    with open(start_fa, 'w') as out:
        out.write('>one\naabb\n>two\nbbaa')
    extract_enigmatic_genes(annotations, start_fa, temp_dir, {'K00001'}, 
                            logger)
    with open(end_fa, 'r') as file:
        result = file.read()
    assert result == '>one\naabb\n'

def make_adjective_apbundace(gene_adj, adjective_abundance):
    adj_abu = gene_adj.groupby('adjective').apply(lambda x: pd.concat([x.drop('adjective', axis=1).sum(), pd.Series({"genes": len(x.index.unique())})]))
    adj_abu.to_csv(adjective_abundance, sep='\t',)


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
@click.option('-o', 'annotations_out', type=click.Path(exists=False), required=False, help="The output annotations file")
@click.option('-c', 'cores', type=int, required=False, help="The number of cores to use", default=10)
@click.option('--annotate_all', is_flag=True, show_default=True, default=False, 
              help="Don't place uncertain genes, place all of them, this will involve splitting and multiprocessing")
def dram_tree_kit(
    dram_annotations:str=None,
    gene_fasta:str=None,
    dram_directory:str=None,
    annotations_out:str=None,
    annotate_all:bool=False,
    cores:int=10
    ):
    logger = logging.getLogger('dram_tree_log')
    tree = NXR_NAR_TREE
    tree.set_logger(logger)
    annotations = pd.read_csv(dram_annotations, sep='\t', index_col=0)
    annotation_ids = get_ids_from_annotations_by_row(annotations, logger)
    with TemporaryDirectory(dir=output_dir) as work_dir:
        trimed_fa = extract_enigmatic_genes(annotation_ids, gene_fasta, work_dir, 
                                tree.target_ids, logger)
        jplace_file = tree.pplacer_place_sequences(trimed_fa, work_dir, threads=cores)
        breakpoint()

#gup
#
#output_dir = './'
#query_faa = "./../second_try/query.faa"
#query_faa = '../first_try/example_one/all_bins_combined_3217db_genes.faa'
dram_tr


# pplacer/guppy columns in one place for reference
PLACEMAT_NAME_1_COL = 'name1'
PLACEMAT_NAME_2_COL = 'name2'
PLACEMAT_DISTANCE_COL = 'distance'


def distance_to_call(row):
    place_seq:str = None
    call_id:str = None
    for name in [row[name1], row[name2]]:
        if name in annotations.index and place_seq is None:
            place_seq = name
            continue
        if name in tree.mapping.index and call_id is None:
            call_id = tree.mapping.loc[name]
            continue
    return  pd.Series({'seq': place_seq, 'col': call_id,
                       'distance': row[dist]})


tree = NXR_NAR_TREE
logger = logging.getLogger('dram_tree_log')
annotations = pd.read_csv('example_one/all_bins_combined_3217db_ACTIVE_GENES_annotations.txt', sep='\t', index_col=0)
jplace_file = 'example_place_output.jplace'
def lable_tree_distance(tree, jplace:str, annotations:pd.DataFrame,
                         logger:logging.Logger):
    placemat_txt = run_process(['guppy', 'distmat', jplace_file], 
                               logger, capture_stdout=True)
    name1 = PLACEMAT_NAME_1_COL
    name2 = PLACEMAT_NAME_2_COL
    dist = PLACEMAT_DISTANCE_COL
# run_process(['guppy', 'classify', jplace_file], logger)
# collect data
# distances = pd.read_csv(StringIO(placemat_txt), sep='\t')
# place_csv_txt = run_process(['guppy', 'to_csv', '--no-csv', jplace_file,], logger, capture_stdout=True)
# pd.read_csv(StringIO(place_csv_txt), delim_whitespace=True).iloc[100]

def lable_tree_placement(tree, jplace:str, annotations:pd.DataFrame,
                         logger:logging.Logger):
    pass
place_tree_txt = run_process(['guppy', 'tog', '-o', 'jplace_to_tree', jplace_file,], logger, capture_stdout=True)
' '.join(['guppy', 'tog', '-o', 'jplace_to_tree', jplace_file,])

tree = NXR_NAR_TREE
lable_tree_placement(tree, jplace, annotations, logger)
treeph = phy.read('jplace_to_tree', format="newick")
import graphviz
phy.draw_ascii(treeph)
color_map = {'nxr/nar-N utilization': 'purple', 
             'other-None': 'blue', 
             'narG-N reducer': 'green', 
             'nxr-Nitrifier': 'red',
             'nxr-None':  'orange'}
tree.mapping['color'] = tree.mapping['call'].apply(lambda x: color_map[x])
# tree.mapping['call'].unique()
for cl in treeph.get_terminals():
    cl.color = None if cl.name not in tree.mapping.index else tree.mapping.loc[cl.name, 'color']

# lets do this

phy.write(treeph, 'color_tree.xml', 'phyloxml')

paths = [treeph.get_path(cl) for cl in treeph.get_terminals() if cl.name in tree.mapping.index]

def split_out(pat):
    if len(pat) == 1:
        i = pat[0]
        return [(tree.mapping.loc[i[-1].name, 'call'], set(i))]
    heads = {i[0] for i in pat}
    if 1 == len(heads):
        if 1 == len(set(calls:=[tree.mapping.loc[i[-1].name, 'call']  for i in pat])):
            return [(calls[0], {i for j in pat for i in j})]
    new_pats =  [[i if len(i) ==1 else i[1:] for i in pat if i[0] == j ] for j in heads]
    iteration = [j for i in new_pats for j in split_out(i) ]
    return iteration

confident_paths = split_out(paths)
for i in confident_paths:
    color = color_map[i[0]]
    for cl in i[1]:
        cl.color = color

phy.write(treeph, 'color_tree_branch.xml', 'phyloxml')
for cl in treeph.get_terminals():
    add_color = None
    for i in treeph.get_path(cl):
        print(i.color, i.color is not None) 
        if i.color is not None and i.color != (0, 0, 0):
            add_color = i.color
        i.color = add_color
    cl = add_color
    

phy.write(treeph, 'color_tree_all.xml', 'newick')

if __name__ == '__main__':
    dram_tree_kit()

import json
with open(jplace_file) as jp_to_load:
    place_dict = json.load(jp_to_load)


phy.draw_ascii(treeph)
