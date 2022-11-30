"""Setup file for package"""
from setuptools import setup, find_packages
from os import path

__version__="b0.0.1"
m_author__ = 'rmflynn'

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="dram2_tree_kit",
    version=__version__,
    packages=['dram2.tree_kit'],
    entry_points={
        'console_scripts': [
            "dram2-tree = dram_tree_kit.dram_phylo_pipe:dram_tree_kit"
        ]
    },
    data_files=[],
    zip_safe=False,
    license='GPL3',
    description="Provides Phylogenetic tree disambiguation for genes in the Distilled and Refined Annotation of Metabolism (DRAM) toolkit",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    python_requires='>=3',
    install_requires=['pandas', 'dram2_bio', 'click', 'biopython'],
    author="Rory Flynn",
    author_email='rory.flynn@colostate.edu',
    url="https://github.com/rmflynn/DRAM_tree_kit/",
    download_url="https://github.com/wrightonlab/DRAM_tree_kit/tarball/%s" % __version__
)
