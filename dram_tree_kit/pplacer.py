
import os
import json
from tempfile import TemporaryDirectory
from mag_annotator.utils import run_process
from mag_annotator.utils import setup_logger

import logging
QUERY_TMP_NAME = 'query'

class DramTree:

    def __init__(self, pplacer_profile:str, reference_seq:str, 
                 logger:logging.Logger):
        self.pplacer_profile = pplacer_profile
        self.reference_seq = reference_seq 
        self.logger = logger

    def _filter_and_align(self, fasta, working_dir):

        logging.info("Aligning sequences with MAFFT")
        # The file at the end of the output_path have extention 'fasta'!
        output_path = os.path.join(working_dir, f"{QUERY_TMP_NAME}.fasta")
        if os.path.exists(output_path):
            logger.warn(f"The alignment file {output_path} already "
                        "exists, so the program is being cald twice in "
                        "the same working directory.")
        run_process( ['mafft', '--add', fasta, '--reorder', self.reference_seq], 
                    self.logger,  capture_stdout=False, save_output=output_path)

        logging.info("Sequences aligned")
        return output_path


    def pplacer_place_sequences(self, query_faa:str, working_dir:str,
                                threads:int):
        logging.info("Sequence placement with pplacer started")
        output_path = os.path.join(working_dir, f"{QUERY_TMP_NAME}.jplace")
        run_process(['pplacer', '-c', self.pplacer_profile, 
                     self._filter_and_align(query_faa, working_dir),
                     '--out-dir', working_dir,
                      '-j', str(threads)],
                    self.logger)
        logging.info("Sequence placement with pplacer finished")
        return json.load(open(output_path))








"""
One day we may want to be able to make these profiles

taxit_comand= ['taxit', 'create', '-l', '16s_rRNA', '-P', 'nxr_nar.refpkg', 
 '--aln-fasta', '$NXR_NAR_FAA', '--tree-stats', 
 './RAxML_info.nxr_nar_raxml', '--tree-file', './nxr_nar.tre']
"""

