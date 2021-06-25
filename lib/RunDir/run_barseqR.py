#!python3

import os
from datetime import datetime
import logging
from RunDir.BarSeqR import prepare_all_poolcount_etc 
from BarSeqPy.RunFEBA import RunFEBA 

# Following Function primarily used to set location of config for running BarSeq
def RunBarSeq(arg_list):
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    cfg_fp = os.path.join(this_file_dir, "barseqr_config_dict.json")

    # Getting time for log before program
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    logging.info("Beginning to run BarSeqR {}:-----------".format(dt_string))

    # Run Preparation Here
    prep_vars = prepare_all_poolcount_etc(cfg_fp, arg_list, this_file_dir)

    RunFEBA(prep_vars['org'], 
            prep_vars['outdir'], 
            prep_vars['FEBA_dir'], 
            1,
            cfg_fp=None,
            debug_bool=False, breakpoints_bool=False,
            meta_ix=7)

    # Getting time for log after program
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    logging.info("Finished running BarSeqR {}:----------".format(dt_string))


    return all_vars 
    



