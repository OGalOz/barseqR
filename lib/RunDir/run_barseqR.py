#!python3

import os
from datetime import datetime
import logging
from RunDir.BarSeqR import main_run 


# Following Function primarily used to set location of config for running BarSeq
# All important files and directories for running must be the same as this file.
def RunBarSeq(arg_list):
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    cfg_fp = os.path.join(this_file_dir, "barseqr_config_dict.json")

    # Getting time for log before program
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    logging.info("Beginning to run BarSeqR {}:-----------".format(dt_string))

    # Run Program Here
    all_vars = main_run(cfg_fp, arg_list, this_file_dir)


    # Getting time for log after program
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    logging.info("Finished running BarSeqR {}:----------".format(dt_string))


    return all_vars 
    



