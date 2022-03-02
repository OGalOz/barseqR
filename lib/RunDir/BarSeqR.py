#! /?/python3
# Translation of BarSeqR.pl from:
# https://bitbucket.org/berkeleylab/feba/src/master/bin/BarSeqR.pl
import os
import shutil
import json
import logging
import argparse
import subprocess
import re
from RunDir.compounds import LoadCompounds, LoadMedia, FindCompound, \
        GetMediaComponents,  GetMixComponents
from RunDir.feba_utils import read_table, read_column_names
from RunDir.FindGene import LocationToGene, CheckGeneLocations






# Called by RunBarSeq in file 'run_barseqR.py' in this dir
def prepare_all_barcodecount_etc(config_fp, inp_arg_list, this_file_dir):
    """
    This is run from the file run_barseqR in this directory.

    barcodecount = poolcount

    What this program does:
       Checks all the input files
       Writes all.poolcount file (Combines all the barcodecounts)

    Args:
        config_fp: Path to config file
        inp_arg_list: list of arguments used by BarSeq
            [-org, org_name, -indir, Scratch_Dir_Input, -metadir, Fixed meta dir,
             -outdir, scratch_dir_output, -sets_dir, within scratch_dir_input, 
             -sets, set1 (sets_dir), set2 (sets_dir), set3 (sets_dir), ... ]
        this_file_dir: the path to directory 'RunDir'

    Description: Overall, we store all the important variables into a dict
                called all_vars, which we pass into every function in which
                we take out/create a new variable. Note that all the barcodecount files,
                including all.poolcount, should have the same exact number of
                lines, which should be equal to the number of lines in the 
                mutantpool file itself. How the program combines the barcodecount
                files is based on this fact.
                The variable 'prespec_sets' MUST be set to True
                within the KBase context.
                The flow of the program is best explained by explaining 
                each consecutive function:
        get_config_dict:
            We use the 'json' library to get the config loaded
            as a python dict. Also, we add the variables 
            'R_path', 'FEBA_dir' and 'this_dir' to all_vars.
            The config file is 'barseqr_config_dict.json'
        get_usage_str:
            We add the variable "usage_str"
            to all_vars. We get the usage string
            from the file under the 'usage_txt_fn'
            key in the config dict. The config file
            must be in the same directory as 'this_dir'
        get_args:
            We use python's argparser library to get the arguments
            for all these values as though they are command line arguments.
            We add these variables to all_vars:
            "org", "indir", "metadir", "exps", "genesfile", "mutantpool", 
            "outdir", "sets_dir","sets", "noR", "test", "feba_strain_usage"
            Note that the following are booleans (if True) or None (if False)
            "noR", "test", "feba_strain_usage" (the last three from above).
            The rest of the values are strings.
        check_input_dirs:
            In this function we check existence of dirs and files- 
            We check that indir, metadir and outdir exist as directories.
            We check that FEBA_Barseq.tsv, genes.GC, and pool.n10 are in
            indir. We check that the mutantpool has the right columns.
            Should we check existence of all 'barcodecount' files?
            We add the variables 'expsfile', 'genesfile' and 'mutantpool'
            to all_vars
        run_load_compounds_load_media:
            First, we use the function 'LoadCompounds' 
            to create dicts 'compounds_dict' and 'synonyms_dict'.
            Where 'compound_dict' maps a compound name to a list
            containing its name, its CAS (ID), and its molecular weight,
            and synonyms dict maps all other names of compounds
            to the original name you would find in 'compound_dict'.
            Then we run the function LoadMedia, which gives us
            three new dicts: "media_dict", "mix_dict", "mixAttr". 
            In both "media_dict" and "mix_dict", the keys are
            the names of the media or mix (str), and the values
            are lists of lists. Where the sublists contain
            compounds with amounts and units. In other words
            the two dicts contain the information for given
            media and mixes. "mixAttr", on the other hand,
            has keys that are the names of mixes, and those
            point to dicts with attributes and values.
            We add all of these new dicts to all_vars,
        read_table:
            Takes a TSV file and returns a list of dicts,
            one per non-header row, with header name
            pointing to value at that row
        clean_exps:
            We remove experiments (rows from exps file) which don't have 
            good entries in the experiments file. Experiments that 
            aren't good enough have, for example, no Description, 
            or their 'SetName' isn't part of the input sets. 
        check_unknown_media_and_compounds:
            We keep track of unknown Media and compounds,
            we don't raise Errors if these exist,
            we keep track of them and create variables
            "noMedia", "unknownMedia", and "unknownCompound"
            to track them.
        set_up_gene_vars
            We add the variables 'genes', which is the list
            resulting from the function 'read_table', i.e.
            a list of dicts, one per row, with keys being 
            column names pointing to the values.
            We also add the variable 'genesSorted', which is 
            a dict of scaffoldIds pointing to a sorted
            list of genes by their beginning position
            within the scaffold.
        find_set_files
            Here we take sets and convert them to ".poolcount" files
            we use sets dir and just add set_name to .poolcount
            Then we create two dicts:
                1. setFiles: maps set names to barcodecount filepaths
                2. pcToSet: maps barcodecount filepaths to set names (why?)
        build_and_check_sets_experiments
            First we create a variable called 'setExps',
            which is a dict that goes from set name
            to a list of experiments belonging
            to that set (set is equivalent to 'lane').
            We then double check that every experiment
            name is unique.
            Then we print to the user the total number
            of experiments, total number of genes,
            and total number of sets.
        map_strains_to_genes_open_filehandles
            First we get the number of preliminary metadata columns
            and call it 'nmeta' and store it.
            Then we create two dicts: 'setFh' and 'setIndex',
            which we both store in all_vars.
            'setFh' (python dict): 
                Maps setName -> list<list<Filehandle, Row num (int)>>
                Goes from setName to a list, whose length is the
                number of barcodecount files relating to that set, 
                with each element being a list with two elements:
                the first is the filehandle for that barcodecount file,
                the second is the length of that file.
            'setIndex' (python dict):
                Maps setName -> list<Indexes (str)>
                For each setName, if the file we are looking
                at is the first for this set, then we store the list of 
                indexes (all column names after the metadata) in
                this dict. Otherwise, if we've already seen a barcodecount
                file for this set, then we compare the existing indeces
                and double check that they match the first one.
                If they don't, then we raise an Exception because
                they should match. We also do a redundant comparison
                to exps SetName which should be removed.
        init_all_poolcount_str
            We add the variables "all_poolcount_fp" and
            "all_poolcount_fH" to all_vars.
            We initialize the all_poolcount file and write
            the headers, which are:
            barcode, rcbarcode, scaffold, strand, pos, locusId, f
            'all_poolcount_fH' contains the file handle for
            all.poolcount, and it already will have had the header
            written
        combine_data_rows
            This is where the majority of the computation occurs
            within this part of the program.
            We go through each set from the 'setFh' dict,
            which maps set name to a list of file handles (and file
            lengths) that relate to that set.
            Then we iterate over every line of the barcodecount
            files, all of which should be the same length,
            and add the info from all of them into 
            a single line of all.poolcount. So we end
            adding lines to all.poolcount when we reach
            the last line of the barcodecount files.
            Simply combine the same experiment name
            from all the barcodecount files and sets into 
            a single massive 'all.poolcount' file.
            We also add the columns 'locusId' and 'f'
            which track if an insertion was inside a gene,
            and if it was inside a gene, where in the gene
            was it inserted? If it wasn't inserted in a gene,
            both are the empty string ''.
        close_filehandles
            We simply close all the fileHandles of all the 
            barcodecount files, including the newly written all.poolcount.
        write_exps
            We rewrite the experiments file with only the 
            experiments we use and with the updated values.
            Values were updated in the function 'clean_exps'.
            We write the file to 'exps' in 'outdir'
        copy_pool_genes_strain_usage
            We copy the files 'mutantpool' and 'genesfile' to 
            'outdir' in preparation for the analysis part
            of the program.
        Then we return the variables that are useful for 
        the analysis part of the program
    """

    # Initialize dict which contains all important variables
    all_vars = {}

    config_dict, all_vars = get_config_dict(config_fp, all_vars, this_file_dir)

    all_vars = get_usage_str(config_dict, all_vars)

    all_vars = get_args(all_vars, inp_arg_list)


    # Checking to make sure input files and dirs exist
    all_vars = check_input_dirs(all_vars)


    # We add info regarding media and compounds from compounds.py
    # 5 keys added. Listed in func.
    all_vars = run_load_compounds_load_media(all_vars)

    # Creating exps variable as list
    read_table_list_input = [
        "SetName",
        "Index",
        "Description",
        "Date_pool_expt_started",
    ]

    all_vars["exps"] = read_table(all_vars["expsfile"], read_table_list_input)

    # Cleaning exps list
    all_vars = clean_exps(all_vars)

    # Checking for unknown media and unknown compounds
    all_vars = check_unknown_media_and_compounds(all_vars)

    # Set up genes related variables: genes, geneScaffolds, genesSorted
    all_vars = set_up_gene_vars(all_vars)

    # Can we find all set files?
    all_vars = find_set_files(all_vars)

    # Look for set files that are not in the metadata
    # metadata_set_file_check(all_vars)

    # Build and check list of experiments for each set
    all_vars = build_and_check_sets_experiments(all_vars)

    # Map strains to genes, compute f, etc.
    all_vars = map_strains_to_genes_open_filehandles(all_vars)

    # Test run stops
    if "test" in all_vars and all_vars["test"] is not None:
        logging.info("All barcodecount file headers verified \n")
        return 0

    # Begin writing all.poolcount string - no change to all_vars
    all_vars = init_all_poolcount_str(all_vars)



    # Below is the heaviest workload function
    # Read data rows, combine and check correctness of metadata
    # Continue writing to all.poolcount string
    all_vars = combine_data_rows(all_vars)

    # Write to all.poolcounts and close set file handles
    close_filehandles(all_vars)

    # Write outdir/exps with only the exps for these sets and with
    # non-empty Description
    write_exps(all_vars)

    # Move strainusage and genes files to outdir
    copy_pool_genes_strain_usage(all_vars)

    # Running the R command
    #run_R_cmd(all_vars)

    # Printing vars
    # del all_vars['setFh']
    # logging.info("All Vars: \n" + json.dumps(all_vars, indent=4))

    ret_d = {
            "outdir": all_vars["outdir"],
            "FEBA_dir": all_vars["FEBA_dir"],
            "org": all_vars["org"]
    }


    return ret_d


def get_args(all_vars, inp_arg_list):
    """
    TD:
        We should change this to just be passed as a dict
        within KBase.

    Description:
        We use python's argparser library to get the arguments
        for all these values as though they are command line arguments.
        We add these variables to all_vars:
        "org", "indir", "metadir", "exps", "genesfile", "mutantpool", 
        "outdir", "sets_dir","sets", "noR", "test", "feba_strain_usage"
        Note that the following are booleans (if True) or None (if False)
        "noR", "test", "feba_strain_usage" (the last three from above).
        The rest of the values are strings.
    """

    parser = argparse.ArgumentParser(
        usage=all_vars["usage_str"], description=all_vars["usage_str"]
    )
    parser.add_argument("-org", type=str, required=True)
    parser.add_argument("-indir", type=str, default=None)
    parser.add_argument("-metadir", type=str, default=None)
    parser.add_argument("-exps", type=str, default=None)
    parser.add_argument("-genesfile", type=str, default=None)
    parser.add_argument("-mutantpool", type=str, default=None)
    parser.add_argument("-outdir", type=str, default=None)
    parser.add_argument("-sets_dir", type=str, default=None)
    # FLAGS
    # if noR or test present then their value will be 1
    parser.add_argument("-noR", action="store_const", const="1", default="0")
    parser.add_argument("-test", action="store_const", const="1", default="0")
    parser.add_argument("-feba_strain_usage", action="store_const", const="1",
            default="0")
    # Sets should all be in -sets_dir
    # Sets are the set names without extensions
    parser.add_argument("-sets", type=str, nargs="+")


    # Running argparser
    args = parser.parse_args(inp_arg_list)

    # Getting results
    all_vars["org"] = args.org
    all_vars["indir"] = args.indir
    all_vars["metadir"] = args.metadir
    all_vars["exps"] = args.exps
    all_vars["genesfile"] = args.genesfile
    all_vars["mutantpool"] = args.mutantpool
    all_vars["outdir"] = args.outdir
    all_vars["sets_dir"] = args.sets_dir
    all_vars["sets"] = args.sets
    if args.noR == "1":
        all_vars["noR"] = True
    else:
        all_vars["noR"] = None
    if args.test == "1":
        all_vars["test"] = True
    else:
        all_vars["test"] = None
    if args.feba_strain_usage == "1":
        all_vars["feba_strain_usage"] = True
    else:
        all_vars["feba_strain_usage"] = None

    return all_vars



def copy_pool_genes_strain_usage(all_vars):
    """
    Description:
        We copy the files 'mutantpool' and 'genesfile' to 
        'outdir' in preparation for the analysis part
        of the program.
    """
    # Write to outdir/pool, outdir/genes
    outdir = all_vars["outdir"]
    indir = all_vars["indir"]
    shutil.copyfile(all_vars["mutantpool"], os.path.join(outdir, "pool"))
    shutil.copyfile(all_vars["genesfile"], os.path.join(outdir, "genes"))

    # Copy over the strain usage files if they exist if 
    # FEBA_STRAIN_USAGE is set
    if all_vars['feba_strain_usage']:
        if os.path.isfile(os.path.join(indir, "strainusage.barcodes")) and not os.environ(
            "FEBA_NO_STRAIN_USAGE"
        ):
            if not os.path.isfile(os.path.join(indir, "strainusage.genes")):
                raise Exception("Couldn't find strainusage.genes file even " + \
                                "though strainusage.barcodes file exists.")
            if not os.path.isfile(os.path.join(indir, "strainusage.genes12")):
                raise Exception("Could not find strainusage.genes12 file " + \
                                "even though strainusage.barcodes file exists.")
            logging.info("Copying over strain usage files\n")
            shutil.copyfile(os.path.join(indir, "strainusage.barcodes"), outdir)
            shutil.copyfile(os.path.join(indir, "strainusage.genes"), outdir)
            shutil.copyfile(os.path.join(indir, "strainusage.genes12"), outdir)
        else:
            os.unlink(os.path.join(outdir, "strainusage.barcodes"))
            os.unlink(os.path.join(outdir, "strainusage.genes"))
            os.unlink(os.path.join(outdir, "strainusage.genes12"))


def write_exps(all_vars):
    """
    Description:
        We rewrite the experiments file with only the 
        experiments we use and with the updated values.
        Values were updated in the function 'clean_exps'.
        We write the file to 'exps' in 'outdir'
    """
    # Below write the columns in a reasonable order
    # The reason we take expCols from expsfile instead of
    # from the keys from the first dict in the 'exps' var 
    # list is in case there are no lines.
    expCols = read_column_names(all_vars["expsfile"])
    expsout = os.path.join(all_vars["outdir"], "exps")

    # We only write the experiments used
    with open(expsout, "w") as g:
        g.write("\t".join(expCols) + "\n")
        for exp in all_vars["exps"]:
            out = [exp[x] for x in expCols]
            g.write("\t".join(out) + "\n")
    logging.info("Wrote {} experiments to {}".format(len(all_vars["exps"]), expsout))

def close_filehandles(all_vars):
    """
    Description:
        We simply close all the fileHandles of all the 
        poolcount files, including the newly written all.poolcount.
    """

    all_vars['all_poolcount_fH'].close()
    logging.info(
        "Wrote data for {} barcodes ({} in genes) to ".format(
            all_vars["nLine"], all_vars["nInGene"]
        )
        + "{}\n".format(all_vars["all_poolcount_fp"])
    )

    for k in all_vars['setFh'].keys():
        for f_list in all_vars['setFh'][k]:
            try:
                # each f_list has filehandle in loc [0], len in loc [1]
                f_list[0].close()
            except Exception:
                logging.critical("Error Closing file: " + f_list[0])
                raise Exception


def combine_data_rows(all_vars):
    """
    Description:
        This is where the majority of the computation occurs
        within this part of the program.
        We go through each set from the 'setFh' dict,
        which maps set name to a list of file handles (and file
        lengths) that relate to that set.
        Then we iterate over every line of the poolcount
        files, all of which should be the same length,
        and add the info from all of them into 
        a single line of all.poolcount. So we end
        adding lines to all.poolcount when we reach
        the last line of the poolcount files.
        Simply combine the same experiment name
        from all the poolcount files and sets into 
        a single massive 'all.poolcount' file.
        We also add the columns 'locusId' and 'f'
        which track if an insertion was inside a gene,
        and if it was inside a gene, where in the gene
        was it inserted? If it wasn't inserted in a gene,
        both are the empty string ''.

    SubRoutines:
        LocationToGene:
            Returns [locusId (str), f (float)]
                where locusId and f are empty string '' if the insertion
                did not hit a gene, and they are the corresponding
                values if it did hit a gene.

    """

    logging.info("Combining all.poolcount files to create all.poolcount")
    
    # Not sure what below code would do
    #namesUsed = {(x["SetName"] + "." + x["Index"]): 1 for x in all_vars["exps"]}

    # Below is True if reached EOF in 1st file
    lastline = False

    # nLine refers to the line num
    nLine = 0

    nInGene = 0

    while not lastline:
        # nFile keeps track of total fileHandles passed through
        nFile = 0
        # counts keeps SetName.Index -> list of counts
        counts = {}
        nLine += 1
        if nLine % 5000 == 0:
            logging.info(f"At line # {nLine} out of {all_vars['expected_fl']}")
        metavalues = []
        # setFh dict mapping setName to list of <FileHandle, FileLength (int)>
        # relating to that set
        for k in all_vars["setFh"].keys():
            s, fhlist = k, all_vars["setFh"][k]

            # Indexes is the list of indexes in the poolcount file
            # related to this set.
            indexes = all_vars["setIndex"][s]

            for fh_info in fhlist:
                file_handle = fh_info[0]
                file_length = fh_info[1]
                nFile += 1
                line = file_handle.readline()
                if line == "":
                    if nFile == 1:
                        lastline = True
                    else:
                        if not lastline:
                            logging.info("Failed fhlist:")
                            print(fhlist)
                            raise Exception("Unexpected EOF in file for " + s)
                else:
                    # Next line from the poolcount file
                    line = line.rstrip()
                    F = line.split("\t")
                    if len(F) != (all_vars["nmeta"] + len(indexes)):
                        raise Exception(
                            "Wrong num columns at line "
                            "{} in file for {}".format(nLine, s)
                        )
                    # The meta info related to this poolcount line
                    metaThis = F[0 : (all_vars["nmeta"])]
                    # save or check metadata
                    if nFile == 1:
                        metavalues = metaThis
                    else:
                        for i in range(0, all_vars['nmeta']):
                            if not metavalues[i] == metaThis[i]:
                                raise Exception(
                                    "Non-matching {} = {} ".format(
                                        all_vars["meta"][i], metaThis[i]
                                    )
                                    + "in {} line {}, expected {}\n".format(
                                        s, nLine, metavalues[i]
                                    )
                                    + "You may need to rerun combineBarSeq" \
                                     + "with a new pool\n"
                                )

                # and increment counts
                for i in range(len(indexes)):
                    index = indexes[i]
                    count = F[all_vars["nmeta"] + i]
                    # We check if the count is an integer as expected
                    if not re.search(r"^\d+$", count):
                        raise Exception(
                            "Not a count field: "
                            "{} at line {} for set {}, i: {}".format(
                                count, nLine, s, i)
                        )
                    new_key = s + "." + index
                    if new_key in counts:
                        counts[new_key] += int(count)
                    else:
                        counts[new_key] = int(count)
        if not lastline:
            countsUsed = []
            for exp in all_vars["exps"]:
                key = exp["SetName"] + "." + exp["Index"]
                cu_var = 0 if key not in counts else counts[key]
                countsUsed.append(cu_var)
            try:
                newpos = int(metavalues[all_vars["POS"]])
            except ValueError:
                newpos = ''
            # Note locusId and f are both empty string if insertion doesn't hit 
            # a gene.
            locusId, f = LocationToGene(
                metavalues[all_vars["SCAFFOLD"]],
                newpos,
                all_vars["genesSorted"],
            )
            if locusId != "":
                nInGene += 1

            # We add a tab separated line to all.poolcount
            all_vars["all_poolcount_fH"].write("\t".join(["\t".join(metavalues), 
                locusId, str(f), "\t".join([str(x) for x in countsUsed])]) + \
                        "\n")

    if nInGene == 0:
        raise Exception(
            "\n No insertions found in genes. Please check that "
            "{} contains genes, {} contains strains, ".format(
                all_vars["genesfile"], all_vars["mutantpool"]
            )
            + " and that the scaffold identifiers match."
        )
    all_vars["nLine"] = nLine
    all_vars["nInGene"] = nInGene

    logging.info("Finished combining all barcodecount files to create all.poolcount")
    return all_vars


def init_all_poolcount_str(all_vars):
    """
    Description:
        We add the variables "all_poolcount_fp" and
        "all_poolcount_fH" to all_vars.
        We initialize the all_poolcount file and write
        the headers, which are:
        barcode, rcbarcode, scaffold, strand, pos, locusId, f
        'all_poolcount_fH' contains the file handle for
        all.poolcount, and it already will have had the header
        written
    """
    allfile = os.path.join(all_vars["outdir"], "all.poolcount")
    all_vars["all_poolcount_fp"] = allfile
    allfields = "barcode rcbarcode scaffold strand pos locusId f".split(" ")
    for exp in all_vars["exps"]:
        allfields.append(exp["SetName"] + "." + exp["Index"])

    # We initiate file handle for all.poolcount
    allfile_handle = open(allfile, "w")
    allfile_handle.write("\t".join(allfields) + "\n")
    all_vars["all_poolcount_fH"] = allfile_handle
    return all_vars


def map_strains_to_genes_open_filehandles(all_vars):
    """
    all_vars:
        nmeta: Number of metadata columns (?)
    Description:
        First we get the number of preliminary metadata columns
        and call it 'nmeta' and store it.
        Then we create two dicts: 'setFh' and 'setIndex',
        which we both store in all_vars.
        'setFh' (python dict): 
            Maps setName -> list<list<Filehandle, Row num (int)>>
            Goes from setName to a list, whose length is the
            number of poolcount files relating to that set, 
            with each element being a list with two elements:
            the first is the filehandle for that poolcount file,
            the second is the length of that file.
        'setIndex' (python dict):
            Maps setName -> list<Indexes (str)>
            For each setName, if the file we are looking
            at is the first for this set, then we store the list of 
            indexes (all column names after the metadata) in
            this dict. Otherwise, if we've already seen a poolcount
            file for this set, then we compare the existing indeces
            and double check that they match the first one.
            If they don't, then we raise an Exception because
            they should match. We also do a redundant comparison
            to exps SetName which should be removed.
    TD: Eventually might be better to convert this into pandas

    """

    meta = "barcode rcbarcode scaffold strand pos".split(" ")
    nmeta = len(meta)
    all_vars['nmeta'] = nmeta


    # Below is set to list of file handles reading counts for that set
    setFh = {}
    # Below is set to list of count fields, may be more than is in the
    # counts file
    setIndex = {}

    # We create this to have index loc for below keys within files
    x = ["BARCODE", "RCBARCODE", "SCAFFOLD", "STRAND", "POS"]
    for i in range(len(x)):
        all_vars[x[i]] = i

    expected_file_length = None
    for k in all_vars["setFiles"].keys():
        # Below s is setName, filelist is list of poolcount fps related
        s, filelist = k, all_vars["setFiles"][k]
        for f in filelist:

            # This is just to keep track of length of file
            with open(f, "r") as fh:
                file_length = len(fh.read().split("\n"))
                if expected_file_length is None:
                    expected_file_length = file_length
                if file_length != expected_file_length:
                    logging.warning(f"For file {f}, file size doesn't match the first poolcount."
                                    f" {f} filesize: {file_length}. First: {expected_file_length}")
            
            # Note we open files here! Close after writing all.poolcount
            fh = open(f,"r")

            # We read the first line and store it for later.
            first_line = fh.readline()

            # We add the filehandle and file_length to setFh
            if s in setFh:
                setFh[s].append([fh, file_length])
            else:
                setFh[s] = [[fh, file_length]]

            # We use the first line from before
            first_line = first_line.rstrip()

            # Poolcount header:
            fields = first_line.split("\t")
            if not (len(fields) >= 6):
                raise Exception("Too few fields in " + f)
            for i in range(nmeta):
                if not fields[i] == meta[i]:
                    raise Exception(
                        "Expected {} but field is {}".format(meta[i], fields[i])
                    )
            if not (s in setIndex):
                # first file for this set
                index_list = fields[nmeta:]
                # setIndex of setName is mapped to all the indexes
                setIndex[s] = index_list
                # check that all Indexes in @exps for this set are present
                index_dict = {x: 1 for x in index_list}


                for exp in all_vars["setExps"][s]:
                    # Each exp is a hash of Exp file column to value

                    # This just makes sure earlier processes are correct
                    if not exp["SetName"] == s:
                        raise Exception("SetName not what expected")
                    # Here is where we check if indexes align between 
                    # poolcount file and experiments file
                    if not exp["Index"] in index_dict:
                        logging.warning(
                            "WARNING! No field for indx {}".format(exp["Index"])
                            + " in "
                            + all_vars["indir"]
                            + "/"
                            + f
                            + ".poolcount!n"
                        )
            else:
                # Additional file for this set
                # expect is the list of expected indices
                expect = setIndex[s]
                for i in range(len(expect)):
                    field = fields[nmeta + i]
                    if not field == expect[i]:
                        raise Exception(
                            "Expected {} from first ".format(expect[i])
                            + "file but see {} in ".format(field)
                            + all_vars["indir"]
                            + "/"
                            + f
                            + ".poolcount"
                        )

    all_vars["setFh"] = setFh
    all_vars["setIndex"] = setIndex
    all_vars["expected_fl"] = expected_file_length

    return all_vars


def build_and_check_sets_experiments(all_vars):
    """

    Description:
        First we create a variable called 'setExps',
        which is a dict that goes from set name
        to a list of experiments belonging
        to that set (set is equivalent to 'lane').
        We then double check that every experiment
        name is unique.
        Then we print to the user the total number
        of experiments, total number of genes,
        and total number of sets.
    """
    # Build list of experiments for each set
    setExps = {}
    for exp in all_vars["exps"]:
        # Each exp is a hash of column Names to column values from 
        # the experiments file.
        if exp["SetName"] in setExps:
            setExps[exp["SetName"]].append(exp)
        else:
            setExps[exp["SetName"]] = [exp]

    # And check that each index is unique for each set
    for k in setExps.keys():
        s, exps = k, setExps[k]
        indexSeen = {}
        for exp in exps:
            index = exp["Index"]
            if index in indexSeen:
                raise Exception(
                    "Duplicate experiment entries "
                    "for index {} in set {}".format(index, s)
                )
            else:
                indexSeen[index] = 1

    all_vars["setExps"] = setExps

    # Print Info:
    test_or_proc_str = "Test found: " if "test" in all_vars and \
            all_vars["test"] is not None else "Processing: "


    logging.info(
        "{} {} experiments, {} genes, {} sets for {}\n".format(
            test_or_proc_str,
            len(all_vars["exps"]),
            len(all_vars["genes"]),
            len(all_vars["sets"]),
            all_vars['org']
        )
    )

    return all_vars


def metadata_set_file_check(all_vars):
    """
    Description:
        If prespec_sets is not set to True, then
        we check if each poolcount file is in 'pcToSet' dict
        (why?)
    """
    # prespec sets is a boolean
    if not all_vars["prespec_sets"]:
        for pcfile in all_vars["pcfiles"]:
            if not (pcfile in all_vars["pcToSet"]) and not (re.search(r"test", pcfile)):
                logging.warning(
                    "Warning: poolcount file with no metadata: " "{}\n".format(pcfile)
                )


def find_set_files(all_vars):
    """
    Args:
        all_vars (d)
            sets:
    Description:
        Here we take sets and convert them to ".poolcount" files
        we use sets dir and just add set_name to .poolcount
        Then we create two dicts:
            1. setFiles: maps set names to poolcount filepaths
            2. pcToSet: maps poolcount filepaths to set names (why?)
    """
    
    # Below set to list of prefixes for poolcount files
    setFiles = {}
    # Below pcFile to set
    pcToSet = {}

    # We are assuming there is only one poolcount file with setname
    # + ".poolcount"
    for s in all_vars['sets']:
        fn = s + ".poolcount"
        pcfp = os.path.join(all_vars['sets_dir'], fn)
        setFiles[s] = [pcfp]
        pcToSet[pcfp] = s

    all_vars["setFiles"] = setFiles
    all_vars["pcToSet"] = pcToSet

    return all_vars


def set_up_gene_vars(all_vars):
    """
    Description:
        We add the variables 'genes', which is the list
        resulting from the function 'read_table', i.e.
        a list of dicts, one per row, with keys being 
        column names pointing to the values.
        We also add the variable 'genesSorted', which is 
        a dict of scaffoldIds pointing to a sorted
        list of genes by their beginning position
        within the scaffold.
    """
    genes = read_table(
        all_vars["genesfile"],
        ("locusId scaffoldId sysName begin end strand".split(" ")),
    )
    geneScaffolds = {}
    # scaffold to list of genes sorted by begin
    genesSorted = {}
    for g in genes:
        geneScaffolds[g["scaffoldId"]] = 1
        if g["scaffoldId"] in genesSorted:
            genesSorted[g["scaffoldId"]].append(g)
        else:
            genesSorted[g["scaffoldId"]] = [g]
    for scaffold in genesSorted.keys():
        # genesSorted[scaffold] is a list of dicts, which we'll sort by 'begin'
        sorted_list = sorted(genesSorted[scaffold], key=lambda i: i["begin"])
        genesSorted[scaffold] = sorted_list

    all_vars["genes"], all_vars["genesSorted"] = genes, genesSorted

    CheckGeneLocations(genesSorted)

    return all_vars


def check_unknown_media_and_compounds(all_vars):
    """
    Args:
        all_vars:
            'exps': list<dict>, a 'read' table, with one dict per row
                    of 'exps' file, dict has columns pointing to values
                    at that row.
                    Cols must include "Index", "SetName"

    Description:
        We keep track of unknown Media and compounds,
        we don't raise Errors if these exist, only
        we keep track of them and create variables
        "noMedia", "unknownMedia", and "unknownCompound".
    """
    unknownMedia = {}
    unknownCompound = {}
    noMedia = []
    for exp in all_vars["exps"]:
        # Ignore lines which are not filled out or dropped
        if exp["Index"] == "" and exp["SetName"] == "":
            continue
        if "Drop" in exp:
            drop = exp["Drop"]
            drop = drop.replace(" ", "")
            if drop.upper() == "TRUE" or drop.upper() == "DROP":
                logging.info(
                    "Dropping {} {} {}\n".format(
                        exp["SetName"], exp["Index"], exp["Description"]
                    )
                )
                exp["Drop"] = "TRUE"
            elif drop.upper() == "NA":
                exp["Drop"] = "" # Why not FALSE
            elif drop.upper() != "FALSE":
                logging.warning(
                    "Unknown drop code {} for {} {}\n".format(
                        drop, exp["SetName"], exp["Index"]
                    )
                )
                exp["Drop"] = ""

        # Clean up media and warn if not a known media
        if "Media" in exp and exp["Media"] != "":
            # Removing white space from beginning and end
            exp["Media"] = exp["Media"].strip()
            gmd_results = GetMediaComponents(exp["Media"], 
                    all_vars['media_dict'],
                    all_vars['mixAttr'],
                    all_vars['mix_dict'])
            if (gmd_results) is not None:
                unknownMedia[exp["Media"]] = 1
        else:
            noMedia.append(exp["SetName"] + "." + exp["Index"])

        # Clean up Group
        if "Group" in exp and exp["Group"] is not None:
            exp["Group"] = exp["Group"].strip()
            if not (exp["Group"] == "pH" or exp["Group"] == "Time0"):
                exp["Group"] = exp["Group"].lower()

        # Clean up condition_1,2,3 and warn if not a known compound --
        # but skip this for Time0 samples as sometimes they have condition
        # set to Time0
        for fld in ["Condition_1", "Condition_2", "Condition_3", "Condition_4"]:
            if fld in exp and exp["Group"] != "Time0":
                exp[fld] = fld.replace(chr(206) + chr(177), "a")
                exp[fld] = exp[fld].strip()
                if exp[fld].lower() == "none" or exp[fld] == "NA":
                    exp[fld] = ""
                if exp[fld] != "":
                    compound = FindCompound(exp[fld], 
                            all_vars['compounds_dict'],
                            all_vars['synonyms_dict'])
                    gmd_results =  GetMediaComponents(exp["Media"], 
                    all_vars['media_dict'],
                    all_vars['mixAttr'],
                    all_vars['mix_dict'])
                    gmx_results = GetMixComponents(exp[fld], 
                            all_vars['mix_dict'])

                    if compound == None and (
                        gmd_results is not None
                        and gmx_results  is not None
                    ):
                        compound = exp[fld]
                    if compound == None:
                        if not re.search(r" exudates?$", exp[fld]) or (
                            re.search(r"^supernatant; ", exp[fld], re.IGNORECASE)
                        ):
                            unknownCompound[exp[fld]] = 1
                    else:
                        exp[fld] = compound

    all_vars["noMedia"] = noMedia
    all_vars["unknownMedia"] = unknownMedia
    all_vars["unknownCompound"] = unknownCompound
    if len(noMedia) > 0:
        logging.info("No media for:\n {}\n".format("\t".join(noMedia)))
    if len(unknownMedia.keys()) > 0:
        logging.info(
            "Unknown media:\n {}\n".format("\t".join(sorted(unknownMedia.keys())))
        )
    if len(unknownCompound.keys()) > 0:
        logging.info(
            "Unknown Compound:\n {}\n".format("\t".join(sorted(unknownCompound.keys())))
        )

    return all_vars


def clean_exps(all_vars):
    """
    Args:
        all_vars:
            exps: (list of dicts, one per row, header name pointing to value)

    Description:
        We remove experiments (rows from exps file) which don't have 
        good entries in the experiments file. Experiments that 
        aren't good enough have, for example, no Description, 
        or their 'SetName' isn't part of the input sets. 
    """
    # We get the chars for an alpha symbol
    alpha = chr(206) + chr(177)
    # Each exp is a dict, containing keys: SetName, Description, Index,
    # and Date_pool_expt_started.
    
    # We keep track of new experiments which are cleaned up
    new_exp = []
    for exp in all_vars["exps"]:
        # Removing cases where Description is empty
        if exp["Description"] == "":
            logging.debug("Empty Description")
            continue

        # extra spaces at end are common data entry errors
        exp["Description"] = exp["Description"].rstrip()
        exp["Description"] = exp["Description"].replace(alpha, "a")
        exp["SetName"] = exp["SetName"].rstrip()
        exp["Index"] = exp["Index"].rstrip()
        new_exp.append(exp)

    logging.info("Number of experiments: "  + str(len(all_vars["exps"])))
    logging.info(all_vars['sets'])
    all_vars["exps"] = new_exp
    all_vars['prespec_sets'] = len(all_vars["sets"]) > 0

    if not all_vars['prespec_sets']:
        raise Exception("Running on KBase and there are no poolcount files."
                        " Cannot run.")

    # In KBase we will almost always go into this first if statement
    if all_vars['prespec_sets']:

        # ignore experiments not in pre-specified sets
        sets_dict = {s:1 for s in all_vars["sets"]}

        # Experiments will be a list of all experiments whose 'SetName' is in
        # the sets
        new_exps = []
        for exp in all_vars["exps"]:
            if exp["SetName"] in sets_dict:
                new_exps.append(exp)
        if len(new_exps) == 0:
            logging.warning("No overlap between experiment set names and "
                            "PoolCount setNames. Printing both:")
            logging.warning("Experiments: ")
            logging.warning(all_vars["exps"])
            logging.warning("PoolCount Sets:")
            logging.warning(all_vars["sets"])
            raise Exception(
                "No experiments in specified sets " 
                "having Description filled out)\n"
                "Check your Experiments File and "
                "specifically the SetName column for"
                " the sets (poolcount files) you "
                "are testing. "
                "Make sure the PoolCount file's (set's)"
                " name is exactly the same as the setNames"
                " in the Experiments file."
            )
        else:
            all_vars["exps"] = new_exps
    else:
        # Ignore tests - any experiment with "test" in SetName is ignored.

        sets_dict = {exp["SetName"]: 1 for exp in all_vars["exps"]}

        for s in sets_dict.keys():
            if re.search(r"test", s):
                logging.info("Ignoring SetName = {}\n".format(s))

        new_exps = []
        for exp in all_vars["exps"]:
            if not (re.search(r"test", exp["SetName"])):
                new_exps.append(exp)
        if len(new_exps) == 0:
            raise Exception(
                "No experiments after filtering out empty "
                "Description and test sets.\n"
            )
        all_vars["exps"] = new_exps

        # We set up sets_list
        setsSeen = {}
        sets = []
        for exp in all_vars["exps"]:
            new_set = exp["SetName"]
            if not new_set in setsSeen:
                sets.append(new_set)
                setsSeen[new_set] = 1
        all_vars["sets"] = sets

    return all_vars


def run_load_compounds_load_media(all_vars):
    """
    Here we use metadir and prepare dicts for later

    Keys added:
        media_dict: (d)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mix_dict: (d) (Like media_dict)
            media (str) -> list<compound_l>
                where compound_l list<compound (str), concentration (str), units (str)>
                e.g. [Ammonium chloride, 0.25, g/L]
        mixAttr: (d)
            media_name (str) -> attr_d (d)
                attribute (str) -> value (str) e.g.
                    Description -> Defined minimal media for soil and groundwater 
                                    bacteria with glucose and MES buffer
                    or
                    Minimal -> TRUE
        compounds: (d)
            compound_name -> [compound_name, CAS (str), MW (str)]
        synonyms: (d)
            synonym_name (str) -> compound_name (str)

    Vars:
        LoadMediaResultsDict (python dict with keys:):
            "media_dict": media,
            "mix_dict": mix,
            "mixAttr": mixAttr,
            "compounds_dict": compounds,
            "synonyms_dict": synonyms

        

    Description:
        First, we use the function 'LoadCompounds' 
        to create dicts 'compounds_dict' and 'synonyms_dict'.
        Where 'compound_dict' maps a compound name to a list
        containing its name, its CAS (ID), and its molecular weight,
        and synonyms dict maps all other names of compounds
        to the original name you would find in 'compound_dict'.
        Then we run the function LoadMedia, which gives us
        three new dicts: "media_dict", "mix_dict", "mixAttr". 
        In both "media_dict" and "mix_dict", the keys are
        the names of the media or mix (str), and the values
        are lists of lists. Where the sublists contain
        compounds with amounts and units. In other words
        the two dicts contain the information for given
        media and mixes. "mixAttr", on the other hand,
        has keys that are the names of mixes, and those
        point to dicts with attributes and values.
        We add all of these new dicts to all_vars,

    """

    # Setting up the variables for LoadCompounds and LoadMedia
    # Below dict goes compound to list of [compound_name, cas, MW]
    # where cas is a string ID of compound, MW is Molar Weight
    compounds_dict = {}
    # Below dict goes from SynToKey => compound name
    synonyms_dict = {}

    compounds_dict, synonyms_dict = LoadCompounds(
        all_vars["metadir"], compounds_dict, synonyms_dict
    )

    # Below dict is media components with no match in compounds table
    unknownComponents = {}

    # Below dict is component => media => 1 if it is reused
    reuseComponents = {}

    LoadMediaResultsDict = LoadMedia(
        all_vars["metadir"],
        compounds_dict,
        synonyms_dict,
        unknownComponents,
        reuseComponents,
    )

    #updates the keys "media_dict", "mix_dict", "mixAttr", "compounds_dict",
    # "synonyms_dict"
    all_vars.update(LoadMediaResultsDict)

    return all_vars


def check_input_dirs(all_vars):
    """
    Args:
        all_vars (python dict): Must contain keys:
            indir: (str path to directory)
            metadir: (str path to directory)

    Returns:
        Adds the following args to all_vars:
            expsfile
            genesfile
            mutantpool

    Description:
        In this function we check existence of dirs and files- 
        We check that indir, metadir and outdir exist as directories.
        We check that FEBA_Barseq.tsv, genes.GC, and pool.n10 are in
        indir. We check that the mutantpool has the right columns.
        Should we check existence of all 'poolcount' files?
        We add the variables 'expsfile', 'genesfile' and 'mutantpool'
        to all_vars
    """
    # Here we find experiments file "expsfile" called FEBA_Barseq.tsv
    # in the indir. We also need a metadir directory.

    if not os.path.isdir(all_vars["indir"]):
        raise Exception("No such directory: {}".format(all_vars["indir"]))

    indir_files = os.listdir(all_vars["indir"])
    for x in ["FEBA_Barseq.tsv", "genes.GC"]:
        if x not in indir_files:
            raise Exception(f"File {x} missing from indir. Current files "
                            " in indir:\n" + ", ".join(indir_files))

    all_vars["expsfile"] = os.path.join(all_vars["indir"], "FEBA_Barseq.tsv")
    all_vars["genesfile"] = os.path.join(all_vars["indir"], "genes.GC")
    if "pool.n10" not in indir_files:
        if "pool" not in indir_files:
            raise Exception("Neither 'pool.n10' nor 'pool' found "
                            " in indir files. Current files "
                            " in indir:\n" + ", ".join(indir_files))
        else:
            all_vars["mutantpool"] = os.path.join(all_vars["indir"], "pool")
    else:
        all_vars["mutantpool"] = os.path.join(all_vars["indir"], "pool.n10")


    if not os.path.isdir(all_vars["metadir"]):
        raise Exception("No such directory: {}".format(all_vars["metadir"]))
    else:
        metadir_files = os.listdir(all_vars["metadir"])
        for x in []:
            if x not in metadir_files:
                raise Exception(f"File {x} missing from metadir. Current "
                                " files in metadir:\n" + ", ".join(metadir_files))

    # Checking columns of mutantpool
    poolcols = read_col_names(all_vars["mutantpool"])
    for col in ["barcode", "rcbarcode", "scaffold", "strand", "pos", "n", "nTot"]:
        if col not in poolcols:
            logging.warning(
                "Warning: no column named {} in {}\n".format(col, all_vars["mutantpool"])
            )

    # Out dir check - where HTML ( html ) is written
    if not os.path.isdir(all_vars["outdir"]):
        raise Exception("No such directory: {}".format(all_vars["outdir"]))

    return all_vars


def read_col_names(pool_fp):
    with open(pool_fp, "r") as f:
        first_line = f.readline()
        if first_line == "":
            raise Exception(f"First line of pool file {pool_fp} is empty.")
    first_row = first_line.split("\t")
    return first_row


def get_config_dict(config_fp, all_vars, this_file_dir):
    """
    Description:
        We use the 'json' library to get the config loaded
        as a python dict. Also, we add the variables 
        'R_path', 'FEBA_dir' and 'this_dir' to all_vars.
        The config file is 'barseqr_config_dict.json'
    """
    with open(config_fp, "r") as f:
        config_dict = json.loads(f.read())

    all_vars["this_dir"] = this_file_dir 

    # Setting R run configuration paths
    init_vars = config_dict['init_vars']
    all_vars['R_path'] = os.path.join(all_vars['this_dir'], 
                                      init_vars['R_path'])
    all_vars['FEBA_dir'] = os.path.join(all_vars['this_dir'], 
                                        init_vars['FEBA_dir'])

    return config_dict, all_vars


def get_usage_str(config_dict, all_vars):
    """
    Description:
        We add the variable "usage_str"
        to all_vars. We get the usage string
        from the file under the 'usage_txt_fn'
        key in the config dict. The config file
        must be in the same directory as 'this_dir'
    """
    iv = config_dict["init_vars"]
    usage_fp = os.path.join(all_vars["this_dir"], iv["usage_txt_fn"])
    with open(usage_fp, "r") as f:
        usage_str = f.read()
    all_vars["usage_str"] = usage_str
    return all_vars


def print_vars(all_vars):
    logging.info(json.dumps(all_vars, indent=4))




def runBarSeqPy(all_vars):
    """
    Args:
    Description:
        We run the analysis part of the program
    """
    if 'bsPy_cfg_fp' not in all_vars:
        raise Exception("No BarSeqPy config indicated.")

    RunFEBA(all_vars['org'], 
            all_vars['outdir'], 
            all_vars['FEBA_dir'], 
            1,
            all_vars['bsPy_cfg_f'],
            debug_bool=False, breakpoints_bool=False,
            meta_ix=7)





def ret_test_inp_arg_list():

    inp_arg_list = [
        "-org",
        "SB2B",
        "-indir",
        "Test_Files",
        "-metadir",
        "metadata",
        "-outdir",
        "tmp",
        "-sets_dir",
        "Test_Files",
        "-sets",
        "SB2B_ML5_set1",
        "SB2B_ML5_set2",
        "SB2B_ML5_set3"
    ]

    return inp_arg_list



def main():
    logging.basicConfig(level=logging.DEBUG)
    help_str = "python3 BarSeqR.py inp_arg_list.json op_dir 1"
    args = sys.argv
    if args[-1] != "1":
        print(help_str)
        sys.exit(0)
    else:
        this_dir = os.path.dirname(os.path.realpath(__file__))
        config_fp = os.path.join(this_dir, "barseqr_config_dict.json")
        sample_inp_arg_list = json.loads(open(args[1]).read())
        prepare_all_barcodecount_etc(config_fp, sample_inp_arg_list, this_dir)


if __name__ == "__main__":
    main()


# Deprecated:
def run_R_cmd(all_vars):
    """ This function calls the FEBA R script entry point 'RunFEBA.R'

        all_vars:
            outdir: (str) path
                dir must contain 'all.poolcount', 'exps'
                'genes' 'pool'    files
            noR: bool
            FEBA_dir: (str) path
            org: str

    """
    RscriptFP = os.path.join(all_vars['this_dir'], all_vars['R_path'])
    # Rscript is used on to run R scripts
    R_exec_path = "Rscript"
    # RscriptFP is location of file, org is organism, outdir is data dir
    Rcmds = [R_exec_path, RscriptFP, all_vars['org'], all_vars['outdir'], 
            all_vars['FEBA_dir']]
    R_op_path = os.path.join(all_vars["outdir"], ".FEBA.success")
    if "noR" in all_vars and all_vars["noR"] is True:
        logging.info("Skipping the R step: {}\n".format(" ".join(Rcmds)))
    else:
        logging.info("Running R: {}\n".format(Rcmds))
        if os.path.isfile(R_op_path):
            os.unlink(R_op_path)
        log_output_fp = os.path.join(all_vars['outdir'], "log")
        f = open(log_output_fp, "w")
        subprocess.call(Rcmds, stdout=f)
        f.close()
        if not os.path.exists(R_op_path):
            logging.critical(f"R output was not written at {R_op_path}")
            #raise Exception("R failed, see outdir log\n")
        else:
            logging.info("Succesfully ran R")

