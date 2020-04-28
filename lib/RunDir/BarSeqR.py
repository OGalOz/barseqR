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
# Where is FindCompound GetMediaComponents and GetMixComponents from?


def main_run(config_fp, inp_arg_list, this_file_dir):
    # Initialize dict which contains all important variables
    all_vars = {}

    config_dict, all_vars = get_config_dict(config_fp, all_vars, this_file_dir)

    all_vars = get_usage_str(config_dict, all_vars)

    all_vars = get_args(all_vars, inp_arg_list)

    # Checking to make sure input files and dirs exist
    all_vars = check_input_dirs(all_vars)


    # LoadCompounds and LoadMedia come from lib/Compounds.pm
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
    metadata_set_file_check(all_vars)

    # Build and check list of experiments for each set
    all_vars = build_and_check_sets_experiments(all_vars)

    # Map strains to genes, compute f, etc.
    all_vars = map_strains_to_genes_compute_f(all_vars)

    # Test run stops
    if "test" in all_vars and all_vars["test"] is not None:
        logging.info("All poolcount file headers verified \n")
        return 0

    # Begin writing all.poolcount string - no change to all_vars
    all_vars = init_all_poolcount_str(all_vars)


    #raise Exception("Stop")

    # Below is the heaviest workload function
    # Read data rows, combine and check correctness of metadata
    # Continue writing to all.poolcount string
    all_vars = combine_data_rows(all_vars)

    # Write to all poolcounts and close set file handles
    write_to_all_poolcount_and_close_filehandles(all_vars)

    # Write outdir/exps with only the exps for these sets and with
    # non-empty Description
    write_exps(all_vars)

    # Move strainusage and genes files to outdir
    copy_pool_genes_strain_usage(all_vars)

    # Running the R command
    run_R_cmd(all_vars)

    # Printing vars
    # del all_vars['setFh']
    # logging.info("All Vars: \n" + json.dumps(all_vars, indent=4))

    return all_vars


def get_args(all_vars, inp_arg_list):
    parser = argparse.ArgumentParser(
        usage=all_vars["usage_str"], description=all_vars["usage_str"]
    )
    parser.add_argument("-org", type=str, required=True)
    parser.add_argument("-indir", type=str, default=None)
    parser.add_argument("-metadir", type=str, default=None)
    parser.add_argument("-exps", type=str, default=None)
    parser.add_argument("-genesfile", type=str, default=None)
    parser.add_argument("-poolfile", type=str, default=None)
    parser.add_argument("-outdir", type=str, default=None)
    parser.add_argument("-sets_dir", type=str, default=None)
    # Sets should all be in -sets_dir
    parser.add_argument("-sets", type=str, nargs="+")
    # FLAGS
    # if noR or test present then their value will be 1
    parser.add_argument("-noR", action="store_const", const="1", default="0")
    parser.add_argument("-test", action="store_const", const="1", default="0")
    parser.add_argument("-feba_strain_usage", action="store_const", const="1",
            default="0")

    # Running argparser
    args = parser.parse_args(inp_arg_list)

    # Getting results
    all_vars["org"] = args.org
    all_vars["indir"] = args.indir
    all_vars["metadir"] = args.metadir
    all_vars["exps"] = args.exps
    all_vars["genesfile"] = args.genesfile
    all_vars["poolfile"] = args.poolfile
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


def run_R_cmd(all_vars):
    RscriptFP = os.path.join(all_vars['this_dir'], all_vars['R_path'])
    # Rscript is used on OSX at least to run R scripts
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
            raise Exception("R failed, see outdir log\n")
        logging.info("Succesfully ran R")


def copy_pool_genes_strain_usage(all_vars):
    # Write to outdir/pool, outdir/genes
    outdir = all_vars["outdir"]
    indir = all_vars["indir"]
    shutil.copyfile(all_vars["poolfile"], os.path.join(outdir, "pool"))
    shutil.copyfile(all_vars["genesfile"], os.path.join(outdir, "genes"))

    # Copy over the strain usage files if they exist if 
    # FEBA_STRAIN_USAGE is set
    if all_vars['feba_strain_usage']:
        if os.path.isfile(os.path.join(indir, "strainusage.barcodes")) and not os.environ(
            "FEBA_NO_STRAIN_USAGE"
        ):
            if not os.path.isfile(os.path.join(indir, "strainusage.genes")):
                raise Exception("A111")
            if not os.path.isfile(os.path.join(indir, "strainusage.genes12")):
                raise Exception("A112")
            logging.info("Copying over strain usage files\n")
            shutil.copyfile(os.path.join(indir, "strainusage.barcodes"), outdir)
            shutil.copyfile(os.path.join(indir, "strainusage.genes"), outdir)
            shutil.copyfile(os.path.join(indir, "strainusage.genes12"), outdir)
        else:
            os.unlink(os.path.join(outdir, "strainusage.barcodes"))
            os.unlink(os.path.join(outdir, "strainusage.genes"))
            os.unlink(os.path.join(outdir, "strainusage.genes12"))


def write_exps(all_vars):
    # Below write the columns in a reasonable order
    expCols = read_column_names(all_vars["expsfile"])
    expsout = os.path.join(all_vars["outdir"], "exps")
    expsout_str = "\t".join(expCols) + "\n"
    for exp in all_vars["exps"]:
        out = [exp[x] for x in expCols]
        expsout_str += "\t".join(out) + "\n"
    with open(expsout, "w") as g:
        g.write(expsout_str)
    logging.info("Wrote {} experiments to {}".format(len(all_vars["exps"]), expsout))

def write_to_all_poolcount_and_close_filehandles(all_vars):

    with open(all_vars["all_poolcount_fp"], "w") as f:
        f.write(all_vars["all_poolcount_str"])
    logging.info(
        "Wrote data for {} barcodes ({} in genes) to ".format(
            all_vars["nLine"], all_vars["nInGene"]
        )
        + "{}\n".format(all_vars["all_poolcount_fp"])
    )

    for k in all_vars['setFh'].keys():
        for f_list in all_vars['setFh'][k]:
            try:
                f_list[0].close()
            except Exception:
                logging.critical("Error Closing file: " + f_list[0])
                raise Exception


def combine_data_rows(all_vars):
    # We read all the data rows, combining them and checking that
    # the metadata is correct. Name each field SetName.Index
    # Only include items that are in the experiment list
    namesUsed = {(x["SetName"] + "." + x["Index"]): 1 for x in all_vars["exps"]}
    # Below if reached EOF in 1st file
    lastline = 0
    nLine = 0
    nInGene = 0
    print(all_vars['setFh'])
    #raise Exception("Stop")
    while lastline == 0:
        nFile = 0
        counts = {}
        nLine += 1
        metavalues = []
        for k in all_vars["setFh"].keys():
            # ind refers to the line num
            s, fhlist = k, all_vars["setFh"][k]
            indexes = all_vars["setIndex"][s]
            for fh_info in fhlist:
                file_handle = fh_info[0]
                file_length = fh_info[1]
                nFile += 1
                line = file_handle.readline()
                if line == "":
                    if nFile == 1:
                        lastline = 1
                    else:
                        if not (lastline == 1):
                            logging.info("Failed fhlist:")
                            print(fhlist)
                            raise Exception("Unexpected EOF in file for " + s)
                else:
                    line = line.rstrip()
                    F = line.split("\t")
                    if len(F) != (all_vars["nmeta"] + len(indexes)):
                        raise Exception(
                            "Wrong num columns at line "
                            "{} in file for {}".format(nLine, s)
                        )
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
                                    + "You may need to rerun combineBarSeq.pl "
                                    "with the new pool\n or use bin/resortPoolCount.pl to make it match\n"
                                )
                # and increment counts
                for i in range(len(indexes)):
                    index = indexes[i]
                    count = F[all_vars["nmeta"] + i]
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
        if lastline == 0:
            countsUsed = []
            for exp in all_vars["exps"]:
                key = exp["SetName"] + "." + exp["Index"]
                cu_var = 0 if key not in counts else counts[key]
                countsUsed.append(cu_var)

            try:
                newpos = int(metavalues[all_vars["POS"]])
            except ValueError:
                newpos = ''
            locusId, f = LocationToGene(
                metavalues[all_vars["SCAFFOLD"]],
                newpos,
                all_vars["genesSorted"],
            )
            if locusId != "":
                nInGene += 1

            # We add a tab separated line to all.poolcount
            all_vars["all_poolcount_str"] += "\t".join(["\t".join(metavalues), 
                locusId, str(f), "\t".join([str(x) for x in countsUsed])]) + \
                        "\n"

    if nInGene == 0:
        raise Exception(
            "\n No insertions in genes. Please check that "
            "{} contains genes, {} contains strains, ".format(
                all_vars["genesfile"], all_vars["poolfile"]
            )
            + " and that the scaffold identifiers match."
        )
    all_vars["nLine"] = nLine
    all_vars["nInGene"] = nInGene
    return all_vars


def init_all_poolcount_str(all_vars):
    allfile = os.path.join(all_vars["outdir"], "all.poolcount")
    all_vars["all_poolcount_fp"] = allfile
    allfields = "barcode rcbarcode scaffold strand pos locusId f".split(" ")
    for exp in all_vars["exps"]:
        allfields.append(exp["SetName"] + "." + exp["Index"])
    all_pool_str = "\t".join(allfields) + "\n"
    all_vars["all_poolcount_str"] = all_pool_str
    return all_vars


def map_strains_to_genes_compute_f(all_vars):
    # Here we map strains to genes, compute "f",
    # combine counts, read each file in parallel

    # Below is set to list of file handles reading counts for that set
    # In python we just use file paths
    setFh = {}
    # Below is set to list of count fields, may be more than is in the
    # counts file
    setIndex = {}
    meta = "barcode rcbarcode scaffold strand pos".split(" ")
    nmeta = len(meta)
    all_vars['nmeta'] = nmeta
    x = ["BARCODE", "RCBARCODE", "SCAFFOLD", "STRAND", "POS"]
    for i in range(len(x)):
        all_vars[x[i]] = i
    for k in all_vars["setFiles"].keys():
        s, filelist = k, all_vars["setFiles"][k]
        for f in filelist:
            # We join the set names with the directory and poolcount
            f = os.path.join(all_vars['sets_dir'], f + ".poolcount")
            with open(f, "r") as fh:
                file_length = len(fh.read().split("\n"))
            
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

            # 
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
                index = fields[nmeta:]
                setIndex[s] = index
                # check that all Indexes in @exps for this set are present
                index_dict = {x: 1 for x in index}
                for exp in all_vars["setExps"][s]:
                    if not exp["SetName"] == s:
                        raise Exception("SetName not what expected")
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

    return all_vars


def build_and_check_sets_experiments(all_vars):
    # Build list of experiments for each set
    setExps = {}
    for exp in all_vars["exps"]:
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
    test_or_proc_str = "Test found" if "test" in all_vars and \
            all_vars["test"] is not None else "Processing"
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
    # prespec sets is a boolean
    if not all_vars["prespec_sets"]:
        for pcfile in all_vars["pcfiles"]:
            if not (pcfile in all_vars["pcToSet"]) and not (re.search(r"test", pcfile)):
                logging.warning(
                    "Warning: poolcount file with no metadata: " "{}\n".format(pcfile)
                )


def find_set_files(all_vars):
    """
    Here we take sets and convert them to ".poolcount" files

    we use sets dir and just add set_name to .poolcount
    """
    
    

    # List of candidate poolcount files (pc - poolcount)
    pcfiles = []
    for s in all_vars['sets']:
        fn = s + ".poolcount"
        pcfiles.append(os.path.join(all_vars['sets_dir'], fn))

    all_vars["pcfiles"] = pcfiles
    logging.info("PoolCount Files:")
    logging.info(pcfiles)

    sets_dict = {i: 1 for i in all_vars["sets"]}
    # Below set to list of prefixes for poolcount files
    setFiles = {}
    # Below pcFile to set
    pcToSet = {}

    for s in all_vars["sets"]:
        logging.debug("S: " + s)
        pcfileThis = []
        for pcfile in pcfiles:
            if pcfile == s:
                pcfileThis.append(pcfile)
                pcToSet[pcfile] = s
            elif pcfile not in all_vars["sets"] and (
                pcfile[0 : len(s)].lower() == s.lower()
            ):
                postfix = pcfile[len(s) :]
                if (
                    (postfix == "" or re.search(r"^_?[a-zA-Z]$", postfix))
                    or (re.search(r"^_rep[a-zA-Z0-9]+$", postfix))
                    or (re.search(r"seq[a-zA-Z0-9]+$", postfix))
                    or (re.search(r"^_re$", postfix))
                ):
                    logging.info("Found extra file {} for {}\n".format(pcfile, s))
                    pcfileThis.append(pcfile)
                    pcToSet[pcfile] = s
        if len(pcfileThis) == 0:
            logging.info("No poolcount file for {}, skipping it\n".format(s))
            new_exps = []
            for exp in all_vars["exps"]:
                if exp["SetName"] != s:
                    new_exps.append(exp)
            if len(new_exps) == 0:
                raise Exception("No experiments remaining, giving up\n")
            all_vars["exps"] = new_exps
        else:
            setFiles[s] = pcfileThis

    all_vars["sets"] = [x for x in all_vars["sets"] if x in setFiles]
    all_vars["setFiles"] = setFiles
    all_vars["pcToSet"] = pcToSet
    all_vars["sets_dict"] = sets_dict

    return all_vars


def set_up_gene_vars(all_vars):
    # returns a list of dicts with each of the keys associated to values
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

    # Writes to STDERR
    CheckGeneLocations(genesSorted)

    return all_vars


def check_unknown_media_and_compounds(all_vars):
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
                exp["Drop"] = ""
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
    alpha = chr(206) + chr(177)
    # Each exp is a dict, containing keys: SetName, Description, Index,
    # Date_pool_expt_started.
    
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
    if all_vars['prespec_sets']:
        # ignore experiments not in pre-specified sets
        sets_dict = {}
        for s in all_vars["sets"]:
            sets_dict[s] = 1
        # Experiments will be a list of all experiments whose 'setName' is in
        # the sets
        new_exps = []
        for exp in all_vars["exps"]:
            if exp["SetName"] in sets_dict:
                new_exps.append(exp)
        if len(new_exps) == 0:
            raise Exception(
                "No experiments in specified sets " 
                "(having Description filled out)\n"
            )
        else:
            all_vars["exps"] = new_exps
    else:
        # Ignore tests

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
            new_set = exp["setName"]
            if not new_set in setsSeen:
                sets.append(new_set)
                setsSeen[new_set] = 1
        all_vars["sets"] = sets

    return all_vars


def run_load_compounds_load_media(all_vars):
    """
    Here we use metadir and prepare dicts for later
    """

    # Setting up the variables for LoadCompounds and LoadMedia
    # Below dict goes compound to list of [compound_name, cas, MW]
    compounds_dict = {}
    # Below dict goes from SynToKey => compound name
    synonyms_dict = {}

    compounds_dict, synonyms_dict = LoadCompounds(
        all_vars["metadir"], compounds_dict, synonyms_dict
    )

    # Below dict is media components with no match in compounds table
    unknownComponents = {}
    # Below dict is compoenent => media => 1 if it is reused
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
    In this function we check existence of dirs and files- but not set files

    We check that indir, metadir and outdir exist as directories.
    We check that FEBA_Barseq.tsv, genes.GC, and pool.n10 are in
    indir. We check that the poolfile has the right columns.
    We check that the R_Path is a real file.
    """
    # Here we find experiments file "expsfile" called FEBA_Barseq.tsv
    # in the indir. We also need a metadir directory.

    if not os.path.isdir(all_vars["indir"]):
        raise Exception("No such directory: {}".format(all_vars["indir"]))
    if not os.path.isdir(all_vars["metadir"]):
        raise Exception("No such directory: {}".format(all_vars["metadir"]))

    # Experiments file check
    all_vars["expsfile"] = os.path.join(all_vars["indir"], "FEBA_Barseq.tsv")
    if not os.path.isfile(all_vars["expsfile"]):
        raise Exception("No such experiments file: {}".format(all_vars["expsfile"]))
    # Genes file check
    all_vars["genesfile"] = os.path.join(all_vars["indir"], "genes.GC")
    if not os.path.isfile(all_vars["genesfile"]):
        raise Exception("No such genes file: {}".format(all_vars["genesfile"]))
    # Pool file check
    if all_vars["poolfile"] == None:
        all_vars["poolfile"] = os.path.join(all_vars["indir"], "pool.n10")
        if not os.path.isfile(all_vars["poolfile"]):
            all_vars["poolfile"] = os.path.join(all_vars["indir"], "pool")
        if not os.path.isfile(all_vars["poolfile"]):
            raise Exception("No poolfile {} or {}.n10".format(all_vars["poolfile"]))
    else:
        if not os.path.isfile(all_vars["poolfile"]):
            raise Exception("No such poolfile {}".format(all_vars["poolfile"]))
    # Checking columns of poolfile
    poolcols = read_col_names(all_vars["poolfile"])
    for col in ["barcode", "rcbarcode", "scaffold", "strand", "pos", "n", "nTot"]:
        if col not in poolcols:
            logging.warning(
                "Warning: no column named {} in {}\n".format(col, all_vars["poolfile"])
            )

    # Out dir check - where HTML ( html ) is written
    if not os.path.isdir(all_vars["outdir"]):
        raise Exception("No such directory: {}".format(all_vars["outdir"]))

    # Checking R script
    if (not os.path.isfile(all_vars["R_path"])) and (all_vars["noR"] is not None):
        raise Exception("No such script file: {}".format(all_vars["Rscript"]))

    return all_vars


def read_col_names(pool_fp):
    with open(pool_fp, "r") as f:
        file_str = f.read()
        split_l = file_str.split("\n")
        first_line = split_l[0]
        if len(first_line) == 0:
            first_line = split_l[1]
            logging.warning("FIRST LINE OF POOL FILE EMPTY")
    first_row = first_line.split("\t")
    return first_row


def get_config_dict(config_fp, all_vars, this_file_dir):
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
    iv = config_dict["init_vars"]
    usage_fp = os.path.join(all_vars["this_dir"], iv["usage_txt_fn"])
    with open(usage_fp, "r") as f:
        usage_str = f.read()
    all_vars["usage_str"] = usage_str
    return all_vars


def print_vars(all_vars):
    logging.info(json.dumps(all_vars, indent=4))



def test():
    logging.basicConfig(level=logging.DEBUG)
    this_dir = os.path.dirname(os.path.realpath(__file__))
    config_fp = os.path.join(this_dir, "barseqr_config_dict.json")
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
    main_run(config_fp, inp_arg_list)


if __name__ == "__main__":
    test()
