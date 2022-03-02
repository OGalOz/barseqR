#!python3
import re
import logging
import os
import json


def validate_params(params):
    """What do we expect params to have as keys:
        genome_ref
        mutantpool_ref
        exps_ref
        sets_ref
        output_name
        workspace_name

    Returns: ^ with strings as values for all
    """

    for p in ["genome_ref", "mutantpool_ref", "exps_ref"]:
        if p not in params:
            raise Exception(p + " missing from params.")
        elif len(params[p].split('/')) != 3:
            raise Exception(p + " doesn't have form A/B/C as expected:\n"
                    "{}".format(params[p]))

    if "sets_refs" not in params or len(params["sets_refs"]) < 1:
        raise Exception("No set files (.poolcount) found")

    # Edit output name to have no spaces and check for weird characters
    if 'output_name' not in params:
        params['output_name'] = 'Untitled'
    else:
        params['output_name'] = check_output_name(params['output_name'])

    advanced_params_list = [
            "okControls",
            "okDay",
            "okLane",
            "drop_exps",
            "compute_cofit_bool",
            "compute_spfc_bool",
            "compute_High_bool",
            "minSampleReads" ,
            "minGenesPerScaffold" ,
            "minT0Strain" ,
            "minT0Gene" ,
            "minGenesAllowed",
            "minGenesUsed12" ,
            "norm_median_window",
            "min_gMed" ,
            "max_mad12" ,
            "min_cor12" ,
            "max_gccor" ,
            "max_adjcor" ,
            "nTopCofit" ,
            "minCofitExp" ,
            "Spfc_minT" ,
            "Spfc_minFit" ,
            "Spfc_percentile" ,
            "Spfc_percentileFit" ,
            "Spfc_minDelta" ,
            "High_min_fit" ,
            "High_min_t" ,
            "High_max_se" ,
            "High_min_reads" ,
            "High_min_gMean" ,
            "High_max_below" ,
            "High_min_strains" ,
            "High_min_strain_fraction" ,
            "Strong_lr",
            "Strong_t",
    ]

    for x in advanced_params_list:
        if x not in params:
            raise Exception(f"Key {x} missing from params. "
                            "Params keys: " + ",".join(params.keys()))


    for x in ["okControls",
            "okDay",
            "okLane",
            "drop_exps",
            "compute_cofit_bool",
            "compute_spfc_bool",
            "compute_High_bool"]:
        if params[x] not in ["True", "true", True, "False", "false", False]:
            raise Exception(f"Key {x} in params should be a boolean, instead: "
                            f"{params[x]}")
        else:
            if not isinstance(params[x], bool):
                if params[x].upper() == "TRUE":
                    params[x] = True
                elif params[x].upper() == "FALSE":
                    params[x] = False
                else:
                    raise Exception(f"Cannot process param {x}, value is {params[x]}")


    for x in ["minSampleReads" ,
            "minGenesPerScaffold" ,
            "minT0Strain" ,
            "minT0Gene" ,
            "minGenesAllowed",
            "minGenesUsed12" ,
            "norm_median_window",
            "min_gMed" ,
            "max_mad12" ,
            "min_cor12" ,
            "max_gccor" ,
            "max_adjcor" ,
            "nTopCofit" ,
            "minCofitExp" ,
            "Spfc_minT" ,
            "Spfc_minFit" ,
            "Spfc_percentile" ,
            "Spfc_percentileFit" ,
            "Spfc_minDelta" ,
            "High_min_fit" ,
            "High_min_t" ,
            "High_max_se" ,
            "High_min_reads" ,
            "High_min_gMean" ,
            "High_max_below" ,
            "High_min_strains" ,
            "High_min_strain_fraction" ,
            "Strong_lr",
            "Strong_t"]:
        if not ( isinstance(params[x], int) or isinstance(params[x], float)):
            raise Exception(f"Key {x} should be either an int or a float. "
                        f"Instead it is {type(params[x])}")

    logging.info("Done validating all params - behavior as expected.")
    
    return params


def create_barseqPy_config(params, cfg_fp):
    """
    Args:
        params (d): The input to the entire program
        cfg_fp (str): Filepath to where this config dict will be written,
                      replaces existing config dict.
    Description:
        We create the config dict for BarSeqPy (2nd general part of analysis),
        which is normally to be placed in lib/BarSeqPy/config.json
        The way it should look is in lib/BarSeqPy/config.json before replacement.
    """

    nTopCofit = None if params["nTopCofit"] == -1 else params["nTopCofit"]

    logging.info(f"Preparing barseqPy to be written at {cfg_fp}.")
    barseqPy_cfg = {
        "nDebug_cols": None,
        "starting_debug_col": 0,
        "dp1_cfg": {
            "drop_exps": params["drop_exps"],
            "okControls": params["okControls"]
        },
        "dp2_cfg": {
            "minSampleReads" : params["minSampleReads"],
            "minGenesPerScaffold" : params["minGenesPerScaffold"],
            "minT0Strain" : params["minT0Strain"],
            "minT0Gene" : params["minT0Gene"],
            "minGenesAllowed": params["minGenesAllowed"],
            "minGenesUsed12" : params["minGenesUsed12"],
            "okControls": params["okControls"],
            "okDay" : params["okDay"],
            "okLane" : params["okLane"]
        },
        "an1_cfg": {
            "minGenesPerScaffold": params["minGenesPerScaffold"],
            "base_se": 0.1,
            "norm_median_window": params["norm_median_window"],
            "avgstrn": {
                "minGeneFactorNStrains": 3, 
                "strainFitAdjust": 0,
                "maxWeight": 20 
            }
            
        },
        "an2_cfg": {
            "minT0Strain": params["minT0Strain"],
            "status_d": {
                "min_gMed" : params["min_gMed"],
                "max_mad12" : params["max_mad12"],
                "min_cor12" : params["min_cor12"],
                "max_gccor" : params["max_gccor"],
                "max_adjcor" : params["max_adjcor"]
            }
        },
        "an3_cfg": {
            "compute_cofit_bool" : params["compute_cofit_bool"],
            "compute_High_bool" : params["compute_High_bool"],
            "compute_spfc_bool": params["compute_spfc_bool"],
            "nTopCofit" : nTopCofit,
            "minCofitExp" : params["minCofitExp"],
            "spec_cfg": {
                "minT" : params["Spfc_minT"],
                "minFit" : params["Spfc_minFit"],
                "percentile" : params["Spfc_percentile"],
                "percentileFit" : params["Spfc_percentileFit"],
                "minDelta" : params["Spfc_minDelta"]
            },
            "high_cfg": {
                "min_fit" : params["High_min_fit"],
                "min_t" : params["High_min_t"],
                "max_se" : params["High_max_se"],
                "min_reads" : params["High_min_reads"],
                "min_gMean" : params["High_min_gMean"],
                "max_below" : params["High_max_below"],
                "min_strains" : params["High_min_strains"],
                "min_strain_fraction" : params["High_min_strain_fraction"]
            }
        },
        "fst_cfg": {
            "strong_lr": params["Strong_lr"],
            "strong_t": params["Strong_t"]
        }
    }

    with open(cfg_fp, "w") as g:
        g.write(json.dumps(barseqPy_cfg, indent=2))

    logging.info(f"Wrote barseqPy config at {cfg_fp}.")

    return None



# op_name string, (output_name)
def check_output_name(op_name):
    op_name = op_name.replace(' ', '_')
    rgx = re.search(r'[^\w]', op_name)
    if rgx:
        logging.warning("Non-alphanumeric character in output name: " + rgx[0])
        op_name = "Default_Name_Check_Chars"
    return op_name
