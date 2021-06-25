#!python3
import re
import logging
import os


def validate_params(params):
    """What do we expect params to have as keys:
        genome_ref
        poolfile_ref
        exps_ref
        sets_ref
        output_name
        workspace_name

    Returns: ^ with strings as values for all
    """

    for p in ["genome_ref", "poolfile_ref", "exps_ref"]:
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


# op_name string, (output_name)
def check_output_name(op_name):
    op_name = op_name.replace(' ', '_')
    rgx = re.search(r'[^\w]', op_name)
    if rgx:
        logging.warning("Non-alphanumeric character in output name: " + rgx[0])
        op_name = "Default_Name_Check_Chars"
    return op_name
