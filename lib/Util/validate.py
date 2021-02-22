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
    
    return params


# op_name string, (output_name)
def check_output_name(op_name):
    op_name = op_name.replace(' ', '_')
    rgx = re.search(r'[^\w]', op_name)
    if rgx:
        logging.warning("Non-alphanumeric character in output name: " + rgx[0])
        op_name = "Default_Name_Check_Chars"
    return op_name
