# python3
# Includes read_table, read_column_names

import re


# Following function takes a filename and a list or required fields (from TSV)
# Returns a list of dicts with required field pointing to value
def read_table(fp, required):
    with open(fp, "r") as f:
        file_str = f.read()
    file_list = file_str.split("\n")
    header_line = file_list[0]
    # Check for Mac Style Files
    if re.search(r"\t", header_line) and re.search(r"\r", header_line):
        raise Exception(
            (
                "Tab-delimited input file {} is a Mac-style text file "
                "which is not supported.\n"
                "Use\ndos2unix -c mac {}\n to convert it to a Unix "
                "text file.\n"
            ).format(fp, fp)
        )
    cols = header_line.split("\t")
    cols_dict = {}
    for i in range(len(cols)):
        cols_dict[cols[i]] = i
    for field in required:
        if field not in cols_dict:
            raise Exception("No field {} in {}".format(field, fp))
    rows = []
    for i in range (1, len(file_list)):
        line = file_list[i]
        # if last line empty
        if len(line) == 0:
            continue
        line = re.sub(r'[\r\n]+$', '', line)
        split_line = line.split("\t")
        if not len(split_line) == len(cols):
            raise Exception("Wrong number of columns in:\n{}\nin {} l:{}".format(
                line, fp, i))
        new_dict = {}
        for i in range(len(cols)):
            new_dict[cols[i]] = split_line[i]
        rows.append(new_dict)

    return rows


def read_column_names(fp):
    with open(fp, "r") as f:
        f_line = f.read().split("\n")[0]
    split_line = f_line.rstrip().split("\t")
    return split_line
