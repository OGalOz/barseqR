#python3
import logging
import os
import sys
from Bio import SeqIO



"""
Inputs:
    genbank_filepath: (str) Path to genbank file.
    output_filepath: (str) Path to write genome table to
config_dict:
    scaffold_name: (str) The scaffold id in the genbank file
    description: (str) The description id in the genbank file
    OPTIONAL
    keep_types: (list) If you only want specific types
    

The genes table must include the fields
   scaffoldId, begin, end, strand, desc, type
        'type' with type=1 for protein-coding genes,
    other fields:
        sysName, name, GC, nTA
"""
def convert_genbank_to_gene_table(genbank_filepath, output_filepath, config_dict):

    # This is the output file start line:
    output_file_string = "scaffoldId\tbegin\tend\tstrand\tdesc\t" \
                            + "locusId\ttype\n"

    # We use BioPython SeqIO to parse genbank file:
    # https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
    gb_record = SeqIO.read(open(genbank_filepath, "r"), "genbank")

    genome_name = gb_record.name
    #Genome sequence:
    g_seq = gb_record.seq
    g_len = len(g_seq)

    #Genome features (list of features):
    g_features = gb_record.features
    g_feat_len = len(g_features)

    scaffoldId_exists= False
    if "scaffold_name" in config_dict:
        scaffold_id = config_dict["scaffold_name"]
        scaffoldId_exists = True
    try:
        for i in range(g_feat_len):
            current_row = ""
            current_feat = g_features[i]

            if scaffoldId_exists:
                if scaffold_id in g_features[i].qualifiers:
                    scaffold = g_features[i].qualifiers[scaffold_id]
                else:
                    """
                    logging.debug("Could not find scaffold id "
                            "{} in qualifiers:".format(scaffold_id))
                    logging.debug(g_features[i].qualifiers)
                    """
                    scaffold = "1"
            else:
                scaffold = "1"
            # ScaffoldId
            current_row += scaffold + "\t"
            # Begin
            current_row += str(current_feat.location.start) + "\t"
            # End
            current_row += str(current_feat.location.end) + "\t"
            # Strand
            if current_feat.strand == 1:
                strand = "+"
            elif current_feat.strand == -1:
                strand = "-"
            else:
                logging.critical("Could not recognize strand type.")
                raise Exception("Parsing strand failed.")
            current_row += strand + "\t"

            # Desc (Description)
            if "product" in current_feat.qualifiers.keys():
                current_row += str(current_feat.qualifiers['product'][0]) + "\t"
            else:
                current_row += "Unknown_Description." + "\t"
                logging.critical("Could not find protein in current_feat")

            # Locus ID:
            if "locus_tag" in current_feat.qualifiers.keys():
                current_row += str(current_feat.qualifiers['locus_tag'][0]) + "\t"
            else:
                current_row += "Unknown_Locus_tag." + "\t"
                logging.critical("Could not find locus tag in current_feat")

            # TYPE - Note that we don't like misc_feature or gene
            # May need to skip anything besides values under 10
            types_dict = {"CDS" : 1, "rRNA" : 2, "tRNA" : 5, 
                            "RNA" : 6, "transcript" : 6,
                            "pseudogene": 7, "misc_feature": 20, 
                            "gene": 21}
            typ_str = current_feat.type.strip()
            if typ_str in types_dict:
                typ = str(types_dict[typ_str])
            else:
                logging.info("Could not recognize type from feature: " \
                        + typ_str)
                typ = "0"
            current_row += typ + "\n"

            output_file_string += current_row

    except:
        logging.critical("Could not parse all features in genbank file.")
        raise Exception("Parsing genbank file into gene table failed")

    
    #We remove duplicate gene lines and remove the last new line symbol
    output_file_string = unduplicate_gene_table(output_file_string[:-2])

    if "keep_types" in config_dict:
        types_to_keep = config_dict["keep_types"]
        output_file_string = keep_types_gene_table(output_file_string, 
                                                    types_to_keep)
    with open(output_filepath, "w") as f:
        f.write(output_file_string)
    logging.info("Wrote Gene Table to " + output_filepath)

    return output_filepath



"""
This function removes duplicates from the gene table
Inputs:
    gene_table_string: (str) A string of the entire gene table file
Outputs:
    gene_table_string: (str) A string of the entire gene table file
Process: 
    Compares the location of features and if they are the same removes
        one of the two.
"""
def unduplicate_gene_table(gene_table_string):


    #first we split the gene_table into lines:
    split_list = gene_table_string.split("\n")
    header_line = split_list[0] + "\n" 
    gt_lines = split_list[1:]
    logging.info("Total number of lines besides headers: " + str(len(gt_lines)))

    #Then for each line we check if it's a duplicate or not.
    #We add the indices of duplicate lines and then remove the lines in reverse order.
    # 'loc' means begin to end in sequence
    splitLine = gt_lines[0].split("\t")
    existing_loc = splitLine[1:3]; existing_typ = splitLine[6]
    previous_index = 0
    # We create a set, duplicate_line_indices
    duplicate_line_indices = set()
    for i in range(1, len(gt_lines)):
        splitLine = gt_lines[i].split("\t")
        current_loc = splitLine[1:3]; crnt_typ = splitLine[6]
        if (current_loc[0] == existing_loc[0]) or (
                current_loc[1] == existing_loc[1]):
            if crnt_typ == '1':
                if existing_typ == '1':
                    logging.warning("Two overlapping location Protein "
                            "Features: loc: {}{},{}{} types: {},{}".format(
                                existing_loc[0], existing_loc[1],
                                current_loc[0], current_loc[1],
                                existing_typ, crnt_typ))
                duplicate_line_indices.add(previous_index)
                previous_index = i
            else:
                if existing_typ == '1':
                    duplicate_line_indices.add(i)
                else:
                    duplicate_line_indices.add(previous_index)

        else:
            existing_loc = current_loc
            previous_index = i
    logging.info("Duplicate Lines: " + str(len(duplicate_line_indices)))

    # Sorting indices so they ascend 
    duplicate_line_indices = sorted(list(duplicate_line_indices))
    #removing the indeces in reverse order:
    duplicate_line_indices.reverse()
    for i in range(len(duplicate_line_indices)):
        del gt_lines[duplicate_line_indices[i]]

    logging.info("New total line number (after duplicate line removal): " \
            + str(len(gt_lines)))



    #Converting list into string again:
    gene_table_string = header_line + "\n".join( gt_lines)


    return gene_table_string


"""
Inputs:
    gene_table_string: (str) The gene table string
    types_to_keep: list<str> Each string in list is a type we want.
Outputs:
    gene_table_string: (str) The gene table string.
"""
def keep_types_gene_table(gene_table_string, types_to_keep):

    split_list = gene_table_string.split("\n")
    header_line = split_list[0] + "\n" 
    gt_lines = split_list[1:]

    non_good_type_indices = []

    #For each line, we check if its type is 1. If not, we remove it later.
    for i in range(len(gt_lines)):
        current_type = gt_lines[i].split("\t")[6]
        if current_type not in types_to_keep:
            non_good_type_indices.append(i)

    #removing the indeces in reverse order:
    non_good_type_indices.reverse()
    for i in range(len(non_good_type_indices)):
        del gt_lines[non_good_type_indices[i]]


    logging.info("New total line number (after type 1): " + str(len(gt_lines)))

    #Converting list into string again:
    gene_table_string = header_line + "\n".join(gt_lines)

    return gene_table_string



def test(args):
    logging.basicConfig(level=logging.DEBUG)
    gb_fp = args[1]
    op_fp = args[2]
    config_dict = {"keep_types": ["1","5","6"]}
    convert_genbank_to_gene_table(gb_fp, 
                                    op_fp, 
                                    config_dict)

def main():
    """
    args should be genbank, output
    """
    args = sys.argv
    test(args)

if __name__ == "__main__":
    main()
