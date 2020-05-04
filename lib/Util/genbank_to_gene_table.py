#python3
import logging
import os
import sys
import copy
from shutil import copyfile, move
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
   scaffoldId, begin, end, strand, desc, type, sysName
        'type' with type=1 for protein-coding genes,
    other fields:
        sysName, name, GC, nTA
"""
def convert_genbank_to_gene_table(genbank_filepath, output_filepath, config_dict):


    # Change output file string to writing to output file.
    f = open(output_filepath, "w")


    # This is the output file start line:
    output_file_header = "scaffoldId\tbegin\tend\tstrand\tdesc\t" \
                            + "locusId\ttype\tsysName\n"
    f.write(output_file_header)

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
                current_row += str(current_feat.qualifiers['locus_tag'][0]) + \
                        "\t"
                sysName = str(current_feat.qualifiers['locus_tag'][0])
            else:
                current_row += "Unknown_Locus_tag." + "\t"
                logging.critical("Could not find locus tag in current_feat")
                sysName = "Unknown_Sys_Name" 

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
            current_row += typ + "\t"

            # Getting sysName

            current_row += sysName + "\n"

            f.write(current_row) 

    except:
        logging.critical("Could not parse all features in genbank file.")
        raise Exception("Parsing genbank file into gene table failed")

    # Finish writing to gene_table file at output_filepath
    f.close()

    
    #We remove duplicate gene lines and remove the last new line symbol
    # f is the file handle that is still open
    gene_table_fp = unduplicate_gene_table(output_filepath)

    if "keep_types" in config_dict:
        types_to_keep = config_dict["keep_types"]
        gene_table_fp = keep_types_gene_table(gene_table_fp, 
                                                    types_to_keep)

    logging.info("Wrote Gene Table to " + output_filepath)

    return output_filepath



"""
This function removes duplicates from the gene table
 Note gene_table_fp does not change throughout function!
Inputs:
    gene_table_fp: (str) A filepath to the gene table file
Outputs:
    gene_table_fp: (str) A filepath to the gene table file
Process: 
    Compares the location of features and if they are the same removes
        one of the two.
"""
def unduplicate_gene_table(gene_table_fp):

    # We read the file
    f = open(gene_table_fp, "r")
    header_line = f.readline()
    if header_line == '':
        raise Exception("Gene Table creation process broken A.")


    #Then for each line we check if it's a duplicate or not.
    #We add the indices of duplicate lines and then remove the lines in reverse order.
    # 'loc' means begin to end in sequence
    newline = f.readline()
    splitLine = newline.split("\t")
    existing_loc = splitLine[1:3]; existing_typ = splitLine[6]
    previous_index = 1
    linecount = 1
    # We create a set that keeps track of duplicate_line_indices
    duplicate_line_indices = set()
    newline = f.readline()
    while newline != '':
        linecount += 1
        splitLine = newline.split("\t")
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
                previous_index = linecount 
            else:
                if existing_typ == '1':
                    duplicate_line_indices.add(linecount)
                else:
                    duplicate_line_indices.add(previous_index)

        else:
            existing_loc = current_loc
            previous_index = linecount
        newline = f.readline()
    
    f.close()

    logging.info("Read {} lines from {}".format(linecount, gene_table_fp))

    logging.info("Number of duplicate lines: {}".format(len(
        duplicate_line_indices)))

    # Sorting indices so they ascend 
    duplicate_line_indices = sorted(list(duplicate_line_indices))

    # We create another file handle so we can copy only the good lines
    h = open(gene_table_fp, "r")
    # We write to a copy of the file
    g = open(gene_table_fp + ".copy", "w")

    newline = h.readline()
    linecount = 0
    newFileLength = 0
    while newline != '':
        if linecount not in duplicate_line_indices:
            g.write(newline)
            newFileLength += 1
        newline = h.readline()
        linecount += 1
    h.close()
    g.close()

    # We move the new file to the location of the old file
    move(gene_table_fp + ".copy", gene_table_fp)

    logging.info("New total line number (after duplicate line removal"
                "): {}".format(newFileLength))

    return gene_table_fp


"""
Inputs:
    gene_table_fp: (str) The gene table string
    types_to_keep: list<str> Each string in list is a type we want.
Outputs:
    gene_table_fp: (str) The gene table string.
"""
def keep_types_gene_table(gene_table_fp, types_to_keep):

    k = open(gene_table_fp, "r")
    header_line = k.readline()
    if header_line == '':
        raise Exception("Gene Table creation process broken B.")

    newline = k.readline()
    LineNum = 1
    non_good_type_indices = []
    while newline != '':
        current_type = newline.split("\t")[6]
        if current_type not in types_to_keep:
            non_good_type_indices.append(LineNum)
        newline = k.readline()
        LineNum += 1

    k.close()

    # We repeat the process from before with writing to copy

    # We create another file handle so we can copy only the good lines
    h = open(gene_table_fp, "r")
    # We write to a copy of the file
    g = open(gene_table_fp + ".copy", "w")

    newline = h.readline()
    linecount = 0
    newFileLength = 0
    while newline != '':
        if linecount not in non_good_type_indices:
            g.write(newline)
            newFileLength += 1
        newline = h.readline()
        linecount += 1
    h.close()
    g.close()

    # We move the new file to the location of the old file
    move(gene_table_fp + ".copy", gene_table_fp)

    logging.info("New total line number (after duplicate line removal"
                "): {}".format(newFileLength))

    return gene_table_fp


def test(args):
    logging.basicConfig(filename='ogTestLog.log', level=logging.DEBUG)
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
