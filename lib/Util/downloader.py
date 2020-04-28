#!python3

import logging
import os
from Bio import SeqIO
from Util.genbank_to_gene_table import convert_genbank_to_gene_table

# Both vp and d_d are dicts - validate_params and download_dict, respectively.
def download_files(vp, d_d):
    logging.info("DOWNLOADING FILES-------------------")

    # Data File Util Client:
    dfu = d_d['dfu']
    # Genome File Util Client:
    gfu = d_d['gfu']
    # Workspace Client:
    ws = d_d['ws']
    # Where to download the sets to
    sets_dir = d_d['sets_dir']
    
    # Download genbank file, get organism name 
    GenomeToGenbankResult = gfu.genome_to_genbank({
                                                'genome_ref': vp['genome_ref']
                                                })

    genbank_fp = GenomeToGenbankResult['genbank_file']['file_path']

    organism_name = get_genome_organism_name(ws, vp['genome_ref']) 

    # Download pool file and get related info. Name it pool.n10, place in indir
    # Ensure related to genome through organism_name
    res = download_poolfile(vp['poolfile_ref'], d_d['poolfile_path'], dfu)
    poolfile_path, related_genome_name, related_genome_ref = res

    if not (related_genome_name == organism_name):
        raise Exception("Poolfile organism name does not match genome " \
                + "organism name")

    # Convert genbank file to genes table, name it indir/genes.GC
    gt_cfg_dict = get_gene_table_config_dict(genbank_fp)
    convert_genbank_to_gene_table(genbank_fp, d_d['gene_table_fp'], gt_cfg_dict)


    # Download Experiments File, name it FEBA_Barseq.TSV, place in indir
    download_exps_file(dfu, d_d['exps_file'], vp['exps_ref'])

    # Download Set Files
    set_names_list, set_fps_list = download_sets_from_refs(vp['sets_refs'], dfu, organism_name,
                                                            sets_dir)

    
    
    DownloadResults = {
            "org" : organism_name,
            "set_names_list": set_names_list,
            "set_fps_list": set_fps_list
            }
    return DownloadResults


# Gets poolfile path
def download_poolfile(poolfile_ref, poolfile_path, dfu):

    GetObjectsParams = {
            'object_refs': [poolfile_ref]
            }

    # We get the handle id
    PoolFileObjectData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
    logging.info("DFU Pool File Get objects results:")
    logging.info(PoolFileObjectData)

        
    related_genome_name = PoolFileObjectData['related_organism_scientific_name']
    related_genome_ref = PoolFileObjectData['related_genome_ref']

    poolfile_handle = PoolFileObjectData['poolfile']


    # Set params for shock to file
    ShockToFileParams = {
            "handle_id": poolfile_handle,
            "file_path": poolfile_path,
            "unpack": "uncompress"
            }
    ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
    logging.info(ShockToFileOutput)
    # Poolfile is at location "poolfile_path"

    return [poolfile_path, related_genome_name, related_genome_ref]


# We want scaffold_name and description_name
def get_gene_table_config_dict(genbank_fp):

    record = SeqIO.read(genbank_fp, "genbank") 

    logging.info("Genbank Description: {}".format(record.description))
    logging.info("Genbank Scaffold Name: {}".format(record.id))
    
    gene_table_config_dict = {
            "scaffold_name": record.id,
            "description": record.description
            }

    return gene_table_config_dict


def get_genome_organism_name(ws, genome_ref):
    # Getting the organism name using WorkspaceClient
    res = ws.get_objects2(
        {
            "objects": [
                {
                    "ref": genome_ref,
                    "included": ["scientific_name"],
                }
            ]
        }
    )
    scientific_name = res["data"][0]["data"]["scientific_name"]
    return scientific_name

def download_sets_from_refs(ref_list, dfu, organism_name, sets_dir):

    logging.info("Downloading .poolcount files")

    GetObjectsParams = {
            'object_refs': ref_list
    }

    SetsInfoList = dfu.get_objects(GetObjectsParams)['data']
    logging.info(SetsInfoList)


    set_names_list = []
    set_fps_list = []
    for obj in SetsInfoList:
        SetInfo = obj['data']

        if SetInfo['related_organism_scientific_name'] != organism_name:
            raise Exception("Poolfile organism name does not match genome " \
                + "organism name")
         

        setfile_handle = SetInfo['poolcount']
        setfile_fn = SetInfo['set_name'] + ".poolcount"
        set_names_list.append(SetInfo['set_name'])
        setfile_fp = os.path.join(sets_dir, setfile_fn)

        
        # Set params for shock to file
        ShockToFileParams = {
                "handle_id": setfile_handle,
                "file_path": setfile_fp,
                "unpack": "uncompress"
                }
        ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
        logging.info(ShockToFileOutput)
        set_fp = ShockToFileOutput['file_path']
        set_fps_list.append(set_fp)

    logging.info("Sets Names list: ")
    logging.info(set_names_list)

    logging.info("Sets Filepaths List: ")
    logging.info(set_fps_list)

    return [set_names_list, set_fps_list]

def download_exps_file(dfu, exps_fp, exps_ref):
    """
    exps_fp is path to experiments file final loc
    """

    GetObjectsParams = {
            'object_refs': [exps_ref]
            }

    # We get the handle id
    expsFileObjectData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
    logging.info("DFU Experiment File Get objects results:")
    logging.info(expsFileObjectData)

    expsfile_handle = expsFileObjectData['expsfile']


    # Set params for shock to file
    ShockToFileParams = {
            "handle_id": expsfile_handle,
            "file_path": exps_fp,
            "unpack": "uncompress"
            }
    ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
    logging.info(ShockToFileOutput)
    exps_fp = ShockToFileOutput['file_path']
    # expsfile is at location "expsfile_path"

    return exps_fp
    
