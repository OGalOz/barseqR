#!python3

import logging
import os
import sys
import shutil
from Bio import SeqIO
#from Util.genbank_to_gene_table import genbank_and_genome_fna_to_gene_table 

# Both vp and d_d are dicts - validate_params and download_dict, respectively.
def download_files(vp, d_d):
    """
    In this function we download files using the data file util,
    the workspace, and the genome file util.
    genome_to_genbank spec:
        https://github.com/kbaseapps/GenomeFileUtil/blob/master/GenomeFileUtil.spec
    dfu.get_objects spec:
        https://github.com/kbaseapps/DataFileUtil/blob/master/DataFileUtil.spec
    
    Args:
        vp: validated params (d) (all values str)
            genome_ref
            mutantpool_ref
            exps_ref
            sets_ref
            output_name
            workspace_name
            
        d_d: download dict (d), keys:
           "dfu": dfu, datafile
           "gfu": gfu, genomefile
            "gt_obj": GeneTable Object 
           "ws": ws, workspace
           "smpl_s": sample service util
           "sets_dir": sets_dir,
           "mutantpool_path": mutantpool_path,
           "gene_table_fp": gene_table_fp,
           "exps_file": exps_file
           "scratch_dir": path to scratch directory

    """
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
    logging.info(vp['genome_ref'])

    ggo_list = get_genome_organism_name(ws, vp['genome_ref']) 
    organism_name = ggo_list[0]
    my_results = ggo_list[1]

    logging.info(my_results)


    GenomeToGenbankResult = gfu.genome_to_genbank({
                                        'genome_ref': vp['genome_ref']
                            })


    genbank_fp = GenomeToGenbankResult['genbank_file']['file_path']

    genome_fna_fp = get_fa_from_scratch(d_d["scratch_dir"])

    if genome_fna_fp is None:
        raise Exception("GFU Genome To Genbank did not download Assembly file in expected Manner.")

    # Download pool file and get related info. Name it pool.n10, place in indir
    # Ensure related to genome through organism_name
    res = download_mutantpool(vp['mutantpool_ref'], d_d['mutantpool_path'], dfu)
    mutantpool_path, related_genome_name, related_genome_ref = res

    '''
    if not (related_genome_name == organism_name):
        raise Exception("mutantpool organism name does not match genome " \
                + "organism name")
    '''

    # Note that this gene table will be created at workdir/g2gt_results/genes.GC
    g2gt_results = d_d['gt_obj'].genome_to_genetable({'genome_ref': vp['genome_ref']})
    logging.info(g2gt_results)
    gene_table_fp = os.path.join(d_d['scratch_dir'], 'g2gt_results', 'genes.GC')
    shutil.copy(gene_table_fp, os.path.join(d_d['scratch_dir'], 'indir', 'genes.GC'))
    '''
    # Convert genbank file to genes table, name it indir/genes.GC
    gt_cfg_dict = get_gene_table_config_dict(genbank_fp)
    genbank_and_genome_fna_to_gene_table(genbank_fp, 
                                        genome_fna_fp, 
                                        d_d['gene_table_fp'])
    '''

    #convert_genbank_to_gene_table(genbank_fp, d_d['gene_table_fp'], gt_cfg_dict)


    # Download Experiments File, name it FEBA_Barseq.TSV, place in scratch/indir
    #download_sample_set_to_file(d_d["smpl_s"], d_d['exps_file'], vp['exps_ref'], dfu)
    download_exps_file(dfu, d_d['exps_file'], vp['exps_ref'])


    # Download Set Files
    # set_names_list just contains the names of the sets without extensions
    # set_fps_list is a list of set filepaths
    set_names_list, set_fps_list = download_sets_from_refs(vp['sets_refs'], dfu,
                                                            organism_name,
                                                            sets_dir)

    
    
    DownloadResults = {
            "org" : organism_name,
            "set_names_list": set_names_list,
            "set_fps_list": set_fps_list
    }

    return DownloadResults


# Gets mutantpool path
def download_mutantpool(mutantpool_ref, mutantpool_path, dfu):

    GetObjectsParams = {
            'object_refs': [mutantpool_ref]
            }

    # We get the handle id
    mutantpoolObjectData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
    logging.info("DFU Pool File Get objects results:")
    logging.info(mutantpoolObjectData)

        
    related_genome_name = mutantpoolObjectData['related_organism_scientific_name']
    related_genome_ref = mutantpoolObjectData['related_genome_ref']

    mutantpool_handle = mutantpoolObjectData['mutantpool']


    # Set params for shock to file
    ShockToFileParams = {
            "handle_id": mutantpool_handle,
            "file_path": mutantpool_path,
            "unpack": "uncompress"
            }
    ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
    logging.info(ShockToFileOutput)
    # mutantpool is at location "mutantpool_path"

    return [mutantpool_path, related_genome_name, related_genome_ref]


# We want scaffold_name and description_name
def get_gene_table_config_dict(genbank_fp):

    record_generator = SeqIO.parse(open(genbank_fp), "genbank") 

    # the first record is the entire scaffold (?)
    record = next(record_generator)
    
    #logging.info(record)

    logging.info("Genbank Description: {}".format(record.description))
    logging.info("Genbank Scaffold Name: {}".format(record.id))
    my_id = record.id
    my_desc = record.description

    gene_table_config_dict = {
            "scaffold_name": my_id,
            "description": my_desc 
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
    return [scientific_name, res]

# Here we download the "poolcount" files, known as sets
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
            raise Exception("mutantpool organism name does not match genome " \
                + "organism name")
         

        barcodecount_handle = SetInfo['barcodecount']
        barcodecount_fn = SetInfo['set_name'] + ".poolcount"
        set_names_list.append(SetInfo['set_name'])
        barcodecount_fp = os.path.join(sets_dir, barcodecount_fn)

        
        # Set params for shock to file
        ShockToFileParams = {
                "handle_id": barcodecount_handle,
                "file_path": barcodecount_fp,
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
    """We download an experiments file

    Args:
        dfu: DataFileUtil class object
        exps_fp: (str) Path to download exps file to
        exps_ref: (str) Reference to file
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


# DEPRECATED
def download_sample_set_to_file(smpl_s, smpl_set_fp, exps_ref, dfu):
    """We download an experiments file as sample set 

    Args:
        smpl_s: sample service util
        exps_fp: (str) Path to download sample set file to:
                    indir/FEBA_BarSeq.tsv
        exps_ref: (str) Reference to file
        dfu: DataFileUtil class object
    """

    GetObjectsParams = {
            'object_refs': [exps_ref]
    }

    # We get the handle id
    sampleSetRes = dfu.get_objects(GetObjectsParams)
    logging.info(sampleSetRes)
    sampleSetObjectData = sampleSetRes['data'][0]['data']
    logging.info("DFU Sample Set Get objects results:")
    logging.info(sampleSetObjectData)
    samples = sampleSetObjectData["samples"]
    #samples is a list of dicts with 'id', 'name', 'version'
    for smpl in samples:
        sm_ret = smpl_s.get_sample({
            "id": smpl['id'],
            "version": smpl['version'],
            "as_admin": False
            })
        logging.info(sm_ret)


    raise Exception("242")
    # This will raise an error:
    smpl_stfile_handle = sampleSetObjectData['expsfile']


    # Set params for shock to file
    ShockToFileParams = {
            "handle_id": smpl_stfile_handle,
            "file_path": smpl_set_fp,
            "unpack": "uncompress"
            }
    ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
    logging.info(ShockToFileOutput)
    smpl_st_fp = ShockToFileOutput['file_path']
    # expsfile is at location "expsfile_path"

    return smpl_st_fp


def get_fa_from_scratch(scratch_dir):
    """
    Careful... May not work in the Future
    Inputs:
        scratch_dir: (str) Path to work dir/ tmp etc..
    Outputs:
        FNA fp: (str) Automatic download through GenbankToGenome
    """
    
    fna_fp = None
    scratch_files = os.listdir(scratch_dir)
    for f in scratch_files:
        if f[-2:] == "fa":
            fna_fp = os.path.join(scratch_dir, f)
            break
    
    if fna_fp is None:
        logging.warning("Could not find Assembly FNA file in scratch (work) dir")

    return fna_fp
