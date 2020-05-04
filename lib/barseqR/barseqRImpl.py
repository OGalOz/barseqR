# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from RunDir.run_barseqR import RunBarSeq
from Util.validate import validate_params
from Util.downloader import download_files
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.DataFileUtilClient import DataFileUtil
#END_HEADER


class barseqR:
    '''
    Module Name:
    barseqR

    Module Description:
    A KBase module: barseqR
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.ws_url = config['workspace-url']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_barseqR(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_barseqR

        # SETUP - Unrelated to inputs --------
        
        logging.basicConfig(level=logging.DEBUG)

        # We create important classes
        dfu = DataFileUtil(self.callback_url)
        gfu = GenomeFileUtil(self.callback_url)
        myToken = os.environ.get('KB_AUTH_TOKEN', None)
        ws = Workspace(self.ws_url, token=myToken)
        ws_id = ws.get_workspace_info({'workspace': params['workspace_name']})[0]
        

        # We create indir, outdir, sets_dir (Input, Output, Sets)
        indir = os.path.join(self.shared_folder, "indir")
        os.mkdir(indir)

        outdir = os.path.join(self.shared_folder, "outdir")
        os.mkdir(outdir)

        sets_dir = os.path.join(indir, "sets_dir")
        os.mkdir(sets_dir)

        metadir = '/kb/module/lib/RunDir/metadata'
        if not (os.path.isdir(metadir)):
            raise Exception("metadata directory not found at: {}".format(metadir))


        # We prepare locations of input files

        poolfile_path = os.path.join(indir, "pool.n10")
        gene_table_fp = os.path.join(indir, "genes.GC")
        exps_file = os.path.join(indir, "FEBA_Barseq.TSV")

        # END SETUP


        # VALIDATE PARAMS:
        # Needs 'poolfile_ref' as poolfile ref
        logging.info("PARAMS:")
        logging.info(params)
        val_par = validate_params(params)
        val_par['username'] = ctx['user_id']


        # DOWNLOAD FILES
        download_dict = {
                "dfu": dfu,
                "gfu": gfu,
                "ws": ws,
                "sets_dir": sets_dir,
                "poolfile_path": poolfile_path,
                "gene_table_fp": gene_table_fp,
                "exps_file": exps_file
                }
        # We copy input files to proper directories.
        # vp must contain genome_ref, poolfile_ref, exps_ref, sets_refs (list)
        # DownloadResults must contain keys 'org', 'set_names_list', 'set_fps_list'
        # set_names_list value contains the names of the sets without extensions
        DownloadResults = download_files(val_par, download_dict)





        # INDIR:
        # Indir contains Poolcount Files, FEBA_Barseq.tsv, genes.GC, pool.n10
        # Outdir must be a directory.
        #
        #

        # Get args in this format:
        # [-org, org_name, -indir, Scratch_Dir_Input, -metadir, Fixed meta dir,
        # -outdir, scratch_dir_output, -sets_dir, within scratch_dir_input, 
        # -sets, set1 (sets_dir), set2 (sets_dir), set3 (sets_dir), ... ]
        # Note meta dir is called metadata and is in RunDir

        # Running the entire program:
        arg_list = ["-org", DownloadResults['org'], '-indir', indir,
                '-metadir', metadir, '-outdir', outdir, 
                '-sets_dir', sets_dir, '-sets']
        arg_list += DownloadResults['set_names_list']
        RunBarSeq(arg_list)


        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['parameter_1']},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_barseqR

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_barseqR return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
