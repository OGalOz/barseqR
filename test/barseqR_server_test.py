# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from barseqR.barseqRImpl import barseqR
from barseqR.barseqRServer import MethodContext
from barseqR.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class barseqRTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('barseqR'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'barseqR',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = barseqR(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_BARSEQR_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')




    def test_ci_1(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        # Needs genome_ref, poolfile_ref, exps_ref, sets_refs (list ), 
        # output_name

        # Burk 376
        genome_ref = "58816/3/1"
        exps_ref = "58816/9/1"

        poolfile_ref = "49371/19/1"
        sets_refs = ["49371/21/1"]
        output_name = "Test_1"

        ret = self.serviceImpl.run_barseqR(self.ctx, 
                                            {
                                             'workspace_name': self.wsName,
                                             'genome_ref': genome_ref,
                                             'poolfile_ref': poolfile_ref,
                                             'exps_ref': exps_ref,
                                             'sets_refs': sets_refs,
                                             'output_name': output_name
                                             }

                                           )
    """
    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_appdev_1(self):
        
        genome_ref = "52416/6/1"
        # actual new GFF + genome genome_ref = "52416/3/1"
        poolfile_ref = "49371/19/1"
        exps_ref = "58816/5/1"
        sets_refs = ["49371/21/1"]
        output_name = "Test_1"

        ret = self.serviceImpl.run_barseqR(self.ctx, 
                                            {
                                             'workspace_name': self.wsName,
                                             'genome_ref': genome_ref,
                                             'poolfile_ref': poolfile_ref,
                                             'exps_ref': exps_ref,
                                             'sets_refs': sets_refs,
                                             'output_name': output_name
                                             }
                                           )
    """
                                                             
