# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import shutil
import requests

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from kb_gffupload.kb_gffuploadImpl import kb_gffupload
from kb_gffupload.kb_gffuploadServer import MethodContext


class kb_gffuploadTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        user_id = requests.post(
            'https://kbase.us/services/authorization/Sessions/Login',
            data='token={}&fields=user_id'.format(token)).json()['user_id']
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_gffupload',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_gffupload'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_gffupload(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_gffupload_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_simple_upload(self):
        gff_upload = self.getImpl()

        data_dir="data/"
        scratch_data_dir = os.path.join(self.scratch, data_dir)
        shutil.copytree(data_dir, scratch_data_dir)

        fasta_file = "Test_v1.0.fa"
        gff_file = "Test_v1.0.gene.gff3"

        fasta_path = scratch_data_dir+"/"+fasta_file
        gff_path = scratch_data_dir+"/"+gff_file

        shutil.copy(data_dir+"/"+fasta_file, fasta_path)
        shutil.copy(data_dir+"/"+gff_file, gff_path)

        ws_obj_name = 'MyGenome'
        ws_name = self.getWsName()
        scientific_name = "Populus trichocarpa"

        ### Test for a Local Function Call
        print('attempting upload via local function directly')

        result = gff_upload.fasta_gff_to_genome(self.getContext(), 
            {
                'fasta_file' : fasta_path,
                'gff_file' : gff_path,
                'workspace_name':ws_name,
                'genome_name':ws_obj_name,
                'scientific_name':scientific_name
            })[0]
        pprint(result)
        self.assertIsNotNone(result['genome_ref'])

        shutil.rmtree(scratch_data_dir)

        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        pass
