import os
import sys
import shutil
import traceback
import uuid
import urllib2

from urlparse import urlparse
from pprint import pprint, pformat

from biokbase.workspace.client import Workspace
from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport

from kb_gffupload.FastaGFFUploaderScript import upload_genome

class FastaGFFToGenome:

    def __init__(self, sdk_config):
        self.cfg = sdk_config

    def import_file(self, ctx, params):

        # 1) validate parameters and extract defaults
        self.validate_params(params)

        # 2) define default parameters
        default_params = {
            'taxon_wsname': 'ReferenceTaxons',
            'scientific_name': "unknown_taxon",
            'taxon_reference': None,
            'source':'User',
            'release': None,
            'type': 'User upload'
        }

        # 3) Add defaults if they don't exist
        for field in default_params:
            if field not in params:
                params[field] = default_params[field]

        # 4) Do the upload
        result = upload_genome(
                shock_service_url     = self.cfg['shock-url'],
                handle_service_url    = self.cfg['handle-service-url'],
                workspace_service_url = self.cfg['workspace-url'],
                callback_url = os.environ['SDK_CALLBACK_URL'],

                input_fasta_file = params["fasta_file"],
                input_gff_file   = params["gff_file"],
        
                workspace_name   = params['workspace_name'],
                core_genome_name = params['genome_name'],
                scientific_name  = params['scientific_name'],

                taxon_wsname    = params['taxon_wsname'],
                taxon_reference = params['taxon_reference'],

                source          = params['source'],
                release         = params['release'],
                genome_type     = params['type']
            )

        # 4) Generate Report
        output_data_ref = params['workspace_name']+"/"+params['genome_name']
        reportObj = {'objects_created':[{'ref':output_data_ref, 'description':'KBase Genome object'}],
                     'text_message':result['report_string']}

        reportClient = KBaseReport(os.environ['SDK_CALLBACK_URL'])
        report_info = reportClient.create({'report':reportObj, 'workspace_name':params['workspace_name']})

        # 5) return the result
        info = result['genome_info']
        details = {
            'genome_ref':  str(info[6]) + '/' + str(info[0]) + '/' + str(info[4]),
            'report_name': report_info['report_name'],
            'report_ref':  report_info['report_ref']
        }

        return details

    def validate_params(self, params):
        for key in ('workspace_name','genome_name','fasta_file','gff_file'):
            if key not in params:
                raise ValueError('required "'+key+'" field was not defined')

        # files must always be paths
        for key in ('fasta_file','gff_file'):
            file = params[key]
            if not isinstance(file, basestring):
                raise ValueError('required file "'+file+'" for "'+key+'" field must be a string')
            if not os.path.isfile(file):
                raise ValueError('required file "'+file+'" for "'+key+'" is not a valid file')

        valid_types = ['Reference','User upload','Representative']
        if 'type' in params and params['type'] not in valid_types:
            raise ValueError('Entered value for type is not one of the valid entries of "Reference", "Representative" or "User upload"')
