# -*- coding: utf-8 -*-
#BEGIN_HEADER
from pprint import pprint

from kb_gffupload.FastaGFFToGenome import FastaGFFToGenome
#END_HEADER

class kb_gffupload:
    '''
    Module Name:
    kb_gffupload

    Module Description:
    A KBase module: kb_gffupload
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
        self.cfg=config
        #END_CONSTRUCTOR
        pass


    def fasta_gff_to_genome(self, ctx, params):
        """
        :param params: instance of type "FastaGFFToGenomeParams" (genome_name
           - becomes the name of the object workspace_name - the name of the
           workspace it gets saved to. source - Source of the file typically
           something like RefSeq or Ensembl taxon_ws_name - where the
           reference taxons are : ReferenceTaxons taxon_reference - if
           defined, will try to link the Genome to the specified taxonomy
           object insteas of performing the lookup during upload release -
           Release or version number of the data per example Ensembl has
           numbered releases of all their data: Release 31
           generate_ids_if_needed - If field used for feature id is not
           there, generate ids (default behavior is raising an exception)
           genetic_code - Genetic code of organism. Overwrites determined GC
           from taxon object type - Reference, Representative or User upload)
           -> structure: parameter "fasta_file" of type "File" -> structure:
           parameter "path" of String, parameter "shock_id" of String,
           parameter "ftp_url" of String, parameter "gff_file" of type "File"
           -> structure: parameter "path" of String, parameter "shock_id" of
           String, parameter "ftp_url" of String, parameter "genome_name" of
           String, parameter "workspace_name" of String, parameter "source"
           of String, parameter "taxon_wsname" of String, parameter
           "taxon_reference" of String, parameter "release" of String,
           parameter "genetic_code" of Long, parameter "type" of String,
           parameter "metadata" of type "usermeta" -> mapping from String to
           String
        :returns: instance of type "GenomeSaveResult" -> structure: parameter
           "genome_ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN fasta_gff_to_genome

        print('fasta_gff_to_genome -- paramaters = ')
        pprint(params)

        importer = FastaGFFToGenome(self.cfg)
        result = importer.import_file(ctx, params)

        #END fasta_gff_to_genome

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method fasta_gff_to_genome return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
