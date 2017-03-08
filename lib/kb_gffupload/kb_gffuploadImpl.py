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
    GIT_URL = "https://github.com/kbaseapps/kb_gffupload.git"
    GIT_COMMIT_HASH = "dbb1e6d3040aadc9e38c42e3d7c1fad42802a047"

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
        :param params: instance of type "FastaGFFToGenomeParams" (fasta_file
           - file containing assembled contigs/chromosomes gff_file - file
           containing gene models (_must_ contain 'gene', 'mRNA', and 'CDS')
           genome_name - becomes the name of the object workspace_name - the
           name of the workspace it gets saved to. source - Source of the
           file typically something like RefSeq or Ensembl taxon_ws_name -
           where the reference taxons are : ReferenceTaxons taxon_reference -
           if defined, will try to link the Genome to the specified taxonomy
           object release - Release or version number of the data
           (i.e.Ensembl Release 31 or Phytozome Release V11) type -
           Reference, Representative or User) -> structure: parameter
           "fasta_file" of type "File_Path", parameter "gff_file" of type
           "File_Path", parameter "genome_name" of String, parameter
           "workspace_name" of String, parameter "source" of String,
           parameter "taxon_wsname" of String, parameter "taxon_reference" of
           String, parameter "release" of String, parameter "type" of String
        :returns: instance of type "GenomeSaveResult" -> structure: parameter
           "genome_ref" of String, parameter "report_name" of String,
           parameter "report_ref" of String
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
