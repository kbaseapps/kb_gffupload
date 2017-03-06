/*
A KBase module: kb_gffupload
*/

module kb_gffupload {

    typedef structure {
        string path;
        string shock_id;
        string ftp_url;
    } File;

    typedef mapping<string, string> usermeta;

    typedef structure {
        string genome_ref;
    } GenomeSaveResult;

    /* 
        genome_name - becomes the name of the object
	workspace_name - the name of the workspace it gets saved to.
	source - Source of the file typically something like RefSeq or Ensembl
	taxon_ws_name - where the reference taxons are : ReferenceTaxons
	taxon_reference - if defined, will try to link the Genome to the specified
        taxonomy object insteas of performing the lookup during upload
	release - Release or version number of the data 
          per example Ensembl has numbered releases of all their data: Release 31
	generate_ids_if_needed - If field used for feature id is not there, 
          generate ids (default behavior is raising an exception)
        genetic_code - Genetic code of organism. Overwrites determined GC from 
          taxon object
	type - Reference, Representative or User upload
    */
    typedef structure {
        File fasta_file;
	File gff_file;

        string genome_name;
        string workspace_name;

        string source;
        string taxon_wsname;
        string taxon_reference;
	string release;

	int    genetic_code;
	string type;
	usermeta metadata;
    } FastaGFFToGenomeParams;

    funcdef fasta_gff_to_genome(FastaGFFToGenomeParams params)
                returns (GenomeSaveResult result) authentication required;

};
