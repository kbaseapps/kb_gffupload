/*
A KBase module: kb_gffupload
*/

module kb_gffupload {

    typedef string File_Path;

    typedef structure {
        string genome_ref;
	string report_name;
	string report_ref;
    } GenomeSaveResult;

    /* 
    fasta_file - file containing assembled contigs/chromosomes
    gff_file - file containing gene models (_must_ contain 'gene', 'mRNA', and 'CDS')
    genome_name - becomes the name of the object
    workspace_name - the name of the workspace it gets saved to.
    source - Source of the file typically something like RefSeq or Ensembl
    taxon_wsname - where the reference taxons are : ReferenceTaxons
    scientific_name - Full name of organism, including genus, to be used with taxon lookup
    taxon_reference - if defined, will try to link the Genome to the specified taxonomy object
    release - Release or version number of the data (i.e.Ensembl Release 31 or Phytozome Release V11)
    type - Reference, Representative or User upload
    */
    typedef structure {
        File_Path fasta_file;
	File_Path gff_file;

        string genome_name;
        string workspace_name;

        string taxon_wsname;
	string scientific_name;
        string taxon_reference;

        string source;
	string release;
	string type;
    } FastaGFFToGenomeParams;

    funcdef fasta_gff_to_genome(FastaGFFToGenomeParams params)
                returns (GenomeSaveResult result) authentication required;

};
