
package us.kbase.kbgffupload;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: FastaGFFToGenomeParams</p>
 * <pre>
 * fasta_file - file containing assembled contigs/chromosomes
 * gff_file - file containing gene models (_must_ contain 'gene', 'mRNA', and 'CDS')
 * genome_name - becomes the name of the object
 * workspace_name - the name of the workspace it gets saved to.
 * source - Source of the file typically something like RefSeq or Ensembl
 * taxon_ws_name - where the reference taxons are : ReferenceTaxons
 * taxon_reference - if defined, will try to link the Genome to the specified taxonomy object
 * release - Release or version number of the data (i.e.Ensembl Release 31 or Phytozome Release V11)
 * type - Reference, Representative or User
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fasta_file",
    "gff_file",
    "genome_name",
    "workspace_name",
    "source",
    "taxon_wsname",
    "taxon_reference",
    "release",
    "type"
})
public class FastaGFFToGenomeParams {

    @JsonProperty("fasta_file")
    private String fastaFile;
    @JsonProperty("gff_file")
    private String gffFile;
    @JsonProperty("genome_name")
    private String genomeName;
    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("source")
    private String source;
    @JsonProperty("taxon_wsname")
    private String taxonWsname;
    @JsonProperty("taxon_reference")
    private String taxonReference;
    @JsonProperty("release")
    private String release;
    @JsonProperty("type")
    private String type;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("fasta_file")
    public String getFastaFile() {
        return fastaFile;
    }

    @JsonProperty("fasta_file")
    public void setFastaFile(String fastaFile) {
        this.fastaFile = fastaFile;
    }

    public FastaGFFToGenomeParams withFastaFile(String fastaFile) {
        this.fastaFile = fastaFile;
        return this;
    }

    @JsonProperty("gff_file")
    public String getGffFile() {
        return gffFile;
    }

    @JsonProperty("gff_file")
    public void setGffFile(String gffFile) {
        this.gffFile = gffFile;
    }

    public FastaGFFToGenomeParams withGffFile(String gffFile) {
        this.gffFile = gffFile;
        return this;
    }

    @JsonProperty("genome_name")
    public String getGenomeName() {
        return genomeName;
    }

    @JsonProperty("genome_name")
    public void setGenomeName(String genomeName) {
        this.genomeName = genomeName;
    }

    public FastaGFFToGenomeParams withGenomeName(String genomeName) {
        this.genomeName = genomeName;
        return this;
    }

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public FastaGFFToGenomeParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("source")
    public String getSource() {
        return source;
    }

    @JsonProperty("source")
    public void setSource(String source) {
        this.source = source;
    }

    public FastaGFFToGenomeParams withSource(String source) {
        this.source = source;
        return this;
    }

    @JsonProperty("taxon_wsname")
    public String getTaxonWsname() {
        return taxonWsname;
    }

    @JsonProperty("taxon_wsname")
    public void setTaxonWsname(String taxonWsname) {
        this.taxonWsname = taxonWsname;
    }

    public FastaGFFToGenomeParams withTaxonWsname(String taxonWsname) {
        this.taxonWsname = taxonWsname;
        return this;
    }

    @JsonProperty("taxon_reference")
    public String getTaxonReference() {
        return taxonReference;
    }

    @JsonProperty("taxon_reference")
    public void setTaxonReference(String taxonReference) {
        this.taxonReference = taxonReference;
    }

    public FastaGFFToGenomeParams withTaxonReference(String taxonReference) {
        this.taxonReference = taxonReference;
        return this;
    }

    @JsonProperty("release")
    public String getRelease() {
        return release;
    }

    @JsonProperty("release")
    public void setRelease(String release) {
        this.release = release;
    }

    public FastaGFFToGenomeParams withRelease(String release) {
        this.release = release;
        return this;
    }

    @JsonProperty("type")
    public String getType() {
        return type;
    }

    @JsonProperty("type")
    public void setType(String type) {
        this.type = type;
    }

    public FastaGFFToGenomeParams withType(String type) {
        this.type = type;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((("FastaGFFToGenomeParams"+" [fastaFile=")+ fastaFile)+", gffFile=")+ gffFile)+", genomeName=")+ genomeName)+", workspaceName=")+ workspaceName)+", source=")+ source)+", taxonWsname=")+ taxonWsname)+", taxonReference=")+ taxonReference)+", release=")+ release)+", type=")+ type)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
