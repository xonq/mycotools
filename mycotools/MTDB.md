# Mycotools Database (MTDB)

MTDBs are locally assimilated, uniformly curated databases of JGI MycoCosm (preferred) and NCBI
fungal genomic data. MTDBs are represented in tab delimitted database `.mtdb` reference files,
which serve as the scaleable input to Mycotools scripts. 

Herein are the objectives, standards, and expectations of MTDB and associated files.

<br /><br />

## OBJECTIVE

Enable broadscale comparative genomics via a systematically curated, automatically
assembled/updated, scaleable genomics database. MTDB primarily seeks to resolve several outstanding
issues in comparative genomics: 

1. Uniformly curate genomic data within and across multiple databases, i.e. the notorious
inconsistency of the gene coordinates, `gff` file
2. Promote ease-of-use for scalable and large scale analyses, i.e. transitioning between datasets
in a phylogenetic analysis
3. Keep-up with the accelerating deposition of public genomic data via automatic updates
4. Implement a modular comparative genomic analyses toolkit to enable automated pipelining and 
make routine comparative genomic analyses accessible

<br /><br />

## MYCOTOOLSDB standard
### `.mtdb` file format standard
Tab-delimited file, with no column headers, organized into columns:
'''
'ome', 'genus', 'species', 'strain', 'taxonomy',
'version', 'source', 'biosample', 'assembly_acc',
'published', 'acquisition_date', 'fna', 'faa', 'gff3'
'''

- `ome`: MTDB accession "ome" - first three letters of genus, first three letters of species
(or "sp.") upon generation with a unique number, e.g. `cryneo24` `[a-zA-Z0-9.]`
- `genus`: Genus name from NCBI/Mycocosm/inputted master table; `[a-zA-Z]`
- `species`: Species name from NCBI/Mycocosm/inputted master table; `[a-zA-Z]`
- `strain`: Strain name with only letters and numbers from NCBI/Mycocosm master table `[a-zA-z0-9]`
- `taxonomy`: NCBI taxonomy `JSON` object derived from genus
- `version`: Mycocosm version/NCBI modification date
- `source`: Genome source, e.g. 'ncbi'/'jgi'/'lab'; `[a-z0-9]`
- `biosample`: NCBI biosample accession
- `assembly_acc`: NCBI GenBank/RefSeq assembly accession or Mycocosm genome_code
- `published`: Publication metadata or binary publication response
- `acquisition_date`: Data of input into master database `YYYYmmdd`
- `fna`: optional assembly `.fna`, required for manual inputs, default `$MYCOFNA/fna/<ome>.fa`
- `faa`: optional proteome `.faa`, required for manual inputs, default `$MYCOFAA/faa/<ome>.aa.fa`
- `gff3`: optional gene coordinate `.gff3`, required for manual inputs, default `$MYCOGFF3/gff3/<ome>.gff3`

<br />

### accession formatting
All MTDB aliases will be formatted as `<ome>_<acc>` where `ome` is `ome` in the MTDB.
For JGI, accessions will pull from the `protein_id` field in the gene coordinates file.
For NCBI, accessions will pull from the `product_id` field in the gene coordinates file.
For entries without a detected protein ID, an alias will be assigned with the
prefix, 'mtdb'. Pseudogenes, tRNAs, and rRNAs aliases will format as
`<ome>_<type><type_count>`

<br />

### `gff3` gene coordinates file standard
MTDB attempts to curate, assimilate, and modernize ALL JGI and NCBI gene coordinates files, including
legacy data. All proteins, associated transcripts, exons, gene, and CDSs will include `;Alias=<ome>_<acc>`
at the end of the attributes column.

Permitted sequence type fields: `{"exon", "gene", "mRNA", "rRNA", "tRNA", "CDS", "pseudogene"}`; introns
will be curated to exons. 

<br />

### JGI and NCBI data assimilation
`updateDB.py` will prioritize Mycocosm (JGI) genomes over NCBI if and *only if* a
single BioSample accession can be retrieved by cross-referencing the Mycocosm
entry with the GOLD master table. Because multiple BioSamples may be associated
with a single species, and because there is no direct tie from the Mycocosm table to
GOLD, the only way to reliably retrieve BioSamples for Mycocosm entries is if
there is only one BioSample associated with a specific Mycocosm strain. Because
Mycocosm is continually pushing their results to NCBI, failure to associate a
Mycocosm accession with a BioSample code will likely be resolved by downloading
from NCBI if the data is not usage restricted.

<br />

### Local data assimilation
Locally annotated genomes can be added to the database by filling out and
submitting a `.predb` file using `predb2db.py`. `predb2db.py` will curate the
inputted data and output into the current directory. Once complete,
`updateDB.py -a <PREDB_RESULT>` will add the `.mtdb` from `predb2db.py` to the
master database.

<br /><br />

## Database management

The master MTDB should be generated from one user, and privileges should be
distributed using `chmod`. Note, the master user/group are the only ones with 
privileges to merge manually curated `predb` files into the master database.
