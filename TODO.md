- [ ] Conserved log class
    - [ ] Must be capable of determining if the run is congruent or new
- [ ] Annotate code (all)
- [ ] Conform old scripts to PEP8 (all)
- [ ] Build all-in-one stable conda package
- [ ] Transition code-base to Rust
- [ ] implement kon_log class, paying attention to verbosity control
- [ ] kon_log class output to local configuration directory
- [ ] support compressed databases

### crap
- [ ] *Outgroup manager for clusters that fit within min and maximum sequences*
- [ ] Percent positives filter
- [ ] Integrate agglomerative clustering
- [ ] Allow for inputing a specific run order
- [ ] Log-based resume
- [ ] Do not reiterate running a gene in the same homology group
- [ ] Allow converting HG runs' names
- [ ] *Better root inference*
- [ ] Assembly query method, i.e. through tblastn
- [ ] Allow changing the clustering variable

### curGFF3
- [ ] pseudogenes can have RNAs, and CDSs from NCBI may reference those
  pseudogene parents or their RNAs (GCA_900074715.1_LAME0)
- [ ] some pseudogenes fail because they are given an "Alias" without being
  completed (GCA_004920355.1)

### db2microsyntree
- [ ] Allow log removal

### db2search
- [ ] Distinguish between nt and aa mmseqs dbs
- [ ] Allow for blastdb construction
- [ ] Streamline mmseqs parsing
- [ ] mmseqs save db option
- [ ] profile mmseqs search
- [ ] concatenate mmseqs query dbs
- [ ] optional fail upon any failures
- [ ] Log hmmer runs
- [ ] *nhmmer option*
- [ ] create all outputs as temp files and move when complete
- [ ] extract covered portion of hits
- [ ] max hits post blast compilation
- [ ] hsp option

### dbtools
- [ ] Vectorize MTDB class
- [ ] make mtdb compiled class
- [ ] remove Entrez email login, simplify API access
- [ ] Get taxonomy of non-genus names

### extract_mtdb
- [ ] allow for a lineage list from command line (may already be integrated)
- [ ] stdin argument input
- [ ] *Fix when lineages have multiple ranks, e.g. Tremellales sp. will be
  extracted from Tremellales input, when the order is likely what's requested*

### fa2clus
- [ ] sort log by default, and only unique run parameters
- [ ] percent positive mode
- [ ] integrate MCL
- [ ] rerun aggclus on new data

### fa2hmmer2fa
- [ ] Move from extracthmm to simplified output parsing

### fa2tree
- [ ] Implement fa2clus
- [ ] ignore non-fasta inputs

### gff2svg
- [ ] find a prettier way to create SVGs
- [ ] parse for in gene coordinates and annotations
- [ ] create a single file output option for multiple inputs

### jgiDwnld
- [ ] remove gff v gff3 option

### manage_mtdb
- [ ] delete database feature
- [ ] *fix local password encryption*
- [ ] overwrite old password
- [ ] move database feature
- [ ] archive and unarchive genomes

### mtdb
- [ ] add a log option of connected MTDBs
- [ ] remove standalone scripts from PATH
- [ ] look for old ome versions in query
- [x] add a version querying option
- [ ] add an option to query log of analyses

### ncbiDwnld
- [ ] db check to ensure log is relevant to input

### predb2mtdb
- [ ] source to reference the annotation source
- [ ] *integrate prokka/bakta*
- [ ] error check FAA

### update_mtdb
- [x] optimize dereplication, currently too slow
- [ ] initial JGI predb2mtdb fails because assemblyPath doesn't exist as a
  column, but restarts are fine
- [ ] update introduction output
- [ ] need a verbose option
- [ ] reversion option
- [ ] *reference a manually curated duplicate check*
- [ ] prohibit specific IDs implementation
- [ ] finish --save
- [ ] singular strain download option
- [ ] *pull failed JGI downloads from NCBI*
- [ ] remove overlap when rerunning failed genomes
- [ ] central MTDB repository and reference option
- [ ] Improve MD5 check efficiency (update_mtdb)
- [ ] print organism name with genome accession
- [ ] don't remove files until after predb2mtdb (requires update_mtdb specific
  function)
- [ ] Need a manually curated file to correct errors in naming, e.g. Vararia v
  Vavraia, Fibularhizoctonia v Fibulorhizoctonia
- [ ] initialize from a predb
