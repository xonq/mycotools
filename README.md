# PURPOSE
Mycotools is a compilation of computational biology tools and database (MycotoolsDB) software designed to increase throughput analyzing fungal genomic data (JGI & NCBI). MycotoolsDB is a database schema with uniform file curation, scalability, and automation as guiding principles. Installation as simple as `updateDB.py --init <DIR>`. Updates by assimilating the most recent JGI & NCBI data to manipulate locally via `updateDB.py -u`. Scalability enables analyses with datasets of 1000s of fungi to as few as 1 or 2 simply by inputting a full database `.db` file or an extracted `.db` portion. MycotoolsDB is not available outside of the Ohio Supercomputer Center until we publish the tools. Please email `konkelzach@protonmail.com` if you are interested.

By integrating with the curated MycotoolsDB, Mycotools cuts time on routine tasks like retrieving `gff` or `fasta` accessions, running and compiling `fasta`s of MycotoolsDB BLAST/hmmsearches, automated phylogenetic pipelines, and database manipulation tools. MycoTools includes sets of utilities that also enable easy acquisition of batches of sequence data using `ncbiDwnld.py` and `jgiDwnld.py`.

<br />

# USAGE
Check out [README.md](https://gitlab.com/xonq/mycotools/-/tree/master/mycotools) for install and the [USAGE.md](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md) for a guide. 

<br />

# CITING
If Mycotools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the Mycotools version in line.
