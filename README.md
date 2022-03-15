![mycotools](https://gitlab.com/xonq/mycotools/-/raw/master/mycotools.png)

# PURPOSE
Bring the power of broadscale genomic analyses to the masses. Mycotools is a compilation of computational biology tools and database (MycotoolsDB) software designed to increase throughput analyzing fungal genomic data (JGI, NCBI, and user-inputted). MycotoolsDB is a database schema with uniform file curation, scalability, and automation as guiding principles. Installation is as simple as `updateDB.py --init <DIR>`. Updates assimilate the most recent JGI & NCBI locally via `updateDB.py -u`. Scalability enables seamless hassle-free transitions from analyses with datasets of 1000s of fungi to as few as 1 or 2 simply by inputting a full database `.db` file or an extracted `.db` obtained via `extractDB.py`. 

Mycotools is currently available as a subset of the full suite, excluding the database assimilation tools. MycotoolsDB is not available outside of the Ohio Supercomputer Center until we publish a manuscript using the tools (expected summer 2022). If you are interested in early access, please email `konkelzach@protonmail.com`.

By integrating with the curated MycotoolsDB, Mycotools aids routine-complex tasks like retrieving `gff` or `fasta` accessions; running and compiling `fasta`s of MycotoolsDB BLAST/hmmsearches; biology-based database manipulation tools; automated phylogenetic analysis pipelines from blast to Pfam extraction to tree prediction, etc etc. Mycotools includes sets of utilities that also enable easy acquisition of batches of sequence data using `ncbiDwnld.py` and `jgiDwnld.py`.

<br />

# USAGE
Check out [README.md](https://gitlab.com/xonq/mycotools/-/tree/master/mycotools) for install and the [USAGE.md](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md) for a guide. 

<br />

# CITING
If Mycotools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the Mycotools version in line.
