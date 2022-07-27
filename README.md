<p align="center">
    <img
        src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/logo.png"
    >
</p>

# PURPOSE
Bring the power of broadscale genomic analyses to the masses. 

Mycotools is a compilation of computational biology tools and database [MycotoolsDB](https://gitlab.com/xonq/mcotools/-/blob/master/mycotools/MTDB.md) software that facilitate large-scale prokaryotic and fungal comparative genomics. MycotoolsDB locally assimilates all NCBI and MycoCosm (Joint Genome Institute) into a database schema with uniform file curation, scalability, and automation as guiding principles. 

- Installation is as simple as `updateDB.py --init <DIR>`. 
- `updateDB.py --update` brings the database to the current date. 
- The `.mtdb` database format enables swift transitions from analyses with datasets of 100,000s genomes to as few as a lineage of interest. No more scrambling to acquire files for an analysis.

Mycotools is currently available as a beta subset of the full suite, excluding the database assimilation tools. MycotoolsDB is currently restricted to Ohio Supercomputer Center - if you are interested in early access, please email `konkelzach@protonmail.com`.

By integrating with the curated MycotoolsDB, Mycotools aids routine-complex tasks like retrieving `gff` or `fasta` accessions; running and compiling `fasta`s of MycotoolsDB BLAST/hmmsearches; biology-based database manipulation tools; automated phylogenetic analysis pipelines from blast to Pfam extraction to tree prediction, etc etc. Mycotools includes sets of utilities that also enable easy acquisition of batches of sequence data using `ncbiDwnld.py` and `jgiDwnld.py`.

<br />

# USAGE
Check out [README.md](https://gitlab.com/xonq/mycotools/-/tree/master/mycotools) for install and the [USAGE.md](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md) for a guide. 

<br />

# CITING
If Mycotools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the Mycotools version in line.
