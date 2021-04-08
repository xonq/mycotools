# PURPOSE
MycoTools is a compilation of computational biology tools and database (MycoDB) software designed to increase throughput analyzing fungal genomic data (JGI & NCBI). MycoDB is a database schema with uniform gene coordinate (`gff`) and fasta header curation, systematic code naming, and modularity as guiding principles, which enable analyses with datasets of 1000s of fungi to as few as 1 or 2. MycoDB is not available outside of the Ohio Supercomputer Center until we publish the tools. Please email `konkelzach@protonmail.com` if you are interested in becoming an early adopter.

By integrating with the curated MycoDB, MycoTools cuts time on routine tasks like retrieving `gff` or `fasta` accessions, running and compiling `fasta`s of MycoDB BLAST/hmmsearches, and other analyses like creating a phylogeny with `fa2tree.py` or wrangling your dataset with agglomerative clustering via `aggClus.py`. MycoTools includes sets of utilities that also enable easy acquisition of batches of sequence data using `ncbiDwnld.py` and `jgiDwnld.py`.

<br />

# USAGE
Check out the [USAGE.md](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md) for a guide on the possibilities MycoTools can enable in your research. 

<br />

# CITING
If MycoTools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the MycoTools version in line.
