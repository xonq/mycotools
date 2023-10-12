The following commands were executed with respect to the reference MycotoolsDB
contained within this folder. Mycotools must be installed locally with the necessary 
dependencies, described in the Mycotools README (github.com/xonq/mycotools).
Full descriptions of the scripts can be reviewed using the help menu
(`-h`/`--help`) and/or reviewing the [USAGE
guide](https://github.com/xonq/mycotools/blob/master/USAGE.md).

<br /><br />

### Installing the reference MycotoolsDB locally:
```mtdb u -i <ANALYSIS_DIRECTORY> -r ust.mtdb```

<br />

### Obtaining annotation statistics:
```annotationStats $(mtdb) > annotation_stats.tsv```

<br />

### Obtaining assembly statistics:
```assemblyStats $(mtdb) > assembly_stats.tsv```

<br />

### Reconstruct phylogenies and synteny diagrams of the nitrate assimilation gene cluster using the Cluster Reconstruction and Phylogenetic Analysis Pipeline (CRAP):
```acc2locus -a ustbro1_1795 -p 1 | crap -q - -s blastp -d $(mtdb) -c <CPUS>```

<br />

### Identify 14 single-copy orthologs in each genome:
```
db2search -a hmmsearch -d $(mtdb) \
-q scos.hmm -e 2 -c <CPUS> -m 1 -o test_hmmsearch
```

<br />

### Reconstruct a phylogenomic tree referencing the retrieved single-copy orthologs:
```fa2tree -i test_hmmsearch/fastas/ -c <CPUS> -p -o test_phylo/```

<br />

### Compile the accession names of the single-copy orthologs in preparation for microsynteny tree reconstruction:
```
for i in test_hmmsearch/fastas/*fa; do grep ‘>’ $i >> scos.txt; done
sed -Ei ‘s/>//’ scos.txt
```

<br />

### Reconstruct a microsynteny tree referencing gene neighborhoods of single-copy orthologs and constrained by the topology of the phylogenomic tree:
```
db2microsyntree -c <CPUS> -d $(mtdb) -t \
test_phylo/concatenated.nex.contree -f scos.txt
```
