## Scripts for exploring the NCBI invertebrate collection

### find_the_flies.py

A script to explore the genomes of suborder Brachycera. This script retrieves the list of genome assemblies in the NCBI RefSeq
invertebrate collection and uses the NCBI taxonomy database to retrieve the taxonomic lineage of each organism.
It then prints out a tab separated list of scientific names, infraorder and genome assembly accession.

Usage: `find_the_flies.py my@email.org >output.txt`

where **my@email.org** is the your email address (this is used to register the query with the Entrez service as described in
the [Biopython documentation](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html)).

### extract_repeats.py

A script to process the .dat file output from TRF when it is run with the -d option and, using the genome sequence,
output a FASTA format file including the repeats and a flanking region around the repeat (by default 200 base pairs).
The script includes some commented out code that illustrates how the repeat collection can be filtered before being
written out as FASTA.

Usage: `extract_repeats.py --flank_size 100 trf_output/genome.fasta.2.5.7.80.10.50.5.dat genome.fasta repeat_regions.fasta`
