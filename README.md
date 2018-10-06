## Scripts for exploring the NCBI invertebrate collection

### find_the_flies.py

A script to explore the genomes of suborder Brachycera. This script retrieves the list of genome assemblies in the NCBI RefSeq
invertebrate collection and uses the NCBI taxonomy database to retrieve the taxonomic lineage of each organism.
It then prints out a tab separated list of scientific names, infraorder and genome assembly accession.

Usage: `find_the_flies.py my@email.org >output.txt`

where **my@email.org** is the your email address (this is used to register the query with the Entrez service as described in
the [Biopython documentation](https://biopython.org/DIST/docs/api/Bio.Entrez-module.html)).