#!/usr/bin/env python

from __future__ import print_function, division
import argparse
import ftplib
import sys
import time
import xml.etree.ElementTree as ET
from Bio import Entrez  # requires biopython - conda install biopython


class Genome(object):
    accession = ''
    taxon_id = ''
    name = ''
    lineage = ''
    
    def __init__(self, accession='', taxon_id='', name=''):
        self.accession = accession
        self.taxon_id = taxon_id
        self.name = name
        
    def __repr__(self):
        return '{} ({}: {})'.format(self.name, self.taxon_id, self.accession)


def get_assembly_summary(filename):
    conn = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
    conn.login()
    conn.cwd('/genomes/refseq/invertebrate')
    with open(filename, 'wb') as output_file:
        conn.retrbinary('RETR ' + filename, output_file.write)


def load_genome_data(filename):
    genomes = []
    with open(filename) as input_file:
        for line in input_file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            genome = Genome(accession=parts[0], taxon_id=parts[6], name=parts[7])
            genomes.append(genome)
    return genomes


def add_genome_lineages(genomes):
    for genome in genomes:
        handle = Entrez.efetch(db="taxonomy", id=genome.taxon_id)
        time.sleep(0.4)
        record = handle.read()
        parsed_record = ET.fromstring(record)
        taxon = parsed_record.find('Taxon')
        lineage = taxon.find('Lineage').text.split('; ')
        genome.lineage = lineage


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find Brachycera genomes from NCBI')
    parser.add_argument('email', help='Email address to use for Entrez')
    args = parser.parse_args()
    Entrez.email = args.email
    
    filename = 'assembly_summary.txt'
    print('Retrieving RefSeq invertebrate assembly summary', file=sys.stderr)
    get_assembly_summary(filename)
    genomes = load_genome_data(filename)
    print('Looking up lineages', file=sys.stderr)
    add_genome_lineages(genomes)
    for genome in genomes:
        if 'Brachycera' in genome.lineage:
            print(genome.name, genome.lineage[20], genome.accession, sep='\t')