#!/usr/bin/env python

from __future__ import print_function, division
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Repeat(object):
    sequence_name = ''
    start = 0
    end = 0
    copy_number = 0.0
    period_size = 0
    consensus_size = 0
    match_percent = 0
    repeat_consensus = ''
    repeat_sequence = ''
    score = 0
    indel_percent = 0
    a_count = 0
    c_count = 0
    g_count = 0
    t_count = 0
    entropy = 0.0

    def __init__(self, sequence_name='', start=0, end=0, copy_number=0.0,
                 period_size=0, consensus_size=0, match_percent=0.0,
                 repeat_consensus='', repeat_sequence='',
                 score=0, indel_percent=0.0, a_count=0, c_count=0,
                 g_count=0, t_count=0, entropy=0.0):
        self.sequence_name = sequence_name
        self.start = start
        self.end = end
        assert self.end > self.start
        self.copy_number = copy_number
        self.period_size = period_size
        self.consensus_size = consensus_size
        self.match_percent = match_percent
        self.indel_percent = indel_percent
        self.repeat_consensus = repeat_consensus
        self.repeat_sequence = repeat_sequence
        self.score = score
        self.a_count = a_count
        self.c_count = c_count
        self.g_count = g_count
        self.t_count = t_count
        self.entropy = entropy

    def __repr__(self):
        return '{} {} {} {} ({} {}-{})'.format(self.repeat_consensus, self.repeat_sequence, self.score,
                                               self.match_percent, self.sequence_name,
                                               self.start, self.end)


def read_repeats(filename):
    repeats = []
    with open(filename) as input_file:
        expect_repeat = False
        for line in input_file:
            if line.startswith('Sequence:'):
                parts = line.strip().split()
                sequence_name = parts[1]
                expect_repeat = False
            elif line.startswith('Parameters'):
                expect_repeat = True
            elif expect_repeat and line.strip() != '':
                (start_str, end_str, period_size_str, copy_number_str, consensus_size_str,
                 match_percent_str, indel_percent_str, score_str, a_count_str,
                 c_count_str, g_count_str, t_count_str, entropy_str, repeat_consensus,
                 repeat_sequence) = line.strip().split()
                # intervals reported by TRF are closed intervals in 1 based coordinates, in other
                # words if 1-5 is reported, the 1st, 2nd, 4th and 5th bases in the sequence are included
                #
                # in python string intervals are half open in 0 based coordinate, in other words
                # 1-5 means that the 2nd, 3rd, and 4th bases are included
                #
                # to convert from TRF to Python conventions, the start position is reduced by 1
                # and the end position is retained. Thus 1-5 in TRF is equivalent to 0-5 in Python
                repeat = Repeat(sequence_name=sequence_name, start=int(start_str) - 1, end=int(end_str),
                                period_size=int(period_size_str), consensus_size=str(consensus_size_str),
                                copy_number=float(copy_number_str),
                                match_percent=int(match_percent_str), indel_percent=int(indel_percent_str),
                                score=int(score_str), a_count=int(a_count_str), c_count=int(c_count_str),
                                g_count=int(g_count_str), t_count=int(t_count_str), entropy=float(entropy_str),
                                repeat_consensus=repeat_consensus, repeat_sequence=repeat_sequence)
                repeats.append(repeat)
    return repeats


def index_genome(genome_filename):
    sequence_index = SeqIO.index(genome_filename, 'fasta')
    return sequence_index


def write_repeats(sequence_index, flank_size, output_file):
    repeat_index = dict()
    for repeat in repeats:
        start = repeat.start - flank_size if repeat.start >= flank_size else 0
        end = repeat.end + flank_size  # does not matter if end > sequence length
        sequence = sequence_index[repeat.sequence_name].seq[start:end]
        index = repeat_index.get(repeat.sequence_name, 0) + 1
        repeat_index[repeat.sequence_name] = index
        seqrecord = SeqRecord(id='{}_{}'.format(repeat.sequence_name, index),
                              description='{} {}-{} {} {} {} {}'.format(repeat.sequence_name,
                                                                        repeat.start, repeat.end,
                                                                        flank_size, index,
                                                                        repeat.repeat_consensus,
                                                                        repeat.copy_number),
                              seq=sequence)
        SeqIO.write(seqrecord, output_file, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read in TRF output and extract repeats (with flanking sequences) as FASTA')
    parser.add_argument('--flank_size', default=200, type=int, help='Size of flanking region to extract before / after repeat')
    parser.add_argument('dat_filename', help='TRF output .dat filename')
    parser.add_argument('genome_filename', help='Genome filename (FASTA format)')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='FASTA format output filename')
    args = parser.parse_args()

    repeats = read_repeats(args.dat_filename)

    # here is where you could filter repeats with e.g.
    #
    # the code below only retains repeats with > 90% match percentage
    #
    # repeats_to_keep = []
    # for repeat in repeats:
    #     if repeat.match_percentage > 95:
    #         repeats_to_keep.append(repeat)
    # repeats = repeats_to_keep

    sequence_index = index_genome(args.genome_filename)
    write_repeats(sequence_index, args.flank_size, args.output_file)
