import click
import pandas as pd
from pathlib import Path
from ahocorasick import Automaton
from typing import List, Dict, Set, Union, Tuple
from Bio import SeqIO

KmerDict = Dict[str, Dict[str, Dict[str, Union[Tuple[int, int], int]]
#  {sequence: {gene_name: {"allele": allele_num,
#                          "position": (start, stop)
#                         }
#   }

def revcomp(sequence: str) -> str:

    complements = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G'
        }

    return ''.join(complements[x] for x in reversed(seq))


def initialize_ac_automaton(kmers: KmerDict):

    A = Automaton()

    for kmer in kmers.keys():
        A.add_word(kmer)

    A.make_automaton()

    return A

def match_kmers_to_reads(A: Automaton, *fastqs):


    kmer_counts = {}

    for fastq in fastqs:

        with fastq.open('r') as f:

            for record in SeqIO.parse(f, 'fastq'):

                sequence = str(record.seq)

                for _, kmer in A.iter(sequence)

                    try:
                        kmer_counts[kmer] += 1
                    except KeyError:
                        kmer_counts[kmer] = 1

    return kmer_counts


def coverage(covered_ranges, expected_length):

    coverage_counts = [0 for _ in range(expected_length)]

    for interval in covered_ranges:
        for position in range(*interval):
            coverage_counts[position] += 1

    return coverage_counts

