import click
import json
import pandas as pd
from pathlib import Path
from ahocorasick import Automaton
from typing import List, Dict, Set, Union, Tuple
from Bio import SeqIO
from collections import defaultdict

Kmer = str
GeneName = str
AlleleName = str
Position = Tuple[int, int]

KmerDict = Dict[Kmer, Dict[GeneName, Dict[AlleleName, List[Position]]]]

def kmer_dict():
    return defaultdict(kmer_dict)

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

                for _, kmer in A.iter(sequence):

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


def end_to_end(loci_directory: Path, kmer_size: int) -> KmerDict:

    # TODO refactor this

    kmers = kmer_dict()

    loci_fasta = loci_directory.glob('*.fasta')

    for locus in loci_fasta:

        locus_name = locus.stem

        with locus.open('r') as f:

            for record in SeqIO.parse(f, 'fasta'):

                allele = str(record.id)


                sequence = str(record.seq)
                seq_length = len(sequence)
                breaks = range(0, seq_length, kmer_size)

                for start in breaks:
                    stop = start + kmer_size

                    if stop <= seq_length:

                        kmer = sequence[start : stop]
                        interval = (start, stop)

                    else:
                        # If the sequence doesn't divide evenly by the kmer
                        # length, add a new full-length kmer to the very end
                        kmer = sequence[-stop : ]
                        interval = (seq_length - stop, seq_length)

                    try:
                        kmers[kmer][locus_name][allele].append(interval)
                    except AttributeError:
                        kmers[kmer][locus_name][allele] = [interval]

                    rc_kmer = revcomp(kmer)
                    try:
                        kmers[rc_kmer][locus_name][allele].append(interval)
                    except AttributeError:
                        kmers[rc_kmer][locus_name][allele] = [interval]

    kmers = json.loads(json.dumps(kmers))
