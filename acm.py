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

GeneLengths = Dict[GeneName, Dict[AlleleName, int]]

def kmer_dict():
    return defaultdict(kmer_dict)

def revcomp(sequence: str) -> str:

    complements = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G'
        }

    return ''.join(complements[x] for x in reversed(sequence))


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


def coverage(alignment, covered_ranges, count):

    for interval in covered_ranges:
        for position in range(*interval):
            alignment[position] += count

    return alignment


def end_to_end(loci_directory: Path, kmer_size: int) -> KmerDict:

    # TODO refactor this

    kmers = kmer_dict()
    gene_expected_lengths = {}

    loci_fasta = loci_directory.glob('*.fasta')

    for locus in loci_fasta:

        locus_name = locus.stem
        gene_expected_lengths[locus_name] = {}

        with locus.open('r') as f:

            for record in SeqIO.parse(f, 'fasta'):

                allele = str(record.id)


                sequence = str(record.seq)
                seq_length = len(sequence)
                breaks = range(0, seq_length, kmer_size)

                gene_expected_lengths[locus_name][allele] = seq_length

                for start in breaks:
                    stop = start + kmer_size

                    if stop <= seq_length:

                        kmer = sequence[start : stop]
                        interval = (start, stop)

                    else:
                        # If the sequence doesn't divide evenly by the kmer
                        # length, add a new full-length kmer to the very end
                        kmer = sequence[-kmer_size: ]
                        interval = (seq_length - kmer_size, seq_length)

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

    return kmers, gene_expected_lengths


def align_kmers(
        kmer_counts: Dict[str, int],
        kmer_scheme: KmerDict,
        gene_expected_lengths: GeneLengths
        ):

    # For each kmer found, look up which gene(s) it came from, which alleles,
    # and the number of kmers aligning to those positions

    alignments = {}

    for kmer, count in kmer_counts.items():

        for locus in kmer_scheme[kmer]:

            for allele in kmer_scheme[kmer][locus]:

                try:
                    alignment = alignments[locus][allele]

                except KeyError:
                    expected_length = gene_expected_lengths[locus][allele]
                    alignment = [0 for _ in range(expected_length)]

                covered_ranges = kmer_scheme[kmer][locus][allele]

                alignment = coverage(alignment, covered_ranges, count)

                alignments[locus][allele] = alignment

    return alignments

def match_alleles(alignments):

    # Must identify missing incomplete

    potentially_novel = set()

    allele_matches = {}

    for locus in alignments:

        full_length_matches = {allele: alignment
                               for allele, alignment in
                               alignments[locus].items()
                               if all(alignment)}

        allele_matches[locus] = full_length_matches


    return allele_matches

def call_alleles(allele_matches, gene_expected_lengths):

    calls = {}

    for locus in gene_expected_lengths:

        # Locus either incomplete or not found
        if locus not in allele_matches:
            # Need logic for absent/truncated loci
            continue

        potential_matches = list(allele_matches[locus].keys())

        # The happy case - one unambiguous allele match
        if len(allele_matches[locus]) == 1:

            calls[locus] = potential_matches[0]

        # Multiple potential matches
        else:
            # Need logic to pick a winner, probably by the number of aligned
            # kmers
            continue

    return calls

@click.group()
def cli():
    pass


@cli.command()
@click.option('-k', '--kmer-length', type=int, required=True)
@click.argument('loci', type=click.Path(exists=True), nargs=1)
@click.argument('genome', type=click.Path(exists=True), nargs=-1)
def call(loci, genome, kmer_length):

    # TODO: more elegant handling of str to Path conversion
    kmer_scheme, gene_expected_lengths = end_to_end(Path(loci), kmer_length)
    print(kmer_scheme)
