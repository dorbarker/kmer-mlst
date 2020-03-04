import click
import gzip, bz2
import json
import pandas as pd
import numpy as np
import magic
from pathlib import Path
from ahocorasick import Automaton
from functools import partial, reduce
import itertools
from typing import List, Dict, Set, Union, Tuple, Generator
from Bio import SeqIO
from collections import defaultdict
from statistics import mean
import pickle


# Diagnostics
import sys
import cProfile

Kmer = str
GeneName = str
AlleleName = str
Position = Tuple[int, int]

KmerDict = Dict[Kmer, Dict[GeneName, Dict[AlleleName, List[Position]]]]
AlleleKmerDict = Dict[AlleleName, Set[Kmer]]

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


def initialize_ac_automaton(kmers: pd.DataFrame):

    A = Automaton()

    for idx, kmer in enumerate(set(kmers['kmer'])):
        A.add_word(kmer, (idx, kmer))

    A.make_automaton()

    return A


def yield_reads(reads: Path) -> Generator[str, None, None]:

    mime_types = {
            'text/plain': open,
            'application/x-gzip':   partial(gzip.open, mode='rt'),
            'application/gzip':     partial(gzip.open, mode='rt'),
            'application/x-bzip2':  partial(bz2.open, mode='rt'),
            'application/bzip2':    partial(bz2.open, mode='rt')
    }

    _open = mime_types[magic.from_file(str(reads), mime=True)]

    with _open(reads) as f:

        for record in SeqIO.parse(f, 'fastq'):

            sequence = str(record.seq)

            yield sequence

def match_kmers_to_reads(A: Automaton, *reads_paths) -> Dict[str, int]:


    kmer_counts = {}

    for reads in reads_paths:

        for sequence in yield_reads(reads):

            for _, (_, kmer) in A.iter(sequence):

                try:
                    kmer_counts[kmer] += 1
                except KeyError:
                    kmer_counts[kmer] = 1

    return pd.DataFrame(pd.Series(kmer_counts, name='count', dtype=int))


def coverage(locus_allele_df: pd.DataFrame, expected_length: int):

    alignment = np.zeros(expected_length)

    for _, row in locus_allele_df.dropna().iterrows():

        alignment[row['start'] : row['stop']] += row['count']

    return alignment


def end_to_end(loci_directory: Path,
               kmer_size: int
               ) -> Tuple[pd.DataFrame, pd.Series]:

    # TODO refactor this

    gene_expected_lengths = {}

    scheme_columns = ('kmer', 'allele', 'locus', 'start', 'stop')

    scheme_dataframe = []

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

                    rc_kmer = revcomp(kmer)
                    fwd_rec = [kmer, allele, locus_name, start, stop]
                    rev_rec = [rc_kmer, allele, locus_name, start, stop]

                    scheme_dataframe.append(fwd_rec)
                    scheme_dataframe.append(rev_rec)

    kmers = pd.DataFrame(scheme_dataframe, columns=scheme_columns)

    return kmers, gene_expected_lengths


def align_kmers(
        kmer_counts: pd.DataFrame,
        kmer_scheme: pd.DataFrame,
        gene_expected_lengths: GeneLengths
        ):

    # For each kmer found, look up which gene(s) it came from, which alleles,
    # and the number of kmers aligning to those positions

    alignments = {}
    start_stop = ['start', 'stop']

    kmers = kmer_scheme.join(kmer_counts, on='kmer')


    # Aggregate by locus and allele and join counts
    locus_allele_groups = (
            kmer_scheme
                .join(kmer_counts, on='kmer')
                .groupby(['locus', 'allele'])
            )


    for (locus, allele), kmer_df in locus_allele_groups:

        alignment = coverage(kmer_df, gene_expected_lengths[locus][allele])

        try:

            alignments[locus][allele] = alignment

        except KeyError:

            alignments[locus] = {allele: alignment}

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

def refine_matches(allele_matches,
                   kmer_scheme: pd.DataFrame,
                   kmer_counts: pd.DataFrame):

    # Take the allele_matches counts

    # TODO: Refactor mutating processing of allele_matches

    for locus, alleles in allele_matches.items():

        kmers_in_alleles = itertools.chain.from_iterable(alleles_kmers[locus].values())
        union_kmers_in_alleles = set(kmers_in_alleles)

        common_kmers = set(reduce(set.intersection,
                                  alleles_kmers[locus].values(),
                                  union_kmers_in_alleles))

        # For each allele, loop over alleles_kmers and subtract any kmers that hit
        # all targets
        for allele in alleles:
            for kmer in common_kmers:

                alignment = allele_matches[locus][allele]
                covered_ranges = kmer_scheme[kmer][locus][allele]

                # because we're subtracting
                negative_count = -1 * kmer_counts[kmer]

                alignment = coverage(alignment, covered_ranges, negative_count)


                allele_matches[locus][allele] = alignment


    return allele_matches


def call_alleles(allele_matches, gene_expected_lengths):

    calls = {}

    for locus in gene_expected_lengths:

        # Locus either incomplete or not found
        if locus not in allele_matches:
            # Need logic for absent/truncated loci
            calls[locus] = '0'

        potential_matches = list(allele_matches[locus].keys())

        # The happy case - one unambiguous allele match
        if len(allele_matches[locus]) == 1:

            calls[locus] = potential_matches[0]

        # Multiple potential matches
        else:

            mean_coverage = [(mean(reads), allele)
                             for allele, reads
                             in allele_matches[locus].items()]

            _, allele = sorted(mean_coverage)[-1]

            calls[locus] = allele
    return calls


def load_model(scheme_path: Path):

    kmer_scheme_path = scheme_path.joinpath('kmer_scheme.pickle')
    automaton_path = scheme_path.joinpath('automaton.pickle')
    expected_length_path = scheme_path.joinpath('expected_lengths.pickle')

    with automaton_path.open('rb') as a:
        automaton = pickle.load(a)

    with expected_length_path.open('rb') as i:
        gene_expected_lengths = pickle.load(i)

    with kmer_scheme_path.open('rb') as k:
        kmer_scheme = pickle.load(k)

    return automaton, gene_expected_lengths, kmer_scheme


def call_(scheme: Path, genomes: List[Path]):

    automaton, gene_expected_lengths, kmer_scheme = load_model(scheme)

    kmer_counts = match_kmers_to_reads(automaton, *genomes)

    alignments = align_kmers(kmer_counts, kmer_scheme, gene_expected_lengths)

    # allele_matches = refine_matches(match_alleles(alignments),
    #                                 alleles_kmers,
    #                                 kmer_scheme,
    #                                 kmer_counts)


    allele_matches = match_alleles(alignments)
    calls = call_alleles(allele_matches, gene_expected_lengths)

    diag_name = genomes[0].parent.name
    pd.DataFrame(calls, index=[diag_name]).to_csv(sys.stdout) # temporary hack solution


@click.group()
def cli():
    pass


@cli.command()
@click.argument('scheme', type=click.Path(exists=True), nargs=1)
@click.argument('genome', type=click.Path(exists=True), nargs=-1)
def call(scheme, genome):

    # click shenannigans
    scheme_path = Path(scheme)
    genome_path = [Path(g) for g in genome]

    call_(scheme_path, genome_path)


@cli.command()
@click.option('-k', '--kmer-length', type=int, required=True)
@click.argument('loci', type=click.Path(exists=True), nargs=1)
@click.argument('scheme', type=click.Path(exists=False), nargs=1)
def build(loci, scheme, kmer_length):

    loci_path = Path(loci)
    scheme_path = Path(scheme)

    scheme_path.mkdir(exist_ok=False)
    automaton_path = scheme_path.joinpath('automaton.pickle')
    expected_length_path = scheme_path.joinpath('expected_lengths.pickle')
    kmer_scheme_path = scheme_path.joinpath('kmer_scheme.pickle')

    kmer_scheme, gene_expected_lengths = end_to_end(loci_path, kmer_length)
    automaton = initialize_ac_automaton(kmer_scheme)

    with automaton_path.open('wb') as a:
        pickle.dump(automaton, a)

    with expected_length_path.open('wb') as o:
        pickle.dump(gene_expected_lengths, o)

    with kmer_scheme_path.open('wb') as k:
        pickle.dump(kmer_scheme, k)


