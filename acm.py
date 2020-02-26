import click
import gzip, bz2
import json
import pandas as pd
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


def initialize_ac_automaton(kmers: KmerDict):

    A = Automaton()

    for idx, kmer in enumerate(kmers.keys()):
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

    return kmer_counts


def coverage(alignment: List[int],
             covered_ranges: List[Position],
             count: int
             ) -> List[int]:

    for interval in covered_ranges:
        for position in range(*interval):
            alignment[position] += count

    return alignment


def end_to_end(loci_directory: Path,
               kmer_size: int
               ) -> Tuple[KmerDict, GeneLengths, AlleleKmerDict]:

    # TODO refactor this

    kmers = kmer_dict()
    alleles_kmers = kmer_dict()

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

                        alleles_kmers[locus_name][allele].add(kmer)

                    except AttributeError:
                        kmers[kmer][locus_name][allele] = [interval]
                        alleles_kmers[locus_name][allele] = set([kmer])


                    rc_kmer = revcomp(kmer)

                    try:
                        kmers[rc_kmer][locus_name][allele].append(interval)
                        alleles_kmers[locus_name][allele].add(rc_kmer)
                    except AttributeError:
                        kmers[rc_kmer][locus_name][allele] = [interval]
                        alleles_kmers[locus_name][allele] = set([rc_kmer])

    kmers = json.loads(json.dumps(kmers))
    alleles_kmers = json.loads(json.dumps(alleles_kmers))

    return kmers, gene_expected_lengths, alleles_kmers


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

def refine_matches(allele_matches, alleles_kmers, kmer_scheme, kmer_counts):

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
                             in calls[locus].items()]

            _, allele = sorted(mean_coverage)[0]

            calls[locus] = allele

    return calls


def load_model(scheme_parent):

    scheme_path = Path(scheme_parent)
    kmer_scheme_path = scheme_path.joinpath('kmer_scheme.json')
    automaton_path = scheme_path.joinpath('automaton.pickle')
    expected_length_path = scheme_path.joinpath('expected_lengths.json')

    with automaton_path.open('rb') as a:
        automaton = pickle.load(a)

    with expected_length_path.open('r') as i:
        gene_expected_lengths = json.load(i)

    with kmer_scheme_path.open('r') as k:
        kmer_scheme = json.load(k)

    return automaton, gene_expected_lengths, kmer_scheme


@click.group()
def cli():
    pass


@cli.command()
@click.argument('scheme', type=click.Path(exists=True), nargs=1)
@click.argument('genome', type=click.Path(exists=True), nargs=-1)
def call(scheme, genome):

    automaton, gene_expected_lengths, kmer_scheme = load_model(scheme)

    kmer_counts = match_kmers_to_reads(automaton, *[Path(g) for g in genome])

    alignments = align_kmers(kmer_counts, kmer_scheme, gene_expected_lengths)

    allele_matches = refine_matches(match_alleles(alignments),
                                    alleles_kmers,
                                    kmer_scheme,
                                    kemr_counts)


    calls = call_alleles(allele_matches, gene_expected_lengths)


@cli.command()
@click.option('-k', '--kmer-length', type=int, required=True)
@click.argument('loci', type=click.Path(exists=True), nargs=1)
@click.argument('scheme', type=click.Path(exists=False), nargs=1)
def build(loci, scheme, kmer_length):

    loci_path = Path(loci)
    scheme_path = Path(scheme)

    scheme_path.mkdir(exist_ok=False)
    automaton_path = scheme_path.joinpath('automaton.pickle')
    expected_length_path = scheme_path.joinpath('expected_lengths.json')
    kmer_scheme_path = scheme_path.joinpath('kmer_scheme.json')

    kmer_scheme, gene_expected_lengths, alleles_kmers = end_to_end(loci_path, kmer_length)
    automaton = initialize_ac_automaton(kmer_scheme)

    with automaton_path.open('wb') as a:
        pickle.dump(automaton, a)

    with expected_length_path.open('w') as o:
        json.dump(gene_expected_lengths, o, indent=4)

    with kmer_scheme_path.open('w') as k:
        json.dump(kmer_scheme, indent=4)
