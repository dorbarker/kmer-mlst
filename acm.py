import click
import pandas as pd
from pathlib import Path
from ahocorasick import Automaton
from typing import List, Dict, Set
from Bio import SeqIO

KmerDict = Dict[str, Dict[str, int]]  # {sequence: {gene: allele}}

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

    for kmer in kmers:
        A.add_word(kmer)

    A.make_automaton()

    return A


