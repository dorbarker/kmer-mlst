Ideas
=====

- Overall plan is to use kmers to match relevant reads
- Assemble reads or align to known alleles?
- How to identify novel alleles?

Kmer Strategies
===============


Variable sites
---------------------------

- Kmer aligns to variable locations
- Use to pre-include / exclude alleles
- Need to occasionally add or redefine kmers as new alleles discovered

Conserved sites
----------------------------

- Kmers align to non-variable sites shared by all known alleles
- Identify to locus rather than the particular allele
- Need to occasionally redefine kmers due to new variability being discovered in
  previously conserved regions

Flanks
-------------------

- Target beginning and ends of loci
- Least amount of searching
- Variation in kmers could lead to false negatives without periodic redefinition

End-to-end 
----------------

- Like exhaustive kmers
- Probably faster and maybe a few more false negatives

Exhaustive
----------------

- Get all possible kmers for gene
- Retrieve matching reads
- Assemble or align reads

Known Allele Identification
===========================

- Use kmers to identify matching reads

Assemble
--------

- Assemble with fast assembler like skesa
- Match resulting contig to alleles / match alleles to contig

Align
-----

- Align reads to alleles
- Pick allele with most matching reads or kmers


