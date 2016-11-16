## Usage
```
Usage: ./get_gene_sequence.pl --gene=<gene_name> [--flanking=<flankingbases>] [--mutation_position=<position>] [--help]

Example: get_gene_sequence.pl -g MT-ATP6 -p 8993 -f 20

=======================
Gene secuence obtainer.
=======================

Script to get the nucleotide sequence of a gene, or get the sequence around the position of a mutation in the gene.

Command-line arguments:

	-g, --gene
		Name of the gene to query.
		
	-f, --flanking (optional)
		Number of bases to return in the flanking regions of the query sequence.
		
	-p, --position (optional)
		A position in the gene (in chromosomal choordinates).
		If present, returns the base at that position, or dies with
		error if the mutation is outside the boundaries of the gene.
	
	-h, --help (optional)
		Show this text and exit.
	
```
Author: Andrés García García @ Nov 2016.

 The script is written in Perl and use BioPerl and the Ensembl Perl API, detailed instructions to install those are in the [**REQUIREMENTS_INSTALL_README.md**](https://github.com/Ad115/Genetics-UNAM/REQUIREMENTS_INSTALL_README.md) file.