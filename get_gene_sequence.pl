#! /usr/bin/env perl
=begin


=cut

my $doc_str = <<END;

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
	
Author: Andrés García García @ Nov 2016.

END


use Getopt::Long; # To parse command-line arguments
use Bio::EnsEMBL::Registry; # From the Ensembl API, allows to conect to the db.
use Data::Dumper; # To preety print hashes easily
	local $Data::Dumper::Terse = 1;
	local $Data::Dumper::Indent = 1;
	
#===============>> BEGINNING OF MAIN ROUTINE <<=====================
	
## INITIALIZATION

	# Declare variables to hold command-line arguments
	my $gene_name = ''; my $position = 0; my $flanking = 0; my $help;
	GetOptions(
		'g|gene=s' => \$gene_name,
		'p|position=i' => \$position,
		'f|flanking=i' => \$flanking,
		'h|help' => \$help
		);
	
	# Check if user asked for help
	if($help or !$gene_name) { print_and_exit($doc_str); }
	
	## WEB DATA INITIALIZATION
	
	# Initialize a connection to the db.
	my $connection = ensembldb_connect();

    # Query for the sequence
    my $sequence = get_sequence_in($gene_name, $flanking, $position); # As a side effect it prints the gene id
    print "$sequence \n";
	
#===============>> END OF MAIN ROUTINE <<=====================


#-------------------------------------------------------------


#	===========
#	Subroutines
#	===========

sub ensembldb_connect
# Initialize a connection to the db
{
  # Initialize a registry object
  my $registry = 'Bio::EnsEMBL::Registry';

  # Connect to the Ensembl database
  print STDERR "Waiting connection to database...\n";
  $registry->load_registry_from_db(
      -host => 'ensembldb.ensembl.org', # Alternatively 'useastdb.ensembl.org'
      -user => 'anonymous'
      );
  print STDERR "Connected to database\n";
  
  return $registry;
}#------------------------------------------------------ 

sub print_and_exit
# Prints given message and exits
{
	my $message = shift;
	print $message;
	exit;
}#-----------------------------------------------------------

sub get_sequence_in
# Receieves the common name of a gene, a number of flanking bases 
# and an optional mutation position, returns the sequence of the 
# gene and prints it's name and stable id.
# The sequence is returned with the number of flanking bases specified
# and if a mutation poition is given, the returned sequence becomes 
# the nucleotide at the position plus the specified flanking bases
{
  my $gene_name = shift; # Get the passed arguments
  my $flanking = shift;
  my $position = shift;

  # Get the gene's stable id
  my $gene_id = get_geneid($gene_name);
  print "GENE: $gene_name,\tGENE_ID: $gene_id \n"; # Prints the gene id as a side effect
  # Get the gene's sequence from the id
  my $sequence = get_sequence_from_id($gene_id, $flanking, $position);

  return $sequence;
}#------------------------------------------------------

sub get_geneid
# Query the database for the stable id of the given gene
{
  my $gene_name = shift;

  # Declare a gene adaptor to get the gene
  my $gene_adaptor = $connection->get_adaptor( 'Human', 'Core', 'Gene' );
  # Declare a gene handler with the given gene
  my $gene = $gene_adaptor->fetch_by_display_label($gene_name);

  # Get the gene's EnsembleStableID
  return $gene->stable_id();
}#------------------------------------------------------

sub get_sequence_from_id
# From the stable id of a gene, query the db for the nucleotide sequence
{
	my $gene_id = shift; # Get the passed arguments
	my $flanking = shift;
	my $position = shift;

	# Declare a slice to get the sequence
	my $slice_adaptor
		= $connection 
			-> get_adaptor('Human', 'Core', 'Slice');
	  
	# Check whether a position was given or not
	my $sequence = '';
	if ($position)
	{
		# Point a slice to where the gene is located, using the gene's ID
		my $slice
			= $slice_adaptor
				-> fetch_by_gene_stable_id(
						$gene_id
						);
	
		# Decompose the slice to get important data
		my %gene = (
			'COORD_SYS' => $slice -> coord_system_name(),
			'CHROM'	=>	$slice -> seq_region_name(),
			'START'	=>	$slice -> start(),
			'END'	=>	$slice -> end()
			);
		
		print "\t@ $gene{'COORD_SYS'} $gene{'CHROM'} ($gene{'START'}-$gene{'END'})\n";
		
		# Check for any error
		my $error_string = "Position $position is not in gene @ $gene{'COORD_SYS'} $gene{'CHROM'} ($gene{'START'}-$gene{'END'})";
		die $error_string unless ($gene{'START'} <= $position and $position <= $gene{'END'});
		
		my $start = $position - $flanking;
		my $end = $position + $flanking;
		
		# Get the sequence
		$slice = $slice_adaptor->fetch_by_region( $gene{'COORD_SYS'}, $gene{'CHROM'}, $start, $end);
		$sequence = $slice -> seq();
	}
	else
	{
		# Point a slice to where the gene is located, using the gene's ID
		my $slice
			= $slice_adaptor
				-> fetch_by_gene_stable_id(
						$gene_id, 
						$flanking
						);
				
		my %gene = (
			'COORD_SYS' => $slice -> coord_system_name(),
			'CHROM'	=>	$slice -> seq_region_name(),
			'START'	=>	$slice -> start() + $flanking,
			'END'	=>	$slice -> end() - $flanking
			);
		
		print "\t@ $gene{'COORD_SYS'} $gene{'CHROM'} ($gene{'START'}-$gene{'END'})\n";
		
		$sequence = $slice -> seq();  
	}
	return $sequence;
}#------------------------------------------------------


