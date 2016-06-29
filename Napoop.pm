package Napoop;

=head1 NAME

  Napoop.pm - Objects for DNA and RNA sequences

=head1 USAGE

  $ perldoc Napoop

  use Napoop;
  $gene = Napoop->new(original => $seq);
  $gene->renew(original => $seq2);
  $length_in_bp = $gene->len;
  $complementary = $gene->complementary;
  ($CpG_score, $GC_content) = $gene->cpg_score;
  %counter     = $gene->content;
  $subseq1     = $gene->snip(101, 200);
  $subseq2     = $gene->snip(200, 101);	# complementary
  $subseq3     = $gene->snip_fill_n(-1000, 2000);
  $aa0         = $gene->translate(0);
  $aa1         = $gene->translate(1);
  $aa2         = $gene->translate(2);
  $aa3         = $gene->translate(3);
  $aa4         = $gene->translate(4);
  $aa5         = $gene->translate(5);
  @dist        = $gene->periodicity('CpG', 50);
  $freq        = $gene->period_ratio('CpG', 10);
  $sequence    = $gene->get_sequence;
  $format      = $gene->get_format;
  $definition  = $gene->get_definition;
  $gene->set_definition;

=head1 DESCRIPTION

  This module is designed to treat DNA and RNA sequences.
  Currently, FASTA, GenBank, and EMBL formats are supported
  as well as a raw sequence.
  The original attribute get one of these formatted
  sequence as one scalar variable, e.g. $seq,
  as shown in SYNOPSIS.

=head1 AUTHOR

  Coded by Kohji OKAMURA, Ph.D.

=head1 HISTORY

  Mar 30, 2005  The first version of Sequence.pm
  Jan 28, 2008  Transplant from Sequence.pm
  May 02, 2008  CpG socre = 0 when n(C)*n(G) = 0
  May 02, 2008  Last modified
  Aug 04, 2008  Add warn about amino acid sequence
  Oct 27, 2008  POD was	 modified a little
  Jan 15, 2009  Subroutine snip is modified
  Jan 21, 2009  A bug in snip is fixed
  Jan 26, 2009  A method snip_fill_n was newly added
  Jan 26, 2009  cpg_score was modified
  Mar 28, 2009  sub content was added
  Apr 03, 2009  Support complementary for snip method
  May 27, 2009  Fail safe for snip_fill_n
  Mar 13, 2012  A minor change in POD
  Nov 13, 2012  Fail-safe of option -t, forcing frameshift 0
  Jun 20, 2016  Renamed as Napoop.pm from NucleicAcids.pm
  Jun 29, 2016  Released as Napoop 1.0 via GitHub

=cut


use strict;
use warnings;
use Carp;

our $AUTOLOAD;


## class data and methods
{
  my %_attributes = (  _original      => '',
                       _format        => '',
                       _sequence      => '',	# case sensitive
                       _definition    => ''     );

  ## return a list of all attributes
  sub _all_attributes
  { return keys %_attributes; }

  ## global variable to keep count of existing objects
  my $_count = 0;

  ## manage the count of existing objects
  sub get_count
  { return $_count; }

  sub _incr_count
  { ++$_count; }

  sub _decr_count
  { --$_count; }
}


## constructor
##   called from class, e.g. $obj = Napoop->new(original => $genbank);
##   This method calls an internal subroutine, _initialize.
sub new
{
  my ($class, %arg) = @_;
  my $self = bless {}, $class;
  my $original;

  if (exists $arg{original})
  { $original = $arg{original}; }
  else
  { $original = ''; }

  $class->_incr_count;
  $self->_initialize($original);
}


## renew
##   called from an object, e.g. $gene->renew(original=>$seq)
##   This subroutine resets an existing object
sub renew
{
  my ($self, %arg) = @_;
  my $original;

  if (exists $arg{original})
  { $original = $arg{original}; }
  else
  { $original = ''; }

  $self->_initialize($original);
}


## _initialize
##   called from either new or renew inside
sub _initialize
{
  my ($self, $original) = @_;
  my $sequence;

  $self->{_original} = $original;

  if ($self->{_original} =~ m/^LOCUS +(\w+?) .+?\nORIGIN\s+1 (.+?)\/\/\n/s or
      $self->{_original} =~ m/.*?\nLOCUS +(\w+?) .+?\nORIGIN\s+1 (.+?)\/\/\n/s)
  {
    $self->{_format} = 'genbank';
    $self->{_definition} = $1;
    $sequence = $2;
    $sequence =~ tr/A-Za-z//cd;
    $self->{_sequence} = $sequence;
  }
  elsif ($self->{_original} =~ m/^ID +(\w+?) .+?\nSQ\s+[^\n]+\n(.+?)\/\/\n/s or
     $self->{_original} =~ m/.*?\nID +(\w+?) .+?\nSQ\s+[^\n]+\n(.+?)\/\/\n/s)
  {
    $self->{_format} = 'embl';
    $self->{_definition} = $1;
    $sequence = $2;
    $sequence =~ tr/A-Za-z//cd;
    $self->{_sequence} = $sequence;
  }
  elsif ($self->{_original} =~ m/^[a-z\s\-\:\,0-9]+$/i
                               or $self->{_original} eq '')
  {
    $self->{_format} = 'raw';
    $self->{_definition} = 'sequence';
    $sequence = $self->{_original};
    $sequence =~ tr/A-Za-z//cd;
    $self->{_sequence} = $sequence;
  }
  elsif ($self->{_original} =~ m/^>([^\n]+)\n(.+)$/s)
  {
    $self->{_format} = 'fasta';
    $self->{_definition} = $1;
    $sequence = $2;
    $sequence =~ tr/A-Za-z//cd;
    $self->{_sequence} = $sequence;
  }
  else
  { croak "Unknown format: $original"; }

  unless ((lc $self->{_sequence}) =~ m/^[abcdghikmnrstuvwy]*$/i)
  {
    if ((lc $self->{_sequence}) =~ m/^[abcdefghiklmnpqrstuvwxyz]*$/i)
    { warn "Warning: it may be amino acid sequence.\n" }
    else
    { warn "Warning: the sequence contains unknown nucleotides.\n" }
  }

  return $self;
}


## destructor
sub DESTROY
{
  my ($self) = @_;
  $self->_decr_count;
}


## make the complementary sequence
sub complementary
{
  my ($self) = @_;
  my $sequence;

  $sequence = reverse $self->{_sequence};
  $sequence =~ tr/acgtmrwsykvhdbun/tgcakywsrmbdhvan/;
  $sequence =~ tr/ACGTMRWSYKVHDBUN/TGCAKYWSRMBDHVAN/;
  return $sequence;
}


## length
sub len
{
  my ($self) = @_;
  return length ($self->{_sequence});
}


## calculates CpG score and G+C content and return them as an array
##
## Note: If there's no C or G, the returned score is 0.
sub cpg_score
{
  my ($self) = @_;
  my ($n_c, $n_g, $n_s, $n_cpg);
  my $lseq = lc $self->{_sequence};

  $n_c = $lseq =~ tr/c//;
  $n_g = $lseq =~ tr/g//;
  $n_s = $lseq =~ tr/s//;
  $n_cpg = $lseq =~ s/cg/cg/g;
  if ($self->len == 0)
  { return (0, 0) }
  elsif ($n_c == 0 or $n_g == 0)
  { return (0, ($n_c + $n_g + $n_s) / $self->len) }
  else
  { return ($n_cpg / $n_c / $n_g * ($self->len),
           ($n_c + $n_g + $n_s) / $self->len)    }
}


## snip a specified region
sub snip
{
  my ($self, $begin, $end) = @_;
  my $len = $self->len;
  my $flag_complementary = 0;
  my $sequence;

  if ($begin > $end)
  {
    my $room = $begin; $begin = $end; $end = $room;
    $flag_complementary = 1;
  }

  if ($begin < 1)       { $begin = 1 }
  elsif ($begin > $len) { $begin = $len }
  if ($end < 0)         { $end = 1 }
  elsif ($end > $len)   { $end = $len }

  $sequence = substr $self->{_sequence}, $begin - 1, $end - $begin + 1;

  if ($flag_complementary == 1)
  {
    $sequence = reverse $sequence;
    $sequence =~ tr/acgtmrwsykvhdbun/tgcakywsrmbdhvan/;
    $sequence =~ tr/ACGTMRWSYKVHDBUN/TGCAKYWSRMBDHVAN/;
  }

  $sequence
}


## snip a specified region
## If the specified region is beyond the sequence,
## appropriate number of 'n' is filled at either or both ends.
sub snip_fill_n
{
  my ($self, $begin, $end) = @_;
  my $len = $self->len;
  my ($head, $tail) = (0, 0);	# numbers of n to be added at the ends
  my $snipped;

  if ($begin > $end) { croak "Position error in snip_fill_n: $begin-$end" }

  if ($begin < 1)       { $head = 1 - $begin; $begin = 1 }
  if ($len < $end)      { $tail = $end - $len; $end = $len }

  $snipped = substr $self->{_sequence}, $begin - 1, $end - $begin + 1;
  unless (defined $snipped)
  { croak "Error in substr: $begin, $end, $self->{_definition}" }

  ('n' x $head) . $snipped . ('n' x $tail);
}


## print each nucleotide content
sub content
{
  my ($self) = @_;
  my $len = $self->len;
  my $lower = lc $self->{_sequence};
  my @nucleotides = split //, $lower;
  my %counter;

  for (my $i = 0; $i < $len; $i++)
  {
    if (exists $counter{$nucleotides[$i]}) { $counter{$nucleotides[$i]}++ }
    else                                   { $counter{$nucleotides[$i]} = 1 }
  }

  %counter;
}


my %genetic_code =
(
  'ttt' => 'F',  'tct' => 'S',  'tat' => 'Y',  'tgt' => 'C',
  'ttc' => 'F',  'tcc' => 'S',  'tac' => 'Y',  'tgc' => 'C',
  'tta' => 'L',  'tca' => 'S',  'taa' => '*',  'tga' => '*',
  'ttg' => 'L',  'tcg' => 'S',  'tag' => '*',  'tgg' => 'W',

  'ctt' => 'L',  'cct' => 'P',  'cat' => 'H',  'cgt' => 'R',
  'ctc' => 'L',  'ccc' => 'P',  'cac' => 'H',  'cgc' => 'R',
  'cta' => 'L',  'cca' => 'P',  'caa' => 'Q',  'cga' => 'R',
  'ctg' => 'L',  'ccg' => 'P',  'cag' => 'Q',  'cgg' => 'R',

  'att' => 'I',  'act' => 'T',  'aat' => 'N',  'agt' => 'S',
  'atc' => 'I',  'acc' => 'T',  'aac' => 'N',  'agc' => 'S',
  'ata' => 'I',  'aca' => 'T',  'aaa' => 'K',  'aga' => 'R',
  'atg' => 'M',  'acg' => 'T',  'aag' => 'K',  'agg' => 'R',

  'gtt' => 'V',  'gct' => 'A',  'gat' => 'D',  'ggt' => 'G',
  'gtc' => 'V',  'gcc' => 'A',  'gac' => 'D',  'ggc' => 'G',
  'gta' => 'V',  'gca' => 'A',  'gaa' => 'E',  'gga' => 'G',
  'gtg' => 'V',  'gcg' => 'A',  'gag' => 'E',  'ggg' => 'G',

  'tty' => 'F',  'ttr' => 'L',  'tay' => 'Y',  'tgy' => 'C',
  'cay' => 'H',  'car' => 'Q',  'aay' => 'N',  'aar' => 'K',
  'agy' => 'S',  'agr' => 'R',  'gay' => 'D',  'gar' => 'E',

  'ath' => 'I',  'aty' => 'I',  'atw' => 'I',  'atm' => 'I',

  'rat' => 'B',  'rac' => 'B',  'saa' => 'Z',  'sag' => 'Z',

  'tc'  => 'S',  'ct'  => 'L',  'cc'  => 'P',  'cg'  => 'R',
  'ac'  => 'T',  'gt'  => 'V',  'gc'  => 'A',  'gg'  => 'G'
);

## translator
sub translate
{
  my ($self, $frameshift) = @_;
  my $seq, my $codon, my $aa = '';

  ## $frameshift
  ## 0, 1, 2: translates the original strand
  ## 3, 4, 5: translates the complementary strand
  ## Check the value before calling this method.

  if ($frameshift > -1 and $frameshift < 3)
  { $seq = lc substr $self->{_sequence}, $frameshift }
  elsif ($frameshift > 2 and $frameshift < 6)
  { $seq = lc substr $self->complementary, ($frameshift - 3) }
  else
  { $seq = lc $self->{_sequence} }
  $seq =~ tr/u/t/;	# RNA to DNA
  for (my $i = 0; $i < (length $seq) - 2; $i += 3)
  {
    $codon = substr $seq, $i, 3;
    if (exists $genetic_code{$codon}) { $aa .= $genetic_code{$codon} }
    else
    {
      $codon = substr $codon, 0, 2;
      if (exists $genetic_code{$codon}) { $aa .= $genetic_code{$codon} }
      else                              { $aa .= 'X' }
    }
  }
  $aa
}


## find periodicity
sub periodicity
{
  my ($self, $oligo, $max_len) = @_;
  my @distribution = ();
  my $subseq;

  unless (defined $max_len) { $max_len = 50 }
  unless (defined $oligo)   { $oligo = 'CpG' }

  $oligo = lc $oligo;
  $oligo =~ tr/abcdghikmnrstuvwy//cd;	# remove 'p'

  while ($self->{_sequence} =~ m/$oligo/g)
  {
    $subseq = substr $self->{_sequence}, $-[0], ($max_len + length $oligo);
    while ($subseq =~ m/$oligo/g)
    { if ($-[0]) { push @distribution, $-[0] } }
  }
  @distribution
}


## know the frequency of a certain periodicity
## e.g. $gene->period_ratio('CpG', 8)
##      This returns a ratio of the number of CpGs that are involved
##      in 8-bp periodicity to the number of all CpGs.
sub period_ratio
{
  my ($self, $oligo, $period) = @_;
  my ($cnt_period, $cnt_total) = (0, 0);
  my $seq = lc $self->{_sequence};

  $oligo = lc $oligo;
  $oligo =~ tr/abcdghikmnrstuvwy//cd;	# remove 'p'

  while ($seq =~ m/$oligo/g)
  {
    my $pos;
    $cnt_total++;
    if (($pos = $-[0] - $period) > -1)
    {
      if ($oligo eq substr $seq, $pos, (length $oligo))
      { $cnt_period++; next }
    }
    if (($pos = $-[0] + $period) < $self->len)
    {
      if ($oligo eq substr $seq, $pos, (length $oligo))
      { $cnt_period++ }
    }
  }

  if ($cnt_total) { return $cnt_period / $cnt_total }
  0
}


## accessor and mutator
sub AUTOLOAD
{
  my ($self, $newvalue) = @_;

  my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)$/);
  unless ($operation and $attribute)
  { croak "Method name $AUTOLOAD is not in the recognized from"; }

  if ($operation eq 'get')
  { return $self->{$attribute}; }

  ## this could be dangerous due to no reset for other attributes
  if ($operation eq 'set')
  { return ($self->{$attribute} = $newvalue); }
}


1;
