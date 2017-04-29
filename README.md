# Napoop
Nucleic acid sequnce utility coded in Perl

## Napoop.pm
A Perl module to handle DNA and RNA sequences

    $ perldoc Napoop.pm

    use Napoop;

## napoop
The command written in Perl and using Napoop.pm  

## DESCRIPTION
Hopefully, this script is a useful command-line tool to handle varous nucleic acid sequences, i.e. DNAs and RNAs. Originally, Napoop.pm was written to study object-oriented programming in Perl. After coding Napoop.pm, preexisting miscellanious scripts were gathered to complete this script, napoop. This has been designed to be usefull especially when it is multiply pipe-lined. The name NA, P, and OOP stand for Nucleic Acids, Perl, and Object-Oriented Programming, respectively. Currently, FASTA, GenBank, and EMBL formats are supported as input data as well as a raw sequence in a simple test.  

## INSTALL
Move napoop into one of the PATH directories, such as /usr/local/bin, and chmod +rx. Check which napoop.  
Move Napoop.pm into one of the PERL5LIB directories, such as /usr/local/lib, and chmod +r. Unless echo $PERL5LIB prints the directory, export PERL5LIB=$PERL5LIB:/usr/local/lib or something like this.  

## USAGE

    $ napoop               Is equivalent to adding "-h"
    $ napoop -a Impact     Changes the FASTA header, be longer than 1
    $ napoop -c            Prints the complementary sequence
    $ napoop -d aa         Prints amino acids
    $ napoop -d DNA        Prints the DNA codon table
    $ napoop -d RNA        Prints the RNA codon table
    $ napoop -e            Prints the FASTA header and continue
    $ napoop -f 1          Reformats the sequence in a single line
    $ napoop -f 66         Reformats the sequence which spaces
    $ napoop -f 80         Reformats the sequence
    $ napoop -g chrY       Select one FASTA sequence form multiple ones
    $ napoop -h            Shows a simple usage
    $ napoop -l            Prints the length
    $ napoop -m            Shows CpG information
    $ napoop -n            Prints all the nucleotide symbols
    $ napoop -o            Prints each nucleotide content
    $ napoop -p ref.fa     Compares two sequences
    $ napoop -r            Reverse sequence (3' to 5')
    $ napoop -s 8-16       Snip a specified region, commas are allowed
    $ napoop -t [0-6]      Translate the sequence
    $ napoop -u            Print the sequence in uppercase
    $ napoop -w            Print the sequence in lowercase
    $ napoop -z            Print each length in bp for a multiple FASTA

## OPTIONS
    -a  The FASTA heaer is replaced and printed,
        so the argument is required after "-a".
        The script continues the process,
        in other words, "-a" can be used with others,
        except for "-e".
    -c  The complementary sequence is printed.
    -d  The condon table is printed.
        An argument is required.
        If any string starting from 'D' or 'd' is followed,
        the DNA codon table is printed.
        If any string starting from 'R' or 'r', is followed,
        the RNA codon table is printed.
        In other cases, information on amino acids is printed.
    -e  The FASTA header is printed prior to the printing
        of the sequence.
        This "-e" and "-a" are mutually exclusive.
    -f  Reformatted sequence is printed.
        If the argument is 1,
        it is printed in a single line without any spaces.
        If the argument is a multiple of 11,
        it is printed with a space in every 10 bp.
        Otherwise, the sequence is printed in multiple lines
        with "argument" bp each line.
    -g  Select a single FASTA format form a multiple FASTA sequences.
        The argument should match the header.
    -h  This PD is printed.
    -l  The length of the sequence is printed.
    -m  The CpG score and G+C content are printed.
        G+C content is shown in percentage and parenthesized.
    -n  All the nucleotide symbols are printed.
    -o  Each nucleotide content is printed.
        The first, second, and third columns show
        nucleotide, count, and ratio.
        -p  Compares two sequences which are input from the standard
        input and specified this option.
        If they are identical, nothing is printed.
        Otherwise, the first mismatch position is printed.
    -r  The reverse sequence is printed.
        It is not the complementary sequence,
        but just a reverse one from 3' to 5'.
        It is useful to draw a double-stranded DNA.
    -t  Translated amino acid sequence is printed.
        0-2 and 3-5 specify the three frames in the original
        and complementary strand, respectively.
    -s  Snipped sequence is printed.
        The region should by specified like "101-200"
        which are 1-based integer coordinates.
        If the order of the two numbers are reversed,
        its complementary sequence is printed.
    -u  The sequence is printed all in uppercase.
    -w  The sequence is printed all in lowercase.
    -z  Each length in bp is printed onto the standard output
        for an input multiple FASTA file

## AUTHOR
Coded by Kohji OKAMURA, Ph.D.  
