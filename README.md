# RestrictionDigest
A powerful Perl module for simulating the genome digestions

DESCRIPTION

RestrictionDigest can simulate the reference genome digestion and
generate comprehensive information of the simulation. It can simulate single-enzyme
digestion, double-enzyme digestion and size selection process. It can also analyze
multiple genomes at one run and generates concise comparison of enzyme(s) performance
across the genomes. 

For more information, please see the academic paper about this module which is posted 
on  ResearchGate (https://www.researchgate.net/publication/297717326_RestrictionDigest_A_powerful_Perl_module_for_simulating_genomic_restriction_digests) or published online (http://www.sciencedirect.com/science/article/pii/S071734581630001X).


### 2022-12-26
**Important updates made to the module!**

I used the Cartesian product function to get all possible combinations for enzyme cut sites with degenerate bases. This update makes any enzyme can be added to the enzyme pool. The use of the Cartesian production function also simplified the programming sentences.

Another benefit of using the Cartesian production function is the compatability of the enzymes often used in 2b-RAD method, like BsaXI. As BsaXI cuts the DNA sequence at two different sites, users can use two pseudo restriction enzymes to simulate the BsaXI enzyme. Below is an example:

1. list the cutting sites of BsaXI:
    N|NNNNNNNNNACNNNNNCTCCNNNNNNNNNN|N
2. the first pseudo enzyme cut at the front cutting-site, for example we name it as BsaXIF, and its cutting site is:
    N|NNNNNNNNNACNNNNNCTCCN
3. the second pseudo enzyme cut at the behind cutting-site, for example we name it as BsaXIB, and its cutting site is:
    NACNNNNNCTCCNNNNNNNNNN|N
4. we add the two pseudo enzymes to the module and use them to perform double-enzyme digest.

The Perl Script can be written like this:
```
  #!/usr/bin/perl -w 
  use strict; 
  use IO::File; 
  use lib "./"; 
  use RestrictionDigest; 
  my $dg = RestrictionDigest::SingleItem::Double->new();  
  $dg->add_ref(-reference=>"example.reference.fa"); 
  $dg->new_enzyme(-enzyme_name=>'BsaXIF', -recognition_site=>'N|NNNNNNNNNACNNNNNCTCCN'); 
  $dg->new_enzyme(-enzyme_name=>'BsaXIB', -recognition_site=>'NACNNNNNCTCCNNNNNNNNNN|N'); 
  $dg->add_enzyme_pair(-front_enzyme=>"BsaXIF",-behind_enzyme=>'BsaXIB'); 
  $dg->add_output_dir(-output_dir=>"output"); 
  $dg->double_digest(); 
```
Users should keep in mind that 2b-RAD enzymes can cut DNA sequences into small pieces, the program will need large memory space. So it is best to run it in a computer with large memory space or in a cluster with enough memory. If the memory space is not available, consider to cut chromosomes or scaffolds into smaller ones.
