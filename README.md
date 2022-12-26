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


# 2022-12-26
Important updates made to the module!
I used the Cartesian product function to get all possible combinations for enzyme cut sites with degenerate bases. This update makes any enzyme can be added to the enzyme pool. The use of the Cartesian production function also simplified the programming sentences.
