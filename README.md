# RestrictionDigest
A powerful Perl module for simulating the genome digestions

DESCRIPTION

RestrictionDigest is used for the simulation in silico of single-enzyme
and double-enzyme GBS/RAD approach. GBS(Genotyping by sequencing) and 
RAD (Restrict site associated  DNA markers) are two popular methods for
reducing the complexity of the whole genome. As sequencing the whole 
genome is expensive and sometimes unnecessary, sequencing the representative
part of the whole genome becomes practicable and attractive. The most 
important thing of representative sequencing is choosing the right enzyme(s)
to digest the whole genome. The appropriate enzyme(s) should satisfy several
criteria including: the recognition sites of the enzyme(s) are evenly 
dispersed on each single chromsome;  the recognition sites of the  enzyme(s)
are not or rarely located at the repetitive region of the genome; the 
fragments produced by the enzyme(s) digestion are neither too long nor 
too short, they should be in the suitable range which is most effective 
for PCR amplification. 

In order to evaluate the effect of particular enzyme(s) on the reference, 
we developed this perl module. RestrictionDigest can simulate the digestions
of the whole reference genomes and provide comprehensive profiles of these 
digestions. RestrictionDigest generates essential information of restriciton 
fragments; the essential information includes the number of restriction 
fragments and the distribution of restriction fragments lengths. The essential
information help researchers to determine whether the candidate enzyme(s)
can be used to digest the genome DNA or not. Besides the essential information, 
RestrictionDigest also provide supplenmentary information of restriciton 
fragments by analysing the GFF file and the SNP coordinate file. 



