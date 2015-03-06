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


INSTALLATION

To install this module, run the following commands:

	perl Makefile.PL
	make
	make test
	make install

METHODS

	new
	add_ref
	add_enzyme_pair
	add_singel_enzyme
	new_enzyme
	change_range
	change_lengths_distribution_parameters
	add_output_dir
	double_digest
	single_digest
	add_SNPs
	count_SNPs_at_fragments
	add_gff
	all_frags_coverage_ratio
	frags_in_range_coverage_ratio



FORMATS OF THE INPUT REFERENCE FILE, GFF FILE, AND THE SNP COORDINATES FILE

1) RestrictionDigest is designed to process the reference file in the "fasta"
format. It should be in the following format:

>ChromsomeName1
ATCGATCGATCGATCGATCGATCGATCGATCG
>ChromsomeName2
ATCGATCGATCGATCGATCGATCGATCGATCG
The chromesome name should be unique. Any kinds of space chracters are 
not allowed, like blank space, tab etc.

In some cases, the reference file is in the following format:
>ChromsomeName1
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
>ChromsomeName2
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCG
.....
Sorry, this format is not allowed. Please transform this format into the
first format described above.

2) For the "GFF" file, it should be in the standard GFF format, here is an example:

<seqname>	<source>	<feature>	<start>	<end>	<score>	<strand>	<frame>	<attribute> 
      

##############################################################################################
scaffold39120	GLEAN	mRNA    36514	50111	.	+	.	ID=OYG_10003188_10008591;   
scaffold39120	GLEAN	Exon	36514	36534	.	+	.	Parent=OYG_10003188_10008591;   
scaffold39120	GLEAN	Exon	43276	43353	.	+	.	Parent=OYG_10003188_10008591;   
scaffold39120	GLEAN	Exon	43766	43868	.	+	.	Parent=OYG_10003188_10008591;   
scaffold39120	GLEAN	Exon	44710	44741	.	+	.	Parent=OYG_10003188_10008591;   
scaffold39120	GLEAN	Exon	49875	50111	.	+	.	Parent=OYG_10003188_10008591;   
scaffold838	GLEAN	mRNA	36580	65627	.	-	.	ID=OYG_10003281_10026064; 
###############################################################################################

The gff file should contain 9 columns and these columns must be seperated by '\t'. 

3) For the SNP coordinates file, it should contain 3 columns, here is an example:

<chromosome name>	<coordinate>	<SNP type>
scaffold18356		19		R
scaffold18356		20		Y
scaffold18356		55		G
scaffold18356		60		Y
scaffold18356		88		Y
scaffold18356		97		M
scaffold18356		114		M


EXPLANATION OF RESULT FILES

1. Files by double digest.


1) "position_frags(_in_range_)?_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND". 
The positions of digested fragments(all fragments or fragments in range) 
are saved in this file. There are 7 columns whose descriptions are as follows:
'>Name of fragments', 'strand','enzyme-front/behind', 'start position', 
'enzyme-front/behind', 'end position', 'length of fragment'.
    
    a) >Name of fragments.  The name of this fragment, the format consists 
    3 parts. The 1st part is '>'. The 2nd part is the name of scaffold in 
    the reference file User provides. The 3rd part is the serial num of this
    fragment     in all fragments produced in this scaffold.
    
    b) strand. The location of this fragment in the reference sequence, + or -.
    
    c) enzyme-front/behind. The corresponding enzyme of the overhang of front
    end of this fragments. If the strand is "+" strand, the enzyme is behind 
    enzyme, else, the enzyme is front enzyme.

    d) start position. The start position of this fragment on the specific
    scaffold. It is the start position of the     overhang, not the recognition
    base position. 

    e) enzyme-front/behind. Like the item 'c', but the opposite one to it.

    f) end postion. The end position of this fragment on the specific scaffold.
    It is also the end position of the     overhang, not the recognition base
    position.

    g) length of fragment. The base length of this fragment.


2) "seq_frags(in_range_)?_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND". 
The seqences of digested fragments(all fragments or fragments in range)
and their names are saved in this file. There are many double-line. The
first line of the double_line is the name of the fragment in the 
">chromsomeName-[0-9]+" format. The second line of the double_line is 
the sequence of this fragment.

3) "reduced_ratio_every_scaffold_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND".
The length ratio of fragments in every scaffold are saved in this file.
There are three columns in this file. The first column contains the names 
of scaffolds. And the second column contains the length ratio of all 
fragments in this scaffold. And the third column contains the length 
ratio of fragments in range in this scaffold.

4) "digestion_summary_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND". 

	Number of all fragments	  
	Number of fragments in range	
	Ratio after reducing of all fragments
	Ratio after reducing of fragments in range
	
	# The lengths' distribution of all fragments. There are three 
	columnds in this part.
	# The first column contains the refined scopes which are shaped
	by three parameters: the front split point,
	# the behind split point and the step size between the two points.
	The default parameters are: 100, 1000 and
	# 50 separately. These parameters can be changed via the method 
	<change_lengths_distribution_parameters>.
	Lengths' scope	Number of fragments in this scope Ratio of these fragments in all  
	0bp-200bp
	201bp-250b
	251bp-300bp
	301bp-350bp
	351bp-400bp
	401bp-450bp
	451bp-500bp
	501bp-550bp
	551bp-600bp
	601bp-650bp
	651bp-700bp
	701bp-750bp
	751bp-800bp
	801bp-850bp
	851bp-900bp
	901bp-950bp
	951bp-1000bp
	1001bp-bigger bp



5) "gene_coverage_ratio(_in_range)?_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND". 
If User apply the 'all_frags_coverage_ratio' or the 
'frags_in_range_coverage_ratio' method, this file will be produced. 
There are 13 columns in this file. They are: 'IntergenicLength',
'IntergenicMapLength','IntergenicMapRatio', 'GenesLength','GenesMapLength',
'GenesMapRatio', 'ExonLength', 'ExonMapLength', 'ExonMapRatio', 'IntronLength', 
'IntronMapLength','IntronMapRatio'. The first line is header. 
   
    a) GeneLength, ExonLength, IntronLength. The length of these different 
    parts, if some part does not exist, N/A is  given.
    b) GeneMapLength, ExonMapLength, IntronMapLength. The length of mapped
    part by all digested fragments or fragments in range of these parts. 
    If some part does not exists, 0 is given.
    c) GeneMapRatio, ExonMapRatio, IntronMapRatio. The ratio of mapped part
    by all digested fragments or fragments in range of these parts. If some
    part does not exits, 0 is given.

6) "SNPs_at_all_frags_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND"
This file is produced by the method 'count_SNPs_at_fragments'. It contains
2 columns: the column of chromosomes and the column of coordinates. These
SNPs are located at all restriciton fragments.

7) "SNPs_at_frags_in_range_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND"
This file is produced by the method 'count_SNPs_at_fragments'. It contains
2 columns: the column of chromosomes and the column of coordinates. These
SNPs are located at restriction fragments in range.

2. Files by single digest. 


1) "position_frags(_in_range_)?_REFERENCE_by_ENZYME". The positions of
digested fragments(all fragments or fragments in range) are saved in this
file. There are 7 columns whose descriptions are as follows:
'>Name of fragments', 'start position',  'end position', 'length of fragment'.
    
    a) >Name of fragments.  The name of this fragment, the format consists
    3 parts. The 1st part is '>'. The 2nd part is the name of scaffold in
    the reference file User provides. The 3rd part is the serial num of 
    this fragment in all fragments produced in this scaffold.
    

    b) start position. The start position of this fragment on the specific
    scaffold. It is the start position of the overhang, not the recognition 
    base position. 

   
    c) end postion. The end position of this fragment on the specific 
    scaffold. It is also the end position of the overhang, not the 
    recognition base position.

    d) length of fragment. The base length of this fragment.


2) "seq_frags(in_range_)?_REFERENCE_by_ENZYME". The seqences of digested
fragments(all fragments or fragments in range) and their names are saved 
in this file. There are many double-line. The first line of the double_line
is the name of the fragment in the ">chromsomeName-[0-9]+" format. The 
second line of the double_line is the sequence of this fragment.

3) "reduced_ratio_every_scaffold_REFERENCE_by_ENZYME". The length ratio 
of fragments in every scaffold are saved in this file. There are three 
columns in this file. The first column contains the names of scaffolds. 
And the second column contains the length ratio of all fragments in this
scaffold. And the third column contains the length ratio of fragments in 
range in this scaffold.

4) "digestion_summary_REFERENCE_by_ENZYME-FRONT_and_ENZYME-BEHIND". 

	Number of all fragments	  
	Number of fragments in range	
	Ratio after reducing of all fragments	
	Ratio after reducing of fragments in range	
	
	
	# The lengths' distribution of all fragments. There are three 
	columnds in this part.
	# The first column contains the refined scopes which are shaped 
	by three parameters: the front split point,
	# the behind split point and the step size between the two points.
	The default parameters are: 100, 1000 and
	# 50 separately. These parameters can be changed via the 
	method <change_lengths_distribution_parameters>.
	Lengths' scope	Number of fragments in this scop Ratio of these fragments in all  
	0bp-200bp
	201bp-250b
	251bp-300bp
	301bp-350bp
	351bp-400bp
	401bp-450bp
	451bp-500bp
	501bp-550bp
	551bp-600bp
	601bp-650bp
	651bp-700bp
	701bp-750bp
	751bp-800bp
	801bp-850bp
	851bp-900bp
	901bp-950bp
	951bp-1000bp
	1001bp-bigger bp

5) "genome_coverage_ratio(_in_range)?_REFERENCE_by_ENZYME". If User 
apply the 'all_frags_coverage_ratio' or the 'frags_in_range_coverage_ratio'
method, this file will be produced. 
There are 13 columns in this file. They are: 
'IntergenicLength','IntergenicMapLength','IntergenicMapRatio', 'GenesLength',
'GenesMapLength', 'GenesMapRatio', 'ExonLength', 'ExonMapLength', 'ExonMapRatio',
'IntronLength', 'IntronMapLength','IntronMapRatio'. The first line is header. 

    a)  
    b) GeneLength, ExonLength, IntronLength, UTRLength. The length of 
    these different parts, if some part does not     exist, N/A is given.
    c) GeneMapLength, ExonMapLength, IntronMapLength, UTRMapLength. The 
    length of mapped part by all digested fragments     or fragments in
    range of these parts. If some part does not exists, 0 is given.
    d) GeneMapRatio, ExonMapRatio, IntronMapRatio, UTRMapRatio. The ratio
    of mapped part by all digested fragments or fragments in range of 
    these parts. If some part does not exits, 0 is given.

6) "SNPs_at_all_frags_REFERENCE_by_ENZYME"
This file is produced by the method 'count_SNPs_at_fragments'. It contains
2 columns: the column of chromosomes and the column of coordinates. 
These SNPs are located at all restriciton fragments.

7) "SNPs_at_frags_in_range_REFERENCE_by_ENZYME"
This file is produced by the method 'count_SNPs_at_fragments'. It 
contains 2 columns: the column of chromosomes and the column of 
coordinates. These SNPs are located at restriction fragments in range.

SUPPORT AND DOCUMENTATION

After installing, you can find documentation for this module with the
perldoc command.

    perldoc RestrictionDigest

You can also look for information at:

    RT, CPAN's request tracker (report bugs here)
        http://rt.cpan.org/NoAuth/Bugs.html?Dist=RestrictionDigest

    AnnoCPAN, Annotated CPAN documentation
        http://annocpan.org/dist/RestrictionDigest

    CPAN Ratings
        http://cpanratings.perl.org/d/RestrictionDigest

    Search CPAN
        http://search.cpan.org/dist/RestrictionDigest/


LICENSE AND COPYRIGHT

Copyright (C) 2015 Jinpeng Wang, Li Li, Haigang Qi, Xuedi Du, and Guofan Zhang

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

