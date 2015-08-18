
package RestrictionDigest;

use 5.8.0;
use strict;
use warnings FATAL => 'all';



=head1 NAME

RestrictionDigest - A simulation tool for reducing the genome with one DNA endonuclease or a pair DNA endonucleases!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 DESCRIPTION


Next Generation Sequencing (NGS) is developing quickly and many cost-saving approaches are developed based on NGS.  As sequencing 
the whole genome is expensive and sometimes unnecessary, sequencing the representative part of the whole genome becomes practicable 
and attractive. The most famous two approaches are Genotyping By Sequencing (GBS) and Restrict site Associated DNA sequencing (RAD).
These two approaches can reduce the  ratio of sequenced part of the whole genome which will save a lot of money. Both of them 
utilize tpye II DNA endonuclease to digest the whole genome. In order to select the correct enzyme or enzyme-pair, we developed 
this perl module: RestrictionDigest.

RestrictionDigest is used for the simulation of single/double-enzyme GBS/RAD approach. The most important thing of representative sequencing 
is choosing the right enzymes to digest the whole genome. The appropriate enzymes should satisfy several criteria including: the 
recognition sites of these two enzymes are evenly dispersed on each single chromsome;  the recognition sites of these two enzymes 
are not or rarely located at the repetitive region of the genome; the fragments produced by two-enzyme digestion are neither too long 
nor too short, they should be in the suitable range which is most  effective for PCR amplification; the complexity reducing rate 
shoule be in the suitable range(1%-10%).  RestrictionDigest can produce several result files including the fragments seqence and their 
locations on the reference. The ratio of these fragments on the reference seqences are also produced. If User can provide the gff 
annotation file, RestrictionDigest can also provide the coverage ratio of these fragments on gene region, even CDS, INTRON, UTR region, 
and intergenic region.

=cut





package RestrictionDigest::Double;

use 5.8.0;
use strict;
use warnings FATAL => 'all';
=head1 NAME

RestrictionDigest::Double - for double-enzyme simulation.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';



=head1 SYNOPSIS

RestrictionDigest::Double is used for the simulation of double-enzyme GBS/RAD approach. 


    use RestrictionDigest;

    $digest = RestrictionDigest::Double->new();
    
    $digest->add_ref(-reference=>'Full path of the reference file');

    $digest->new_enzyme(-enzyme_name=>'Ncol',-recognition_site=>'C|CATGG');

    $digest->add_enzyme_pair(-front_enzyme=>'EcoRI',-behind_enzyme=>'HinfI');

    $digest->add_output_dir(-output_dir=>'Full path of the output directory');

    $digest->double_digest();




=cut

use IO::File;
use File::Basename;

use vars qw(%fields %enzyme_list);

# Define the fields to be used as identifiers.
%fields = (
  ref			=> undef,
  gff =>  undef,
  SNPs => undef,
  enzyme1 => undef,
  enzyme2 =>  undef,
  output_dir=>  undef,
  range_start => 201,
  range_end => 500,
  lengths_distribution_start =>100,
  lengths_distribution_end  => 1000,
  lengths_distribution_step => 50,
  sequence_type => '100PE',
  sequence_end => 'front_enzyme',
);

# Make the endonuclease enzyme-container--a hash.
%enzyme_list=(
      AleI  =>  'CACNN|NNGTG',
      AluI  =>  'AG|CT',
      AsiSI =>  'GCGAT|CGC',
      EcoRI => 'G|AATTC',
      HinfI => 'G|ANTC',
      AluI => 'AG|CT',
      AvaII => 'G|GWCC',
      BamHI => 'G|GATCC',
      HindIII => 'A|AGCTT',
      MseI => 'T|TAA',
      MspI => 'C|CGG',
      NdeI => 'CA|TATG',
      PstI => 'CTGCA|G',
      SalI => 'G|TCGAC',
      XbaI => 'T|CTAGA',
      XhoI => 'C|TCGAG',
      NlaIII => 'NCATG|N',
      AatII =>  'GACGT|C',
      AccI  =>  'GT|MKAC',
      Acc65I  =>  'G|GTACC',
      AciI  =>  'C|CGC',
      AclI  => 'AA|CGTT',
      AfeI  => 'AGC|GCT',
      AflII   => 'C|TTAAGA',
      AflIII  => 'A|CRYGT',
      AgeI  => 'A|CCGGT',
      AlwNI  => 'CAGNNN|CTG',
      ApaI  => 'GGGCC|C',
      ApaLI   => 'G|TGCAC',
      ApeKI =>  'G|CWGC',
      ApoI  => 'R|AATTY',
      AscI  => 'GG|CGCGCC',
      AseI => 'AT|TAAT',
      AvaI => 'C|YCGRG',
      BclI =>  'T|GATCA',
      AvrII => 'C|CTAGG',
      BaeGI => 'GKGCM|C',
      BanI => 'G|GYRCC',
      BanII => 'GRGCY|C',
      BbvCI =>  'CC|TCAGC',
      BsaAI => 'YAC|GTR',
      BsaJI => 'C|CNNGG',
      BsaWI => 'W|CCGGW',
      BfaI  => 'C|TAG',
      BglII => 'A|GATCT',
      BlpI  => 'GC|TNAGC',
      BmgBI => 'CAC|GTC',
      BmtI  => 'GCTAG|C',
      Bpu10I => 'CC|TNAGC',
      BsaAI => 'YAC|GTR',
      BsaBI => 'GATNN|NNATC',



);

=head1 SUBROUTINES/METHODS

=head2 new
     
     $digest=RestrictionDigest::Double->new(); 
          Create a new object.

=cut

sub new {
  my $class = shift;
  my $self = {
    %fields,
  };
  bless ($self, $class);
  return $self;
}

=head2 add_ref
     
     $digest->add_ref(-reference=>'Full path of the reference file');
          Add the reference file  to the RestrictionDigest variable. 
          
=cut

sub add_ref {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-reference$/){
      $self->{ref}=$parameters{$_};
    }
    else{
    die "Unacceptable parameters in the method <add_ref>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }

  # Check the format of reference
  my $ref_fh=IO::File->new("$self->{ref}",'r');
  my $line_num=0;
  print"Examing the format of the reference file...\t";
  while(<$ref_fh>){
    chomp;
    $line_num++;
    if($line_num % 2 != 0){
      unless($_=~/^>\S+$/){
        die "The format of reference file is not allowed, please convert it.\nYou can find the example format in the README!\n";
      }
    }
    if($line_num % 2 ==0 ){
      unless($_=~/^[AGTCNRYMKSWHBVD]+$/i){
        die"The sequences of reference contain some charactors other than AGTCNRYMKSWHBVD, please correct that!\n";
      }
    }
    
  }
  print"the format is good.\n";
  $ref_fh->close;
}

=head2 add_enzyme_pair
     
     $digest->add_enzyme_pair(-front_enzyme=>'EcoRI', -behind_enzyme=>'HinfI');
          Add  two enzymes to the RestrictionDigest variable. 
          The names of the two enzymes are case-insensitive.

=cut



sub add_enzyme_pair {

  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-front_enzyme$/){
      $self->{enzyme1}=$parameters{$_};
    }
    elsif($_=~/^-behind_enzyme$/){
      $self->{enzyme2}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters in the method <add_enzyme_pair>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  } 

  # Check the existence of enzymes in the enzyme-container.
  my ($enzyme1_exists, $enzyme2_exists)=qw(0 0);
  foreach(keys %enzyme_list){
    if($_=~/^$self->{enzyme1}$/i){
      $enzyme1_exists++;
    }
    if($_=~/^$self->{enzyme2}$/i){
      $enzyme2_exists++;
    }
  }
  if($enzyme1_exists >0 ){
    
  }else{
    die "The front restrict endonuclease  User provides does not exist.\nHowever User can add this enzyme and its recognition sites to enzyme-container via the method 'new_enzyme' .\n";
  }
  if($enzyme2_exists >0 ){
    
  }else{
    die "The behind restrict endonuclease  User provides does not exist.\nHowever User can add this enzyme and its recognition sites to enzyme-container via the method 'new_enzyme' .\n";
  }
}

=head2 new_enzyme

     $digest->new_enzyme(-enzyme_name=>'EcoRI', -recognition_site=>'G|AATTC');
     
     If the enzyme User wants to use does not exists in the default enzyme-pool, the enzyme can be added to the enzyme-pool.
     
     For now, the enzymes in the default enzyme pool are  as follows:
           EcoRI,HinfI,AluI,AvaII,BamHI,HindIII,MseI,MspI,NdeI,PstI,SalI,XbaI,XhoI,NlaIII.
     
     Adding new enzyme(s) to the enzyme-pool is easy but the enzyme added to the pool is temporary. The added enzyme(s) is active
     only in the class context to which the 'new_method' belongs.

=cut 

sub new_enzyme {
  my $self=shift;
  my %parameters=@_;
  my $new_enzyme_name;
  my $new_enzyme_site;
  for (keys %parameters){
    if($_=~/^-enzyme_name$/){
      $new_enzyme_name=$parameters{$_};
    }
    elsif($_=~/^-recognition_site$/){
      $new_enzyme_site=$parameters{$_};
      unless($new_enzyme_site=~/|/){
        die "The 'recognition_site' must contain a cut flag '|'. Please perldoc RestrictionDigest for example or read README for moer help\n";
      }
    }
    else{
      die"Unacceptable parameters in the method <new_enzyme>, please perldoc RestrictionDigest for example or read README for more help\n";
    }
  }
  $enzyme_list{$new_enzyme_name}=$new_enzyme_site;
}


=head2 change_range

     $digest->change_range(-start=>201, -end=>500);
          Usually, not all fragments are needed by the scientists. The fragments in particular range are more attractive to them.
          RestrictionDigest provides the sequences and the positions on the reference of these fragments in range. The default range is 
           [201bp, 500]. User can change this range via the method <change_range>. 


=cut

sub change_range {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-start$/){
      $self->{range_start}=$parameters{$_};
    }
    elsif($_=~/^-end$/){
      $self->{range_end}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters in the method <change_range>,please perldoc RestrictionDigest for example or read README for more help.\n";
    }
  }
  if($self->{range_start} < $self->{range_end} ) {
    
  }else{
    die "The first parameter of the method <change_range> must be smaller than the second parameter\n";
  }
}

=head2 change_lengths_distribution_parameters
    
     $digest->change_lengths_distribution_parameters(-front=>100,-behind=>1000,-step=>50);
        This program will give a simple summary of the lengths of all fragments after digestion. By default, three scopes, which are  
         <=100bp, 101bp-1000bp and >1000bp, are given. The number of the fragments whose length fallen in these scopes are given. More 
         fine sorted, the 101bp-1000bp scope is split by the step of 50bp. These three parameters: the front edge and the behind edge 
         of the scope and the step can be changed via the method 'change_lengths_distribution_parameters';
=cut

sub change_lengths_distribution_parameters {
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-front$/){
      $self->{lengths_distribution_start}=$parameters{$_};    
    }
    elsif($_=~/^-behind$/){
      $self->{lengths_distribution_end}=$parameters{$_};
    }
    elsif($_=~/^-step$/){
      $self->{lengths_distribution_step}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters in the method <change_lengths_distribution_parameters>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
}


=head2 add_output_dir
    
     $digest->add_output_dir(-output_dir=>'Full path of the output directory');
          Adds the directory into which the result files are put by the program. 
          The directory must be in the absolute path format. 

=cut


sub add_output_dir {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-output_dir$/){
      $self->{output_dir}=$parameters{$_};
    }
    else{
      die "Unacceptable parameter in the method <add_output_dir>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
  if(-d $self->{output_dir} ) {
  }else{
    die "The output directory User provides does not exist.\n";
  }
}

=head2 double_digest

     $digest->double_digest();
          Execute the digestion process. 
          This process will produce several result files which will be located in the output directory 
          added through the add_output_dir method. 
          This process will take some time. It depends on the size of the reference file.
          During the process, Prompting messages will be output to the STDOUT 
          like "Digestion: >ScaffoldName is done!".

=cut

sub double_digest {

  my $self=shift;
  my $ref=$self->{ref};
  if($ref){
    my $name_reference=basename($ref);
    print "The reference file User provides is:\t",$name_reference,"\n";
  }else{
    die "No reference file provied! Please perldoc RestrictionDigest for example.\n";
  }
  my $front_enzyme=$self->{enzyme1};
  if($front_enzyme){
    print "The front restrict endonuclease User provides is:\t", $front_enzyme, "\n";
  }
  else{
    die"No front restrict endonuclease provided! Please perldoc RestrictionDigest for example.\n";
  }
  my $later_enzyme=$self->{enzyme2};
  if($later_enzyme){
    print "The behind restrict endonuclease User provides is:\t", $later_enzyme, "\n";
  }
  else{
    die "No behind restrict endonuclease provided! Please perldoc RestrictionDigest for example.\n";
  }
  my $range1=$self->{range_start};
  my $range2=$self->{range_end};
  print "The seperate range of the fragments User wants to output is:\t$range1 to $range2 bp.\n";
  my $output_dir=$self->{output_dir};
  if($output_dir){
    print "The output directory User provides is:\t",$output_dir, "\n";
  }
  else{
    die "No output directory provided! Please perldoc RestrictionDigest for example.\n";
  }
  $output_dir=~s/^(.+)\/?$/$1/;
  
  # Make the file handle to the reference file.
  my $ref_fh=IO::File->new("$ref",'r');
  
  # Get the basename of reference,not include the path.
  my $name_reference=basename($ref);
  
  # Make the file handles of result files.
  my $all_seq_file_fh=IO::File->new(">>$output_dir/seq_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $all_loc_file_fh=IO::File->new(">>$output_dir/position_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $seq_in_range_file_fh=IO::File->new(">>$output_dir/seq_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $loc_in_range_file_fh=IO::File->new(">>$output_dir/position_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $reduced_ratio_every_scaffold_file_fh=IO::File->new(">>$output_dir/reduced_ratio_every_scaffold_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $summary_digestion_fh=IO::File->new(">>$output_dir/digestion_summary_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  
  # Get the base-compositions of these two enzymes which User inputs.
  my ($front_seq,$later_seq);
  foreach(keys %enzyme_list){
    if($_=~/^$front_enzyme$/i){
      $front_seq=$enzyme_list{$_};
    }
    if($_=~/^$later_enzyme$/i){
      $later_seq=$enzyme_list{$_};
    }
  }

  #get the relative location of cutting of the recognation sequences; $front_enzyme_cut_loc and $later_enzyme_cut_loc are two relative locations we wanted;
  my @front_enzyme_rec_seq=split //, $front_seq;
  my @later_enzyme_rec_seq=split //, $later_seq;
  my @front_enzyme_cut_seq;
  my @later_enzyme_cut_seq;
  my $front_enzyme_rec_start=0;
  my $later_enzyme_rec_start=0;
  my $front_enzyme_cut_loc;
  my $later_enzyme_cut_loc;
  foreach(@front_enzyme_rec_seq){
    $front_enzyme_rec_start++;
    if($_=~/\|/){
      $front_enzyme_cut_loc=$front_enzyme_rec_start;
    }
    else{
      push @front_enzyme_cut_seq,$_;
    }
  }
  foreach(@later_enzyme_rec_seq){
    $later_enzyme_rec_start++;
    if($_=~/\|/){
      $later_enzyme_cut_loc=$later_enzyme_rec_start;
    }
    else{
      push @later_enzyme_cut_seq, $_;
    }
  }
  
  # Replace the degenerate bases of front enzyme;
  my $num_of_front_enzyme_cut_seq=@front_enzyme_cut_seq;
  for(@front_enzyme_cut_seq){
    if($_=~/N/){$_="[ATCG]";}
    if($_=~/R/){$_="[AG]";}
    if($_=~/Y/){$_="[CT]";}
    if($_=~/M/){$_="[AC]";}
    if($_=~/K/){$_="[GT]";}
    if($_=~/S/){$_="[GC]";}
    if($_=~/W/){$_="[AT]";}
    if($_=~/H/){$_="[ATC]";}
    if($_=~/B/){$_="[GTC]";}
    if($_=~/V/){$_="[GAC]";}
    if($_=~/D/){$_="[GAT]";}
  }
  my $seq_string_front_enzyme=join '', @front_enzyme_cut_seq;



  # Replace the degenerate bases of later enzyme;
  my $num_of_later_enzyme_cut_seq=@later_enzyme_cut_seq;
  foreach(@later_enzyme_cut_seq){
    if($_=~/N/){$_="[ATCG]";}
    if($_=~/R/){$_="[AG]";}
    if($_=~/Y/){$_="[CT]";}
    if($_=~/M/){$_="[AC]";}
    if($_=~/K/){$_="[GT]";}
    if($_=~/S/){$_="[GC]";}
    if($_=~/W/){$_="[AT]";}
    if($_=~/H/){$_="[ATC]";}
    if($_=~/B/){$_="[GTC]";}
    if($_=~/V/){$_="[GAC]";}
    if($_=~/D/){$_="[GAT]";}
  }
  my $seq_string_later_enzyme=join '', @later_enzyme_cut_seq;


  
  # Make the scalars and arrays used in the whole programme;
  my $rate_overall=0;
  my $rate_in_range=0;
  my @lengths_of_fragments=();
  my @lengths_fragments_in_range=();
  my $overall_fragments_length=0;
  my $overall_fragments_in_range_length=0;
  my $overall_length_of_scfd=0;

  ########## CIRCLES START HERE!###########

  # Process one scaffold one circle;
  my $scaffold_name; my $scaffold_seq;
  while(<$ref_fh>){
    chomp;
    if($_=~/^>/){
      $scaffold_name=$_;
    }
    else{
        print"Digesting $scaffold_name...\t";
        $scaffold_seq=$_;
        $reduced_ratio_every_scaffold_file_fh->print("$scaffold_name\t");
        my $length_of_scafd=length $scaffold_seq;
        
        # Set the arrays which contain cut locations of front enzyme and later enzyme;
        my @front_enzyme_locs=();
        
        my @later_enzyme_locs=();
        
        my $seq_for_search=$scaffold_seq;
        # Find locs of recognition of front enzyme in this scaffold;the locs are meaningful in this array; not in the human context;
          my $string=$seq_string_front_enzyme;
          if($string=~/^(\w+)$/){
            my $loc_in_array;
            $loc_in_array=index($seq_for_search,$string);
            unless($loc_in_array==-1){
              push @front_enzyme_locs, $loc_in_array;
            }
            while($loc_in_array!=-1){
              $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
              unless($loc_in_array==-1){
                push @front_enzyme_locs, $loc_in_array;
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $back_bases=$3;
            my $generate_bases=$2;
            my @generate_bases=split //, $generate_bases;
            my $seq_for_search=$scaffold_seq;
            for my $one_base(@generate_bases){
              $string=$breast_bases.$one_base.$back_bases;
              my $loc_in_array;
              $loc_in_array=index($seq_for_search,$string);
              unless($loc_in_array==-1){
                push @front_enzyme_locs, $loc_in_array;
              }
              while($loc_in_array!=-1){
                $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                unless($loc_in_array==-1){
                  push @front_enzyme_locs, $loc_in_array;
                }
              }
            }
          }

          elsif($string=~/^\[(\w+)\](\w+)\[(\w+)\]$/){
            my $g1_bases=$1;
            my $middle=$2;
            my $g2_bases=$3;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                $string=$one_base.$middle.$two_base;
                my $loc_in_array;
                  $loc_in_array=index($seq_for_search,$string);
                  unless($loc_in_array==-1){
                    push @front_enzyme_locs, $loc_in_array;
                  }
                  while($loc_in_array!=-1){
                    $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                    unless($loc_in_array==-1){
                      push @front_enzyme_locs, $loc_in_array;
                    }
                  }                
              }
            }
          }

          elsif($string=~/(\w+)\[(\w+)\](\w+)\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $middle=$3;
            my $back_bases=$5;
            my $g2_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                $string=$breast_bases.$one_base.$middle.$two_base.$back_bases;
                  my $loc_in_array;
                  $loc_in_array=index($seq_for_search,$string);
                  unless($loc_in_array==-1){
                    push @front_enzyme_locs, $loc_in_array;
                  }
                  while($loc_in_array!=-1){
                    $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                    unless($loc_in_array==-1){
                      push @front_enzyme_locs, $loc_in_array;
                    }
                  }
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $g2_bases=$3;
            my $back_bases=$5;
            my $g3_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            my @g3_bases=split //, $g3_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                for my $three_base(@g3_bases){
                  $string=$breast_bases.$one_base.$two_base.$three_base.$back_bases;
                    my $loc_in_array;
                    $loc_in_array=index($seq_for_search,$string);
                    unless($loc_in_array==-1){
                      push @front_enzyme_locs, $loc_in_array;
                    }
                    while($loc_in_array!=-1){
                      $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                      unless($loc_in_array==-1){
                        push @front_enzyme_locs, $loc_in_array;
                      }
                    }
                }
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $g2_bases=$3;
            my $back_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                  $string=$breast_bases.$one_base.$two_base.$back_bases;
                    my $loc_in_array;
                    $loc_in_array=index($seq_for_search,$string);
                    unless($loc_in_array==-1){
                      push @front_enzyme_locs, $loc_in_array;
                    }
                    while($loc_in_array!=-1){
                      $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                      unless($loc_in_array==-1){
                        push @front_enzyme_locs, $loc_in_array;
                      }
                    }
              }
            }
          }
        
  
        
        
        # Find locs of recognition of later enzyme in this scaffold;the locs are meaningful in this array; not in the human context;
   
          $string=$seq_string_later_enzyme;
          if($string=~/^(\w+)$/){
            my $loc_in_array;
            $loc_in_array=index($seq_for_search,$string);
            unless($loc_in_array==-1){
              push @later_enzyme_locs, $loc_in_array;
            }
            while($loc_in_array!=-1){
              $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
              unless($loc_in_array==-1){
                push @later_enzyme_locs, $loc_in_array;
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $back_bases=$3;
            my $generate_bases=$2;
            my @generate_bases=split //, $generate_bases;
            for my $one_base(@generate_bases){
              $string=$breast_bases.$one_base.$back_bases;
                my $loc_in_array;
                $loc_in_array=index($seq_for_search,$string);
                unless($loc_in_array==-1){
                  push @later_enzyme_locs, $loc_in_array;
                }
                while($loc_in_array!=-1){
                  $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                  unless($loc_in_array==-1){
                    push @later_enzyme_locs, $loc_in_array;
                  }
                }
            }
          }

          elsif($string=~/^\[(\w+)\](\w+)\[(\w+)\]$/){
            my $g1_bases=$1;
            my $middle=$2;
            my $g2_bases=$3;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                $string=$one_base.$middle.$two_base;
                  my $loc_in_array;
                  $loc_in_array=index($seq_for_search,$string);
                  unless($loc_in_array==-1){
                    push @later_enzyme_locs, $loc_in_array;
                  }
                  while($loc_in_array!=-1){
                    $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                    unless($loc_in_array==-1){
                      push @later_enzyme_locs, $loc_in_array;
                    }
                  }
              }
            }
          }

          elsif($string=~/(\w+)\[(\w+)\](\w+)\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $middle=$3;
            my $back_bases=$5;
            my $g2_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                $string=$breast_bases.$one_base.$middle.$two_base.$back_bases;
                  my $loc_in_array;
                  $loc_in_array=index($seq_for_search,$string);
                  unless($loc_in_array==-1){
                    push @later_enzyme_locs, $loc_in_array;
                  }
                  while($loc_in_array!=-1){
                    $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                    unless($loc_in_array==-1){
                      push @later_enzyme_locs, $loc_in_array;
                    }
                  }
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $g2_bases=$3;
            my $back_bases=$5;
            my $g3_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            my @g3_bases=split //, $g3_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                for my $three_base(@g3_bases){
                  $string=$breast_bases.$one_base.$two_base.$three_base.$back_bases;
                    my $loc_in_array;
                    $loc_in_array=index($seq_for_search,$string);
                    unless($loc_in_array==-1){
                      push @later_enzyme_locs, $loc_in_array;
                    }
                    while($loc_in_array!=-1){
                      $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                      unless($loc_in_array==-1){
                        push @later_enzyme_locs, $loc_in_array;
                      }
                    }
                }
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $g2_bases=$3;
            my $back_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                  $string=$breast_bases.$one_base.$two_base.$back_bases;
                    my $loc_in_array;
                    $loc_in_array=index($seq_for_search,$string);
                    unless($loc_in_array==-1){
                      push @later_enzyme_locs, $loc_in_array;
                    }
                    while($loc_in_array!=-1){
                      $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                      unless($loc_in_array==-1){
                        push @later_enzyme_locs, $loc_in_array;
                      }
                   }
              }
            }
          }

        
    
     
        # Integrate the front enzyme locs and later enzyme locs into a hash;
        my %hash;
        while(@front_enzyme_locs){
          my $loc=shift @front_enzyme_locs;
          $hash{$loc}=$front_enzyme;
        }

        while(@later_enzyme_locs){
          my $loc=shift @later_enzyme_locs;
          $hash{$loc}=$later_enzyme;
        }
        
        # Sort the locs by the locations of them and restore the fragments' locs with different ends.
        my @locs_sorted=sort{$a<=>$b} keys %hash;
        my @locs_wanted_before;
        my @enzyme_wanted_before;
        my @locs_wanted_after;
        my @enzyme_wanted_after;
        my ($after_loc,$enzyme_before, $enzyme_after);
        my $before_loc=shift @locs_sorted;
        while(@locs_sorted){
          $after_loc=shift @locs_sorted;
          $enzyme_before=$hash{$before_loc};
          $enzyme_after=$hash{$after_loc};
          unless($enzyme_before eq $enzyme_after){
            push @locs_wanted_before, $before_loc;
            push @enzyme_wanted_before, $enzyme_before;
            push @locs_wanted_after,$after_loc;
            push @enzyme_wanted_after, $enzyme_after;
          }
          $before_loc=$after_loc;
        }

        
        # Define the scalars which contains the length of fragments;
        my $nums_of_fragments=@locs_wanted_before;
        my @pairs_of_locs=(0..$nums_of_fragments-1);
        my $length_of_fragments=0;
        my $length_of_fragments_in_range=0;
        my $count_all=0;
        my $count_in_range=0;
        my ($before_cut_coorr,$after_cut_coorr,$bf_ct_crr_hm_rd,$af_ct_crr_hm_rd);

   
        

        while(@locs_wanted_before){
        my  $before_coorr=shift @locs_wanted_before;
        my $before_enzyme=shift @enzyme_wanted_before;
        my  $after_coorr=shift @locs_wanted_after;
        my  $after_enzyme=shift @enzyme_wanted_after;
        my $length_of_fragment;my $seq_selected;
          $count_all++;
          if($front_enzyme=~/$before_enzyme/i){
            $before_cut_coorr=$before_coorr+$front_enzyme_cut_loc-1;
            $after_cut_coorr=$after_coorr+$later_enzyme_cut_loc-2;
            $length_of_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_selected=substr($scaffold_seq, $before_cut_coorr, $length_of_fragment);     
          }
          elsif($later_enzyme=~/$before_enzyme/i){
                $before_cut_coorr=$before_coorr+$later_enzyme_cut_loc-1;
                $after_cut_coorr=$after_coorr+$front_enzyme_cut_loc-2;
                $length_of_fragment=$after_cut_coorr-$before_cut_coorr+1;
                $seq_selected=substr($scaffold_seq, $before_cut_coorr, $length_of_fragment);
          }
         my $bf_ct_crr=$before_cut_coorr+1;
         my $af_ct_crr=$after_cut_coorr+1;
              
              # Output the sequence and fragment positions of all-fragments to files;
              $all_seq_file_fh->print("$scaffold_name-$count_all\n$seq_selected\n");
              my $strand;
              if($before_enzyme=~/$front_enzyme/i){$strand="+";}
              elsif($before_enzyme=~/$later_enzyme/i){$strand="-";}
              $all_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$before_enzyme\t$bf_ct_crr\t$after_enzyme\t$af_ct_crr\t$length_of_fragment\n");
              push @lengths_of_fragments, $length_of_fragment;
                        
              # We treat the fragments in range in the same way of all fragments as above;
              if($length_of_fragment>=$range1){
                if($length_of_fragment<=$range2){
                  $count_in_range++;
                  $seq_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_selected\n");
                  $loc_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$before_enzyme\t$bf_ct_crr\t$after_enzyme\t$af_ct_crr\t$length_of_fragment\n");
                  push @lengths_fragments_in_range,$length_of_fragment;
                  $length_of_fragments_in_range+=$length_of_fragment;
                }
              }
              
              # Count all fragments length no matter in range or not!
              $length_of_fragments+=$length_of_fragment;


        }
        my $reduced_rate=$length_of_fragments/$length_of_scafd;
        my $fragment_in_range_rate_scafd=$length_of_fragments_in_range/$length_of_scafd;
        $overall_fragments_length+=$length_of_fragments;
        $overall_fragments_in_range_length+=$length_of_fragments_in_range;
        $overall_length_of_scfd+=$length_of_scafd;
        $reduced_ratio_every_scaffold_file_fh->print("$reduced_rate\t$fragment_in_range_rate_scafd\n");
        print"Done!\n"; 
    }
  }
  my $length_ratio_all_frags=$overall_fragments_length/$overall_length_of_scfd;
  my $length_ratio_frags_in_range=$overall_fragments_in_range_length/$overall_length_of_scfd;
  my $num_all_frags=@lengths_of_fragments;
  my $num_frags_in_range=@lengths_fragments_in_range;
  my $lengths_fragments_in_range=\@lengths_fragments_in_range;
  

  # Output the summary of digestion.
  $summary_digestion_fh->print("Expected all fragments:\t$num_all_frags\n");
  $summary_digestion_fh->print("Expected fragments in range:\t$num_frags_in_range\n");
  $summary_digestion_fh->print("Expected reducing ratio of all fragments:\t$length_ratio_all_frags\n");
  $summary_digestion_fh->print("Expected reducing ratio of fragments in range:\t$length_ratio_frags_in_range\n");

  my $lengths_distribution_start=$self->{lengths_distribution_start};
  my $lengths_distribution_end=$self->{lengths_distribution_end};
  my $lengths_distribution_step=$self->{lengths_distribution_step};
  my $lengths_distribution_tmp;
  my (@starts,@ends,@counts,%lengths_distribution);

  my $num_pair=0;

  for($lengths_distribution_tmp= $lengths_distribution_start;$lengths_distribution_tmp<$lengths_distribution_end+$lengths_distribution_step;$lengths_distribution_tmp=$lengths_distribution_tmp+$lengths_distribution_step){
    if($lengths_distribution_tmp== $lengths_distribution_start){
      $num_pair++;
      push @starts, '0';
      push @ends , $lengths_distribution_tmp;
      $lengths_distribution{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;

      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      $lengths_distribution{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      
    }
    elsif($lengths_distribution_tmp  < $lengths_distribution_end ){
      if($lengths_distribution_tmp+$lengths_distribution_step <$lengths_distribution_end ){
        $num_pair++;
        my $tmp_start=$lengths_distribution_tmp+1;
        my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
        push @starts, $tmp_start;
        push @ends, $tmp_end;
        $lengths_distribution{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        
    }
    elsif($lengths_distribution_tmp +$lengths_distribution_step >= $lengths_distribution_end){
      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_end;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      $lengths_distribution{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;

      $num_pair++;
      $tmp_start=$tmp_end+1;
      
      push @starts, $tmp_start;
      push @ends, 'bigger';
      $lengths_distribution{$num_pair}{"${tmp_start}bp-bigger bp"}=0;
      
      }
    }
  }
  


  for (@lengths_of_fragments){
    my $length=$_;
    my $num_pair=1;
    my $num_last_pair=@starts;
    for($num_pair=1;$num_pair<=$num_last_pair;$num_pair++){
      
      my $array_index=$num_pair-1;
      my $start=$starts[$array_index];
      my $end=$ends[$array_index];
      if($num_pair<$num_last_pair){
        if($length>=$start && $length<=$end){
          for my $description (keys %{$lengths_distribution{$num_pair}}){
            $lengths_distribution{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
      elsif($num_pair==$num_last_pair){
        if($length>=$start){
          for my $description (keys %{$lengths_distribution{$num_pair}}){
            $lengths_distribution{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
    }
  }
  
  $summary_digestion_fh->print("\nLengths\' scope\tNumber of fragments in this scope\tRatio of these fragments in all\n");
  
  for(1..@starts){
    my $num_pair=$_;
    for my $description(keys %{$lengths_distribution{$num_pair}}){
      my $num_fragments=$lengths_distribution{$num_pair}{$description};
      my $ratio;
      if($num_all_frags>0){ 
        $ratio=$num_fragments/$num_all_frags;
        }else{
          $ratio='N/A';
        }
      $summary_digestion_fh->print("$description\t$num_fragments\t$ratio\n");
    }
  }
  $summary_digestion_fh->print("\n");
}

=head2 add_SNPs
  
    $digest->add_SNPs(-SNPs=>'Full path of the SNPs file');
        Adds the absolute path of the SNPs file to RestrictionDigest variable.
=cut

sub add_SNPs {
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-SNPs$/){
      $self->{SNPs}=$parameters{$_};
    }
    else{
      die "Unacceptable parameters of the method <add_SNPs>, please perldoc RestrictionDigest for example or read README for more help\n";
    }
  }
}

=head2 count_SNPs_at_fragments
  
    $digest->count_SNPs_at_fragments();
        Count the expected SNPs appeared in the framgents.
=cut

sub count_SNPs_at_fragments {
  my $self=shift;
  my %parameters=@_;
  
  if(keys %parameters){
    for(keys %parameters){
      if($_=~/^-sequence_type$/){
        $self->{sequence_type}=$parameters{$_};
      }
      elsif($_=~/^-sequence_end$/){
        $self->{sequence_end}=$parameters{$_};
      }
      else{
        die "Unacceptable parameters in the method <count_SNPs_at_fragments>, please perldoc RestrictionDigest for example or read README for more help\n";
      }
    }
  }
  
  my $ref=$self->{ref};
  my $front_enzyme=$self->{enzyme1};
  my $later_enzyme=$self->{enzyme2};
  my $range1=$self->{range_start};
  my $SNPs=$self->{SNPs};
  if($SNPs){
    my $SNPs_basename=basename($SNPs);
    print "The SNPs file provided is\t", $SNPs_basename, ".\n"; 
  }
  else{
    die"No SNPs file provided! Please perldoc RestrictionDigest for example.\n";
  }

  my $sequence_type=$self->{sequence_type};

  # Check the sequencing type and determine the domains of the fragments to calculate the SNPs number.
  my ($pe_sequence_length, $se_sequence_length);
  my $sequence_type_code;
  if($sequence_type=~/^(\d+)pe$/i){
    $sequence_type_code=2;
    my $sequence_length=$1;
    unless($sequence_length*2 <= $range1){
      die"Unacceptable sequence_type when range's front border is $range1 bp. Please perldoc RestrictionDigest for example or read README for more help!\n";
    }
    else{
      $pe_sequence_length=$sequence_length;
    }
  }
  elsif($sequence_type=~/^(\d+)se$/i){
    $sequence_type_code=1;
    my $sequence_length=$1;
    unless($sequence_length<=$range1){
      die "Unacceptable sequence_type when range's front border is $range1 bp. Please perldoc RestrictionDigest for example or read README for more help!\n";
    }
    else{
      $se_sequence_length=$sequence_length;
    }
  }
  else{
    die "Unacceptable sequence_type. Please perldoc RestrictionDigest for example or read README for more help!\n";
  }

  print "The sequencing type provided is\t", $sequence_type, ".\n";

  # Check the sequencing end of single-end sequencing.
  my $se_sequence_end;
  my $se_sequence_end_code;
  if($sequence_type_code==1){
    $se_sequence_end=$self->{sequence_end};
    if($se_sequence_end=~/^front_enzyme$/i){
      $se_sequence_end_code=1;
      print "The sequencing end provided is\t", $se_sequence_end, ".\n";
    }
    elsif($se_sequence_end=~/^behind_enzyme$/i){
      $se_sequence_end_code=2;
      print "The sequencing end provided is\t", $se_sequence_end, ".\n";
    }
    else{
      die "Unacceptable sequence_end. Please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }


  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  # Get the basename of reference,not include the path.
  my $name_reference=basename($ref);

  # Make the file handle to the summary file;
  my $summary_digestion_fh=IO::File->new(">>$output_dir/digestion_summary_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");

  # Make the file handle to the locs file of all fragments and fragments in range.
  my $all_loc_file_fh=IO::File->new("$output_dir/position_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  my $loc_in_range_file_fh=IO::File->new("$output_dir/position_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  
  # Hold SNPs into hash.
  my $SNPs_fh=IO::File->new("$SNPs",'r');
  my %SNPs;
  while(<$SNPs_fh>){
    chomp;
    my ($SNPs_scaffold_name,$SNPs_pos,$type)=split /\s+/, $_;
    push @{$SNPs{$SNPs_scaffold_name}}, $SNPs_pos;
  }
  print "Counting SNPs at the output fragments...\t";
  
  # Create the SNPs files to hold SNPs of all fragments and fragments in range.
  my $SNPs_at_all_frags_fh=IO::File->new(">>$output_dir/SNPs_at_all_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $SNPs_at_frags_in_range_fh=IO::File->new(">>$output_dir/SNPs_at_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");



  # Scan all locs file to count SNPs and output the result to the summary file.
  my $all_SNPs_number=0;
  if($sequence_type_code == 1){
    if($se_sequence_end_code==1){
      while(<$all_loc_file_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          for my $SNP(@{$SNPs{$scaffold_name}}){
            if($front_enzyme=~/$enzyme1/i){
              if($SNP >=$pos1 && $SNP <=($pos1+$se_sequence_length-1)){
                $all_SNPs_number++;
                $SNPs_at_all_frags_fh->print("$scaffold_name\t$SNP\n");
              }
            }
            elsif($front_enzyme=~/$enzyme2/i){
              if($SNP >= ($pos2-$se_sequence_length+1) && $SNP <=$pos2){
                $all_SNPs_number++;
                $SNPs_at_all_frags_fh->print("$scaffold_name\t$SNP\n");
              }
            }
          }
        }
      }
    }
    elsif($se_sequence_end_code==2){
      while(<$all_loc_file_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          for my $SNP(@{$SNPs{$scaffold_name}}){
            if($later_enzyme=~/$enzyme1/i){
              if($SNP>=$pos1 && $SNP<=($pos1+$se_sequence_length-1)){
                $all_SNPs_number++;
                $SNPs_at_all_frags_fh->print("$scaffold_name\t$SNP\n");
              }
            }
            elsif($later_enzyme=~/$enzyme2/i){
              if($SNP>= ($pos2-$se_sequence_length+1) && $SNP<=$pos2){
                $all_SNPs_number++;
                $SNPs_at_all_frags_fh->print("$scaffold_name\t$SNP\n");
              }
            }
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected SNPs at all fragments via $sequence_type:\t$all_SNPs_number\n");
  }
  elsif($sequence_type_code ==2){
    while(<$all_loc_file_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        for my $SNP(@{$SNPs{$scaffold_name}}){
          if( ($SNP>=$pos1 && $SNP<=($pos1+$pe_sequence_length-1)) || ( $SNP>= ($pos2-$pe_sequence_length+1) && $SNP<=$pos2 ) ){
          $all_SNPs_number++;
          $SNPs_at_all_frags_fh->print("$scaffold_name\t$SNP\n");
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected SNPs at all fragments via $sequence_type:\t$all_SNPs_number\n");
  }

  
  # Scan locs in range file to count SNPs and output the result to the summary file.

  my $SNPs_in_range_number=0;

  if($sequence_type_code == 1){
    if($se_sequence_end_code==1){
      while(<$loc_in_range_file_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          for my $SNP(@{$SNPs{$scaffold_name}}){
            if($front_enzyme=~/$enzyme1/i){
              if($SNP>=$pos1 && $SNP<=($pos1+$se_sequence_length-1)){
                $SNPs_in_range_number++;
                $SNPs_at_frags_in_range_fh->print("$scaffold_name\t$SNP\n");
              }
            }
            elsif($front_enzyme=~/$enzyme2/i){
              if($SNP>= ($pos2-$se_sequence_length+1) && $SNP<=$pos2){
                $SNPs_in_range_number++;
                $SNPs_at_frags_in_range_fh->print("$scaffold_name\t$SNP\n");
              }
            }
          }
        }
      }
    }
    elsif($se_sequence_end_code==2){
      while(<$loc_in_range_file_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          for my $SNP(@{$SNPs{$scaffold_name}}){
            if($later_enzyme=~/$enzyme1/i){
              if($SNP>=$pos1 && $SNP<=($pos1+$se_sequence_length-1)){
                $SNPs_in_range_number++;
                $SNPs_at_frags_in_range_fh->print("$scaffold_name\t$SNP\n");
              }
            }
            elsif($later_enzyme=~/$enzyme2/i){
              if($SNP>= ($pos2-$se_sequence_length+1) && $SNP<=$pos2){
                $SNPs_in_range_number++;
                $SNPs_at_frags_in_range_fh->print("$scaffold_name\t$SNP\n");
              }
            }
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected SNPs at fragments in range via $sequence_type:\t$SNPs_in_range_number\n");
  }
  elsif($sequence_type_code ==2){
    while(<$loc_in_range_file_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        for my $SNP(@{$SNPs{$scaffold_name}}){
          if( ($SNP>=$pos1 && $SNP<=($pos1+$pe_sequence_length-1)) || ( $SNP>= ($pos2-$pe_sequence_length+1) && $SNP<=$pos2 ) ){
          $SNPs_in_range_number++;
          $SNPs_at_frags_in_range_fh->print("$scaffold_name\t$SNP\n");
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected SNPs at fragments in range via $sequence_type:\t$SNPs_in_range_number\n");
  }
  print "Done!\n";
}


=head2 add_gff

    $digest->add_gff(-ref=>'Full path of the gff file');
        Adds the absolute path of gff file to RestrictionDigest variable. 
        If User want to calculate the coverage ratio of different genomic strutures by the digested 
        fragments, this method is obligatory.
=cut

sub add_gff { 
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-gff$/){
      $self->{gff}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters of the method <add_gff>, please perldoc RestrictionDigest for example or read README for more help\n";
    }
  }
}

=head2 all_frags_coverage_ratio

     $digest->all_frags_coverage_ratio();
          Calculate the coverage ratio of different parts of the genome reference, 
          including mRNA, CDS, UTR and intergenic region, by all digested fragments.
          This process will take a very long time, please be patient.  
=cut

sub all_frags_coverage_ratio {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $e1=$self->{enzyme1};
  my $e2=$self->{enzyme2};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $loc_file_all_frags="$output_dir/position_frags_${basename_ref}_by_${e1}_and_${e2}";

  $self->genome_structure_coverage_ratio($loc_file_all_frags);
}

=head2 frags_in_range_coverage_ratio

     $digest->frags_in_range_coverage_ratio();
          Calculate the coverage ratio of different parts of the genome reference, 
          including mRNA, CDS, UTR and intergenic region, by the fragments in range.
          This process will take a very long time, please be patient. 
=cut

sub frags_in_range_coverage_ratio {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $e1=$self->{enzyme1};
  my $e2=$self->{enzyme2};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $loc_file_all_frags="$output_dir/position_frags_in_range_${basename_ref}_by_${e1}_and_${e2}";

  $self->genome_structure_coverage_ratio($loc_file_all_frags);
}

=head2  genome_structure_coverage_ratio

        A subroutine invoked by the 'all_frags_coverage_ratio' and 
        the 'frags_in_range_coverage_ratio' methods. 
        User do not use this subroutine.

=cut

 
sub genome_structure_coverage_ratio {
  my $self=shift;
  my $loc_frags_file=shift;

  # Get the reference name, enzyme names from hash.
  my $ref=$self->{ref};
  my $ref_fh=IO::File->new($ref,'r');

  my $basename_ref=basename($ref);
  my $e1=$self->{enzyme1};
  my $e2=$self->{enzyme2};
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  # Get the full path and name of gff file and make a filehandle to it.
  my $gff=$self->{gff};
   
  if($gff){
    my $gff_basename=basename($gff);
    print "The GFF file provided is\t", $gff_basename, ".\n";
  }
  else{
    die"No GFF file provided! Please perldoc RestrictionDigest for example.\n";
  }
  
  # Make  output filehandles to generic region and intergenic region coverage ratio files. 
  my $coverage_ratio_fh;
  
  my $basename_loc_frags_file=basename($loc_frags_file);

  if($basename_loc_frags_file =~/in_range/){
    $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_ratio_in_range_${basename_ref}_by_${e1}_and_${e2}");
  }
  else{
    $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_ratio_${basename_ref}_by_${e1}_and_${e2}");
  }

  $coverage_ratio_fh->print("ScaffoldName\tIntergenicLength\tIntergenicMapLength\tIntergenicMapRatio\tGenesLength\tGenesMapLength\tGenesMapRatio\tExonsLength\tExonsMapLength\tExonsMapRatio\tIntronsLength\tIntronsMapLength\tIntronsMapRatio\n");


  # Set variables used in this program.
  my($all_intergenic_length,$all_intergenic_map_length,$all_gene_length,$all_gene_map_length,$all_exon_length,$all_exon_map_length,$all_intron_length,$all_intron_map_length)=qw(0 0 0 0 0 0 0 0);
  


  # Fetch the lengths of different parts of this scaffold from the GFF file.
      my $gff_fh=IO::File->new($gff,'r');
      # Create two hashes to check if start or stop position exists.
      my %gff; my %exon_check; my %gene_check;
      print "Scanning the GFF file...\t";

      my($gff_gene_count,$gff_exon_count)=qw(0 0);
      while(<$gff_fh>){
        chomp;
        unless($_=~/^#/){
          my @gff=split /\t/, $_;
          if(@gff == 9){
           my ($gff_scfd_name,$gff_type,$gff_start,$gff_stop)=@gff[0,2,3,4];
              if($gff_type=~/CDS|exon/){
                $gff_type='exon';
                unless($exon_check{$gff_scfd_name}{start}{$gff_start}){
                  $gff_exon_count++;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_exon_count}{start}=$gff_start;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_exon_count}{stop}=$gff_stop;
                  $exon_check{$gff_scfd_name}{start}{$gff_start}=1;
                }
              }
              elsif($gff_type=~/mRNA|gene/){
                $gff_type='gene';
                unless($gene_check{$gff_scfd_name}{start}{$gff_start}){
                  $gff_gene_count++;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_gene_count}{start}=$gff_start;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_gene_count}{stop}=$gff_stop;
                  $gene_check{$gff_scfd_name}{start}{$gff_start}=1;
                }
              }
          }
        }
      }
      $gff_fh->close;
      print "Done!\n";

      # Scan the loc file and fetch the information.
      
      my $loc_frags_file_fh=IO::File->new($loc_frags_file,'r');
      my %all_locs;
      print "Scanning the positions file...\t";
      while(<$loc_frags_file_fh>){
        chomp;
        my($frgt_name,$strand,$front_enzyme,$start_coor,$later_enzyme,$stop_coor,$length)=split /\t/,$_;
        my $scfd_frgt;
        if($frgt_name=~/^>(.+)-\d+$/){
          $scfd_frgt=$1;
          push @{$all_locs{$scfd_frgt}}, ($start_coor, $stop_coor);
        }
      }
      $loc_frags_file_fh->close;
      print "Done!\n";


  # Handle every scaffold one by one.
  my($scfd_name);
  while(<$ref_fh>){
    chomp;
    if($_=~/^>(\S+)$/){
      $scfd_name=$1;
    }
    else{
      my $scaffold_length=length $_;
      print "CoverageRatio: Processing $scfd_name...\t";


      # Calculate the lengths of different parts of this scaffold.
      my ($exon_length,$intron_length,$gene_length,$intergenic_length)=qw(0 0 0 0);
      
      if($gff{$scfd_name}{gene}){
        for my $gene_order(keys %{$gff{$scfd_name}{gene}}){
          my $front=$gff{$scfd_name}{gene}{$gene_order}{start};
          my $behind=$gff{$scfd_name}{gene}{$gene_order}{stop};
          $gene_length+=$behind-$front+1;
        }
        $intergenic_length=$scaffold_length-$gene_length;
      }else{
        $intergenic_length=$scaffold_length;
      }

      if($gff{$scfd_name}{exon}){
        for my $exon_order(keys %{$gff{$scfd_name}{exon}}){
          my $front=$gff{$scfd_name}{exon}{$exon_order}{start};
          my $behind=$gff{$scfd_name}{exon}{$exon_order}{stop};
          $exon_length+=$behind-$front+1;
        }
        $intron_length=$gene_length-$exon_length;
      }

      $all_intergenic_length+=$intergenic_length;
      $all_gene_length+=$gene_length;
      $all_exon_length+=$exon_length;
      $all_intron_length+=$intron_length;     

      # Scan the locs file and fetch the map lengths of different parts.
      my $gene_map_length=0;
      my $exon_map_length=0;
      my $intron_map_length=0;
      my $intergenic_map_length=0;      

      
      if($gff{$scfd_name}{gene}){
        # Sort orders of 'genes' according to the 'start' value;
        my (@sorted_gene_orders,@unsorted_gene_orders,@sorted_exon_orders,@unsorted_exon_orders);
        for my $gene_order(keys %{$gff{$scfd_name}{gene}}){
          push @unsorted_gene_orders, $gene_order;
        }
        for my $exon_order(keys %{$gff{$scfd_name}{exon}}){
          push @unsorted_exon_orders, $exon_order;
        }
        @sorted_gene_orders =sort {$gff{$scfd_name}{gene}{$a}{start} <=> $gff{$scfd_name}{gene}{$b}{start} } @unsorted_gene_orders;
        @sorted_exon_orders =sort {$gff{$scfd_name}{exon}{$a}{start} <=> $gff{$scfd_name}{exon}{$b}{start} } @unsorted_exon_orders;

        # Extract fragments' locs for %all_locs;
        my $scfd_loc_length=0;
        if($all_locs{$scfd_name}){
          my @tmp_locs=@{$all_locs{$scfd_name}};
          while(@tmp_locs){
            my $loc_front=shift @tmp_locs;
            my $loc_behind=shift @tmp_locs;
            # for my $one($front..$behind){
            my $length_loc=$loc_behind-$loc_front+1;
            $scfd_loc_length+=$length_loc;

            while (@sorted_gene_orders){
              my $gene_order=$sorted_gene_orders[0];
              my $gff_front=$gff{$scfd_name}{gene}{$gene_order}{start};
              my $gff_behind=$gff{$scfd_name}{gene}{$gene_order}{stop};
              my $length_gff=$gff_behind-$gff_front+1;
              if($length_gff == $length_loc){
                if($gff_behind <= $loc_behind){
                  shift @sorted_gene_orders; 
                }
                if($gff_front == $loc_front ){
                   $gene_map_length+=$length_gff;
                }elsif($loc_behind >= $gff_front && $gff_front > $loc_front ){
                          $gene_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front < $loc_front  && $loc_front <= $gff_behind){
                          $gene_map_length+=$gff_behind-$loc_front+1;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
              elsif($length_gff < $length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_gene_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $gene_map_length+=$length_gff;
                }elsif($gff_front< $loc_front && $loc_front<=$gff_behind ){
                    $gene_map_length+=$gff_behind-$loc_behind+1;
                }elsif($loc_front< $gff_front && $gff_front <=$loc_behind && $gff_behind > $loc_behind){
                    $gene_map_length+=$loc_behind-$gff_front+1;
                }elsif($loc_front<$gff_front && $gff_behind<$loc_behind){
                    $gene_map_length+=$length_gff;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
              elsif($length_gff>$length_loc){
                if($gff_behind <= $loc_behind){
                  shift @sorted_gene_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $gene_map_length+=$length_loc;
                }elsif($loc_front<$gff_front && $gff_front<=$loc_behind){
                    $gene_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front<$loc_front && $loc_front <=$gff_behind && $loc_behind> $gff_behind){
                    $gene_map_length+=$gff_behind-$loc_front+1;
                }elsif($loc_front>$gff_front && $loc_behind<$gff_behind){
                    $gene_map_length+=$length_loc;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
            }
                  
            
            while (@sorted_exon_orders){
              my $exon_order=$sorted_exon_orders[0];
              my $gff_front=$gff{$scfd_name}{exon}{$exon_order}{start};
              my $gff_behind=$gff{$scfd_name}{exon}{$exon_order}{stop};
              my $length_gff=$gff_behind-$gff_front+1;      
              if($length_gff == $length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_exon_orders;
                }
                if($gff_front == $loc_front ){
                    $exon_map_length+=$length_gff;
                }elsif($loc_behind >= $gff_front && $gff_front > $loc_front ){
                    $exon_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front < $loc_front && $loc_front <= $gff_behind){
                    $exon_map_length+=$gff_behind-$loc_front+1;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
              elsif($length_gff < $length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_exon_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $exon_map_length+=$length_gff;
                }elsif($gff_front< $loc_front && $loc_front <=$gff_behind ){
                    $exon_map_length+=$gff_behind-$loc_behind+1;
                }elsif($loc_front< $gff_front && $gff_front<=$loc_behind && $gff_behind > $loc_behind){
                    $exon_map_length+=$loc_behind-$gff_front+1;
                }elsif($loc_front<$gff_front && $gff_behind<$loc_behind){
                    $exon_map_length+=$length_gff;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }  
              elsif($length_gff>$length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_exon_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $exon_map_length+=$length_loc;
                }elsif($loc_front<$gff_front && $gff_front<=$loc_behind){
                    $exon_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front<$loc_front && $loc_front <=$gff_behind && $loc_behind> $gff_behind){
                    $exon_map_length+=$gff_behind-$loc_front+1;
                }elsif($loc_front>$gff_front && $loc_behind<$gff_behind){
                    $exon_map_length+=$length_loc;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }                       
            }
          }
        }
        $intergenic_map_length=$scfd_loc_length-$gene_map_length;
        $intron_map_length=$gene_map_length-$exon_map_length;
      }

    

      # Output part!!!
      my($intergenic_map_ratio,$gene_map_ratio,$exon_map_ratio,$intron_map_ratio)=qw(0 0 0 0);
      unless($intergenic_length==0){
        $intergenic_map_ratio=$intergenic_map_length/$intergenic_length;
      }
      else{$intergenic_map_ratio="N/A";
      }

      unless($gene_length==0){
        $gene_map_ratio=$gene_map_length/$gene_length;
      }
      else{$gene_map_ratio="N/A";
      }
         
      unless($exon_length==0){
        $exon_map_ratio=$exon_map_length/$exon_length;
      }
      else{$exon_map_ratio="N/A";
      }

      unless($intron_length==0){
        $intron_map_ratio=$intron_map_length/$intron_length;
      }
      else{$intron_map_ratio="N/A";
      }
      $all_intergenic_map_length+=$intergenic_map_length;
      $all_gene_map_length+=$gene_map_length;
      $all_exon_map_length+=$exon_map_length;
      $all_intron_map_length+=$intron_map_length;


      $coverage_ratio_fh->print("$scfd_name\t$intergenic_length\t$intergenic_map_length\t$intergenic_map_ratio\t$gene_length\t$gene_map_length\t$gene_map_ratio\t$exon_length\t$exon_map_length\t$exon_map_ratio\t$intron_length\t$intron_map_length\t$intron_map_ratio\n");
      print "Done!\n";
    }
  }
  my ($all_intergenic_map_ratio,$all_gene_map_ratio,$all_exon_map_ratio,$all_intron_map_ratio)=qw(0 0 0 0);
  if($all_gene_length==0){$all_gene_map_ratio='N/A';}
  else{
  $all_gene_map_ratio=$all_gene_map_length/$all_gene_length;}
  if($all_exon_length==0){$all_exon_map_ratio='N/A';}
  else{
  $all_exon_map_ratio=$all_exon_map_length/$all_exon_length;}
  if($all_intron_length==0){$all_intron_map_ratio='N/A';}
  else{
  $all_intron_map_ratio=$all_intron_map_length/$all_intron_length;}
  if($all_intergenic_length==0){$all_intergenic_map_ratio='N/A';}
  else{
  $all_intergenic_map_ratio=$all_intergenic_map_length/$all_intergenic_length;}
  
  $coverage_ratio_fh->print("Intergenic region map ratio is\t$all_intergenic_map_ratio\nGenes region map ratio is\t$all_gene_map_ratio\nExon region map ratio is\t$all_exon_map_ratio\nIntron region map ratio is\t$all_intron_map_ratio\n");
}




package RestrictionDigest::Single;

use 5.8.0;
use strict;
use warnings FATAL => 'all';

=head1 NAME


RestrictionDigest::Single - the simulation tool for singel enzyme digestion!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';



=head1 SYNOPSIS

RestrictionDigest::Single is used for the simulation of single-enzyme GBS/RAD approaches.


    use RestrictionDigest;

    $digest = RestrictionDigest::Single->new();
   
    $digest->add_ref($reference_file); # $reference_file  refers to the full path of the reference file;
    
    $digest->add_single_enzyme($enzyme); # $enzyme refers to the enzyme User wants to use in the digestion simulation;

    $digest->add_output_dir($output_dir);# $output_dir="A directory in which you want to put the result files";

    $species->single_digest();


=cut

use IO::File;
use File::Basename;

use vars qw(%fields %enzyme_list);

# Define the fields to be used as identifiers.
%fields = (
  ref     => undef,
  gff =>  undef,
  SNPs => undef,
  enzyme => undef,
  output_dir=>  undef,
  range_start => 201,
  range_end => 500,
  sequence_type => '100PE',
  lengths_distribution_start =>100,
  lengths_distribution_end  => 1000,
  lengths_distribution_step => 50,
);

# Make the endonuclease enzyme-container--a hash.
%enzyme_list=(
      EcoRI => 'G|AATTC',
      HinfI => 'G|ANTC',
      AluI => 'AG|CT',
      AvaII => 'G|GWCC',
      BamHI => 'G|GATCC',
      HindIII => 'A|AGCTT',
      MseI => 'T|TAA',
      MspI => 'C|CGG',
      NdeI => 'CA|TATG',
      PstI => 'CTGCA|G',
      SalI => 'G|TCGAC',
      XbaI => 'T|CTAGA',
      XhoI => 'C|TCGAG',
      NlaIII => 'NCATG|N',
      AatII =>  'GACGT|C',
      AccI  =>  'GT|MKAC',
      Acc65I  =>  'G|GTACC',
      AciI  =>  'C|CGC',
      AclI  => 'AA|CGTT',
      AfeI  => 'AGC|GCT',
      AflII   => 'C|TTAAGA',
      AflIII  => 'A|CRYGT',
      AgeI  => 'A|CCGGT',
      AlwNI  => 'CAGNNN|CTG',
      ApaI  => 'GGGCC|C',
      ApaLI   => 'G|TGCAC',
      ApeKI =>  'G|CWGC',
      ApoI  => 'R|AATTY',
      AscI  => 'GG|CGCGCC',
      AseI => 'AT|TAAT',
      AvaI => 'C|YCGRG',
      BclI =>  'T|GATCA',
      AvrII => 'C|CTAGG',
      BaeGI => 'GKGCM|C',
      BanI => 'G|GYRCC',
      BanII => 'GRGCY|C',
      BbvCI =>  'CC|TCAGC',
      BsaAI => 'YAC|GTR',
      BsaJI => 'C|CNNGG',
      BsaWI => 'W|CCGGW',
);

=head1 SUBROUTINES/METHODS

=head2 new
     
     $digest=RestrictionDigest::Single->new(); 
          Create a new object.

=cut

sub new {
  my $class = shift;
  my $self = {
    %fields,
  };
  bless ($self, $class);
  return $self;
}

=head2 add_ref
     
     $digest->add_ref(-reference=>'Full path of the reference file');
          Add the reference file  to the RestrictionDigest variable.

=cut

sub add_ref {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-reference$/){
      $self->{ref}=$parameters{$_};
    }
    else{
    die "Unacceptable parameters in the method <add_ref>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  } 
  
  # Check the format of reference
  my $ref_fh=IO::File->new("$self->{ref}",'r');
  my $line_num=0;
  print"Examing the format of the reference file...\t";
  while(<$ref_fh>){
    chomp;
    $line_num++;
    if($line_num % 2 != 0){
      unless($_=~/^>\S+$/){
        die "The format of reference file is not allowed, please convert it.\nYou can find the example format in the perldoc RestrictionDigest!\n";
      }
    }
    if($line_num % 2 ==0 ){
      unless($_=~/^[AGTCNRYMKSWHBVD]+$/i){
        die"The sequences of reference contain some charactors other than AGTCNRYMKSWHBVD, please correct that!\n";
      }
    }
  }
  print"the format is OK.\n";
  $ref_fh->close;
}

=head2 add_single_enzyme
     
     $digest->add_single_enzyme(-enzyme=>'EcoRI');
          Add the single enzyme's name used to digest the whole genome.
          
=cut


sub add_single_enzyme {
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-enzyme$/){
      $self->{enzyme}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters in the method <add_singel_enzyme>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
  
  # Check the existence of enzymes in the enzyme-container.
  my $enzyme_exists=0;
  foreach(keys %enzyme_list){
    if($_=~/^$self->{enzyme}$/i){
      $enzyme_exists++;
    }
  }
  if($enzyme_exists >0 ){
  }else{
    die "The single restrict endonuclease  User provides does not exist.\nHowever User can add this enzyme and its recognition sites to enzyme-container via the method 'new_enzyme'.\n";
  }
}

=head2 new_enzyme

     $digest->new_enzyme(-enzyme_name=>'EcoRI', -recognition_site=>'G|AATTC');
     
     If the enzyme User wants to use does not exists in the default enzyme-pool, the enzyme can be added to the enzyme-pool.
     
     For now, the enzymes in the default enzyme pool are  as follows:
           EcoRI,HinfI,AluI,AvaII,BamHI,HindIII,MseI,MspI,NdeI,PstI,SalI,XbaI,XhoI,NlaIII.
     
     Adding new enzyme(s) to the enzyme-pool is easy but the enzyme added to the pool is temporary. The added enzyme(s) is active
     only in the class context to which the 'new_method' belongs.

=cut

sub new_enzyme {
  my $self=shift;
  my %parameters=@_;
  my $new_enzyme_name;
  my $new_enzyme_site;
  for (keys %parameters){
    if($_=~/^-enzyme_name$/){
      $new_enzyme_name=$parameters{$_};
    }
    elsif($_=~/^-recognition_site$/){
      $new_enzyme_site=$parameters{$_};
      unless($new_enzyme_site=~/|/){
        die "The 'recognition_site' must contain a cut flag '|'. Please perldoc RestrictionDigest for example or read README for moer help\n";
      }
    }
    else{
      die"Unacceptable parameters in the method <new_enzyme>, please perldoc RestrictionDigest for example or read README for more help\n";
    }
  }
  $enzyme_list{$new_enzyme_name}=$new_enzyme_site;
}

=head2 change_range

     $digest->change_range(-start=>201, -end=>500);
          Usually, not all fragments are needed by the scientists. The fragments in particular range are more attractive to them.
          RestrictionDigest provides the sequences and the positions on the reference of these fragments in range. The default range is 
           [201bp, 500]. User can change this range via the method <change_range>. 

=cut

sub change_range {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-start$/){
      $self->{range_start}=$parameters{$_};
    }
    elsif($_=~/^-end$/){
      $self->{range_end}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters in the method <change_range>,please perldoc RestrictionDigest for example or read README for more help.\n";
    }
  }
  if($self->{range_start} < $self->{range_end} ) {
  }else{
    die "The first parameter of the method <change_range> must be smaller than the second parameter\n";
  }
}

=head2 change_lengths_distribution_parameters
    
     $digest->change_lengths_distribution_parameters(-front=>100,-behind=>1000,-step=>50);
        This program will give a simple summary of the lengths of all fragments after digestion. By default, three scopes, which are  
         <=100bp, 101bp-1000bp and >1000bp, are given. The number of the fragments whose length fallen in these scopes are given. More 
         fine sorted, the 101bp-1000bp scope is split by the step of 50bp. These three parameters: the front edge and the behind edge 
         of the scope and the step can be changed via the method 'change_lengths_distribution_parameters';
=cut

sub change_lengths_distribution_parameters {
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-front$/){
      $self->{lengths_distribution_start}=$parameters{$_};    
    }
    elsif($_=~/^-behind$/){
      $self->{lengths_distribution_end}=$parameters{$_};
    }
    elsif($_=~/^-step$/){
      $self->{lengths_distribution_step}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters in the method <change_lengths_distribution_parameters>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
}



=head2 add_output_dir
    
     $digest->add_output_dir(-output_dir=>'Full path of the output directory');
          Adds the directory into which the result files are put by the program. 
          The directory must be in the absolute path format. 

=cut

sub add_output_dir {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-output_dir$/){
      $self->{output_dir}=$parameters{$_};
    }
    else{
      die "Unacceptable parameter in the method <add_output_dir>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
  if(-d $self->{output_dir} ) {
  }else{
    die "The output directory User provides does not exist.\n";
  }
}

=head2 single_digest

     $digest->single_digest();
          Execute the digestion process. 
          This process will produce several result files which will be located in the output directory 
          added through the add_output_dir method. 
          This process will take some time. It depends on the size of the reference file.
          During the process, Prompting messages will be output to the STDOUT 
          like "Digestion: >ScaffoldName is done!".

=cut

sub single_digest {

  my $self=shift;
  my $ref=$self->{ref};
  if($ref){
    my $name_reference=basename($ref);
    print "The reference file User provides is:\t",$name_reference,"\n";
  }else{
    die "No reference file provied! Please perldoc RestrictionDigest for example.\n";
  }
  my $enzyme=$self->{enzyme};
  if($enzyme){
    print "The restrict endonuclease User provides is:\t", $enzyme, "\n";
  }
  else{
    die"No restrict endonuclease provided! Please perldoc RestrictionDigest for example.\n";
  }
  my $range1=$self->{range_start};
  my $range2=$self->{range_end};
  print "The seperate range of the fragments User wants to output is:\t$range1 to $range2 bp.\n";
  my $output_dir=$self->{output_dir};
  if($output_dir){
    print "The output directory User provides is:\t",$output_dir, "\n";
  }
  else{
    die "No output directory provided! Please perldoc RestrictionDigest for example.\n";
  }
  $output_dir=~s/^(.+)\/?$/$1/;
  
  # Make the file handle to the reference file.
  my $ref_fh=IO::File->new("$ref",'r');
  
  # Get the basename of reference,not include the path.
  my $name_reference=basename($ref);
  
  # Make the file handles of result files.
  my $all_seq_file_fh=IO::File->new(">>$output_dir/seq_frags_${name_reference}_by_${enzyme}");
  my $all_loc_file_fh=IO::File->new(">>$output_dir/position_frags_${name_reference}_by_${enzyme}");
  my $seq_in_range_file_fh=IO::File->new(">>$output_dir/seq_frags_in_range_${name_reference}_by_${enzyme}");
  my $loc_in_range_file_fh=IO::File->new(">>$output_dir/position_frags_in_range_${name_reference}_by_${enzyme}");
  my $reduced_ratio_every_scaffold_file_fh=IO::File->new(">>$output_dir/reduced_ratio_every_scaffold_${name_reference}_by_${enzyme}");
  my $summary_digestion_fh=IO::File->new(">>$output_dir/digestion_summary_${name_reference}_by_${enzyme}");
  
  # Get the base-compositions of the enzyme which User inputs.
  my $enzyme_seq;
  foreach(keys %enzyme_list){
    if($_=~/^$enzyme$/i){
      $enzyme_seq=$enzyme_list{$_};
    }
  }

  #get the relative location of cutting of the recognation sequences;
  my @enzyme_rec_seq=split //, $enzyme_seq;
  my @enzyme_cut_seq;
  my $enzyme_rec_start=0;
  my $enzyme_cut_loc;
  
  foreach(@enzyme_rec_seq){
    $enzyme_rec_start++;
    if($_=~/\|/){
      $enzyme_cut_loc=$enzyme_rec_start;
    }
    else{
      push @enzyme_cut_seq,$_;
    }
  }
  
  # Replace the degenerate bases of the enzyme;
  my $num_of_enzyme_cut_seq=@enzyme_cut_seq;
  for(@enzyme_cut_seq){
    if($_=~/N/){$_="[ATCG]";}
    if($_=~/R/){$_="[AG]";}
    if($_=~/Y/){$_="[CT]";}
    if($_=~/M/){$_="[AC]";}
    if($_=~/K/){$_="[GT]";}
    if($_=~/S/){$_="[GC]";}
    if($_=~/W/){$_="[AT]";}
    if($_=~/H/){$_="[ATC]";}
    if($_=~/B/){$_="[GTC]";}
    if($_=~/V/){$_="[GAC]";}
    if($_=~/D/){$_="[GAT]";}
  }
  my $seq_string_enzyme=join '', @enzyme_cut_seq;


  
  # Make the scalars and arrays used in the whole programme;
  my $rate_overall=0;
  my $rate_in_range=0;
  my @lengths_of_fragments=();
  my @lengths_fragments_in_range=();
  my $overall_fragments_length=0;
  my $overall_fragments_in_range_length=0;
  my $overall_length_of_scfd=0;

  ########## CIRCLES START HERE!###########

  # Process one scaffold one circle;
  my ($scaffold_name, $scaffold_seq);
  while(<$ref_fh>){
    chomp;
    if($_=~/^>/){
      $scaffold_name=$_;
    }
    else{
      print"Digesting $scaffold_name...\t";
        $scaffold_seq=$_;
        $reduced_ratio_every_scaffold_file_fh->print("$scaffold_name\t");
        my $length_of_scafd=length $scaffold_seq;
        # Set the array which contains cut locations of the enzyme.
        my @enzyme_locs=();
        # Find locs of recognition of the enzyme in this scaffold;the locs are meaningful in this array; not in the human context; 
        my $seq_for_search=$scaffold_seq;
        my $string=$seq_string_enzyme;
        # Find the recognition sites of this enzyme in this scaffold.
          if($string=~/^(\w+)$/){
            my $loc_in_array;
            $loc_in_array=index($seq_for_search,$string);
            unless($loc_in_array==-1){
              push @enzyme_locs, $loc_in_array;
            }
            while($loc_in_array!=-1){
              $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
              unless($loc_in_array==-1){
                push @enzyme_locs, $loc_in_array;
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $back_bases=$3;
            my $generate_bases=$2;
            my @generate_bases=split //, $generate_bases;
            my $seq_for_search=$scaffold_seq;
            for my $one_base(@generate_bases){
              $string=$breast_bases.$one_base.$back_bases;
              my $loc_in_array;
              $loc_in_array=index($seq_for_search,$string);
              unless($loc_in_array==-1){
                push @enzyme_locs, $loc_in_array;
              }
              while($loc_in_array!=-1){
                $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                unless($loc_in_array==-1){
                  push @enzyme_locs, $loc_in_array;
                }
              }
            }
          }

          elsif($string=~/^\[(\w+)\](\w+)\[(\w+)\]$/){
            my $g1_bases=$1;
            my $middle=$2;
            my $g2_bases=$3;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                $string=$one_base.$middle.$two_base;
                my $loc_in_array;
                  $loc_in_array=index($seq_for_search,$string);
                  unless($loc_in_array==-1){
                    push @enzyme_locs, $loc_in_array;
                  }
                  while($loc_in_array!=-1){
                    $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                    unless($loc_in_array==-1){
                      push @enzyme_locs, $loc_in_array;
                    }
                  }                
              }
            }
          }

          elsif($string=~/(\w+)\[(\w+)\](\w+)\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $middle=$3;
            my $back_bases=$5;
            my $g2_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                $string=$breast_bases.$one_base.$middle.$two_base.$back_bases;
                  my $loc_in_array;
                  $loc_in_array=index($seq_for_search,$string);
                  unless($loc_in_array==-1){
                    push @enzyme_locs, $loc_in_array;
                  }
                  while($loc_in_array!=-1){
                    $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                    unless($loc_in_array==-1){
                      push @enzyme_locs, $loc_in_array;
                    }
                  }
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $g2_bases=$3;
            my $back_bases=$5;
            my $g3_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            my @g3_bases=split //, $g3_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                for my $three_base(@g3_bases){
                  $string=$breast_bases.$one_base.$two_base.$three_base.$back_bases;
                    my $loc_in_array;
                    $loc_in_array=index($seq_for_search,$string);
                    unless($loc_in_array==-1){
                      push @enzyme_locs, $loc_in_array;
                    }
                    while($loc_in_array!=-1){
                      $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                      unless($loc_in_array==-1){
                        push @enzyme_locs, $loc_in_array;
                      }
                    }
                }
              }
            }
          }
          elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\](\w+)/){
            my $breast_bases=$1;
            my $g1_bases=$2;
            my $g2_bases=$3;
            my $back_bases=$4;
            my @g1_bases=split //, $g1_bases;
            my @g2_bases=split //, $g2_bases;
            for my $one_base(@g1_bases){
              for my $two_base(@g2_bases){
                  $string=$breast_bases.$one_base.$two_base.$back_bases;
                    my $loc_in_array;
                    $loc_in_array=index($seq_for_search,$string);
                    unless($loc_in_array==-1){
                      push @enzyme_locs, $loc_in_array;
                    }
                    while($loc_in_array!=-1){
                      $loc_in_array=index($seq_for_search, $string, $loc_in_array+1);
                      unless($loc_in_array==-1){
                        push @enzyme_locs, $loc_in_array;
                      }
                    }
              }
            }
          }
        
        # Sort the array.
        @enzyme_locs=sort {$a<=>$b} @enzyme_locs;  

        
        # Produce the fragments' loc-pair.
        my @fragments_loc_pair;
        my ($front_rec_loc, $behind_rec_loc)=qw(0 0);
        if(@enzyme_locs) {
          while(@enzyme_locs){
            $behind_rec_loc=shift @enzyme_locs;
            # Push the front end loc;
            if($front_rec_loc==0){
              push @fragments_loc_pair, $front_rec_loc;
            }
            else{
              push @fragments_loc_pair, $front_rec_loc+$enzyme_cut_loc-1;
            }
            # Push the behind end loc;
            push @fragments_loc_pair, $behind_rec_loc+$enzyme_cut_loc-2;
            
            # Check whether the @enzyme_locs is empty.
            if(@enzyme_locs){
              $front_rec_loc=$behind_rec_loc;
            }
            else{
              push @fragments_loc_pair, $behind_rec_loc+$enzyme_cut_loc-1;
              push @fragments_loc_pair, $length_of_scafd-1;
            }
          }
        }

        # Output the fragments' locs of this scaffold.        
        
        # Define the scalars which contains the length of fragments;
        my $nums_of_fragments=@fragments_loc_pair/2;
        my $length_of_fragments=0;
        my $length_of_fragments_in_range=0;
        my $count_all=0;
        my $count_in_range=0;
        my ($before_cut_coorr,$after_cut_coorr,$bf_ct_crr_hm_rd,$af_ct_crr_hm_rd);


          # Attentions here!!! the output sequences of fragments are in accordance with the results of illumian sequencing.
          # For removing the sequencing primers, read1 and read2 can align with the fragments here directly!!! 
          # That is to say, fragments here contain the restriction overhangs!!!
          
          while(@fragments_loc_pair){
            $count_all++;
            my $front_end_loc=shift @fragments_loc_pair;
            my $behind_end_loc=shift @fragments_loc_pair;
            my $length_of_fragment=$behind_end_loc-$front_end_loc+1;
            my $seq_selected=substr($scaffold_seq, $front_end_loc, $length_of_fragment);         
              
          # Output the sequence and fragment positions of all-fragments to files;
          $all_seq_file_fh->print("$scaffold_name-$count_all\n$seq_selected\n");
          my $front_loc_on_scaffold=$front_end_loc+1;
          my $behind_loc_on_scaffold=$behind_end_loc+1;
        
          $all_loc_file_fh->print("$scaffold_name-$count_all\t$front_loc_on_scaffold\t$behind_loc_on_scaffold\t$length_of_fragment\n");
          push @lengths_of_fragments, $length_of_fragment;
                    
          # We treat the fragments in range in the same way of all fragments as above;
          if($length_of_fragment>=$range1){
            if($length_of_fragment<=$range2){
              $count_in_range++;
              $seq_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_selected\n");
              $loc_in_range_file_fh->print("$scaffold_name-$count_all\t$front_loc_on_scaffold\t$behind_loc_on_scaffold\t$length_of_fragment\n");
              push @lengths_fragments_in_range,$length_of_fragment;
              $length_of_fragments_in_range+=$length_of_fragment;
            }
          }
          
          # Count all fragments length no matter in range or not!
          $length_of_fragments+=$length_of_fragment;
        }


        my $reduced_rate=$length_of_fragments/$length_of_scafd;
        my $fragment_in_range_rate_scafd=$length_of_fragments_in_range/$length_of_scafd;
        $overall_fragments_length+=$length_of_fragments;
        $overall_fragments_in_range_length+=$length_of_fragments_in_range;
        $overall_length_of_scfd+=$length_of_scafd;
        $reduced_ratio_every_scaffold_file_fh->print("$reduced_rate\t$fragment_in_range_rate_scafd\n");
        print"Done!\n"; 
    }
  }
  my $length_ratio_all_frags=$overall_fragments_length/$overall_length_of_scfd;
  my $length_ratio_frags_in_range=$overall_fragments_in_range_length/$overall_length_of_scfd;
  my $num_all_frags=@lengths_of_fragments;
  my $num_frags_in_range=@lengths_fragments_in_range;
  my $lengths_fragments_in_range=\@lengths_fragments_in_range;


  # Output the summary of digestion.
  $summary_digestion_fh->print("Number of all fragments\t$num_all_frags\n");
  $summary_digestion_fh->print("Number of fragments in range\t$num_frags_in_range\n");
  $summary_digestion_fh->print("Ratio after reducing of all fragments\t$length_ratio_all_frags\n");
  $summary_digestion_fh->print("Ratio after reducing of fragments in range\t$length_ratio_frags_in_range\n");

  my $lengths_distribution_start=$self->{lengths_distribution_start};
  my $lengths_distribution_end=$self->{lengths_distribution_end};
  my $lengths_distribution_step=$self->{lengths_distribution_step};
  my $lengths_distribution_tmp;
  my (@starts,@ends,@counts,%lengths_distribution);

  my $num_pair=0;

  for($lengths_distribution_tmp= $lengths_distribution_start;$lengths_distribution_tmp<$lengths_distribution_end+$lengths_distribution_step;$lengths_distribution_tmp=$lengths_distribution_tmp+$lengths_distribution_step){
    if($lengths_distribution_tmp== $lengths_distribution_start){
      $num_pair++;
      push @starts, '0';
      push @ends , $lengths_distribution_tmp;
      $lengths_distribution{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;

      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      $lengths_distribution{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      
    }
    elsif($lengths_distribution_tmp  < $lengths_distribution_end ){
      if($lengths_distribution_tmp+$lengths_distribution_step <$lengths_distribution_end ){
        $num_pair++;
        my $tmp_start=$lengths_distribution_tmp+1;
        my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
        push @starts, $tmp_start;
        push @ends, $tmp_end;
        $lengths_distribution{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        
    }
    elsif($lengths_distribution_tmp +$lengths_distribution_step >= $lengths_distribution_end){
      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_end;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      $lengths_distribution{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;

      $num_pair++;
      $tmp_start=$tmp_end+1;
      
      push @starts, $tmp_start;
      push @ends, 'bigger';
      $lengths_distribution{$num_pair}{"${tmp_start}bp-bigger bp"}=0;
      
      }
    }
  }
  


  for (@lengths_of_fragments){
    my $length=$_;
    my $num_pair=1;
    my $num_last_pair=@starts;
    for($num_pair=1;$num_pair<=$num_last_pair;$num_pair++){
      
      my $array_index=$num_pair-1;
      my $start=$starts[$array_index];
      my $end=$ends[$array_index];
      if($num_pair<$num_last_pair){
        if($length>=$start && $length<=$end){
          for my $description (keys %{$lengths_distribution{$num_pair}}){
            $lengths_distribution{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
      elsif($num_pair==$num_last_pair){
        if($length>=$start){
          for my $description (keys %{$lengths_distribution{$num_pair}}){
            $lengths_distribution{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
    }
  }
  
  $summary_digestion_fh->print("\nLengths\' scope\tNumber of fragments in this scope\tRatio of these fragments in all\n");
  
  for(1..@starts){
    my $num_pair=$_;
    for my $description(keys %{$lengths_distribution{$num_pair}}){
      my $num_fragments=$lengths_distribution{$num_pair}{$description};
      my $ratio;
      if($num_all_frags>0){
        $ratio=$num_fragments/$num_all_frags;
      }
      else{
        $ratio='N/A';
      }
      $summary_digestion_fh->print("$description\t$num_fragments\t$ratio\n");
    }
  }
  $summary_digestion_fh->print("\n");
}

=head2 add_SNPs
  
    $digest->add_SNPs(-SNPs=>'Full path of the SNPs file');
        Adds the absolute path of the SNPs file to RestrictionDigest variable.
=cut

sub add_SNPs {
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-SNPs$/){
      $self->{SNPs}=$parameters{$_};
    }
    else{
      die "Unacceptable parameters of the method <add_SNPs>, please perldoc RestrictionDigest for example or read README for more help\n";
    }
  }
}

=head2 count_SNPs_at_fragments
  
    $digest->count_SNPs_at_fragments();
        Count the expected SNPs appeared in the framgents.
=cut

sub count_SNPs_at_fragments {
  my $self=shift;
  my %parameters=@_;
  
  if(keys %parameters){
    for(keys %parameters){
      if($_=~/^-sequence_type$/){
        $self->{sequence_type}=$parameters{$_};
      }
      else{
        die "Unacceptable parameters in the method <count_SNPs_at_fragments>, please perldoc RestrictionDigest for example or read README for more help\n";
      }
    }
  }
  
  my $ref=$self->{ref};
  my $enzyme=$self->{enzyme};
  my $range1=$self->{range_start};
  my $SNPs=$self->{SNPs};
  if($SNPs){
    my $SNPs_basename=basename($SNPs);
    print "The SNPs file provided is\t", $SNPs_basename, ".\n"; 
  }
  else{
    die"No SNPs file provided! Please perldoc RestrictionDigest for example.\n";
  }
  my $sequence_type=$self->{sequence_type};

  # Check the sequencing type and determine the domains of the fragments to calculate the SNPs number.
  my ($pe_sequence_length, $se_sequence_length);
  my $sequence_type_code;
  if($sequence_type=~/^(\d+)pe$/i){
    $sequence_type_code=2;
    my $sequence_length=$1;
    unless($sequence_length*2 <= $range1){
      die"Unacceptable sequence_type when range's front border is $range1 bp. Please perldoc RestrictionDigest for example or read README for more help!\n";
    }
    else{
      $pe_sequence_length=$sequence_length;
    }
  }
  elsif($sequence_type=~/^(\d+)se$/i){
    $sequence_type_code=1;
    my $sequence_length=$1;
    unless($sequence_length<=$range1){
      die "Unacceptable sequence_type when range's front border is $range1 bp. Please perldoc RestrictionDigest for example or read README for more help!\n";
    }
    else{
      $se_sequence_length=$sequence_length;
    }
  }
  else{
    die "Unacceptable sequence_type. Please perldoc RestrictionDigest for example or read README for more help!\n";
  }

  print "The sequencing type provided is\t", $sequence_type, ".\n";

  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  # Get the basename of reference,not include the path.
  my $name_reference=basename($ref);

  # Make the file handle to the summary file;
  my $summary_digestion_fh=IO::File->new(">>$output_dir/digestion_summary_${name_reference}_by_${enzyme}");

  # Make the file handle to the locs file of all fragments and fragments in range.
  my $all_loc_file_fh=IO::File->new("$output_dir/position_frags_${name_reference}_by_${enzyme}",'r');
  my $loc_in_range_file_fh=IO::File->new("$output_dir/position_frags_in_range_${name_reference}_by_${enzyme}",'r');
  
  # Make the file handle to the SNPs file;
  my $SNPs_fh=IO::File->new("$SNPs",'r');
  my %SNPs;
  while(<$SNPs_fh>){
    chomp;
    my ($SNPs_scaffold_name,$SNPs_pos,$type)=split /\s+/, $_;
    push @{$SNPs{$SNPs_scaffold_name}}, $SNPs_pos;
  }

  # Create the SNPs file to hold SNPs of all fragments and fragments in range when sequencing through pair-end.
  my($SNPs_at_all_frags_fh,$SNPs_at_frags_in_range_fh);
  if($sequence_type_code==2){
    $SNPs_at_all_frags_fh=IO::File->new(">>$output_dir/SNPs_at_all_frags_${name_reference}_by_${enzyme}");
    $SNPs_at_frags_in_range_fh=IO::File->new(">>$output_dir/SNPs_at_frags_in_range_${name_reference}_by_${enzyme}");
  }

  print "Counting SNPs at the output fragments...\t";
  # Scan all locs file to count SNPs and output the result to the summary file.
  my $all_SNPs_number=0;
  if($sequence_type_code == 1){
    while(<$all_loc_file_fh>){
      chomp;
      my($fragment_name,$pos1,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        for my $SNP (@{$SNPs{$scaffold_name}}){
          if( ($SNP>=$pos1 && $SNP<=($pos1+$se_sequence_length-1)) || ( $SNP>= ($pos2-$se_sequence_length+1) && $SNP<=$pos2 ) ){
          $all_SNPs_number++;
          }
        }
      }
    }
    $all_SNPs_number=int($all_SNPs_number/2);
    $summary_digestion_fh->print("Expected SNPs at all fragments via $sequence_type:\t$all_SNPs_number\n");
  }
  elsif($sequence_type_code ==2){
    while(<$all_loc_file_fh>){
      chomp;
      my($fragment_name,$pos1,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        for my $SNP (@{$SNPs{$scaffold_name}}){
          if( ($SNP>=$pos1 && $SNP<=($pos1+$pe_sequence_length-1)) || ( $SNP>= ($pos2-$pe_sequence_length+1) && $SNP<=$pos2 ) ){
          $all_SNPs_number++;
          $SNPs_at_all_frags_fh->print("$scaffold_name\t$SNP\n");
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected SNPs at all fragments via $sequence_type:\t$all_SNPs_number\n");
  }
  # Scan locs in range file to count SNPs and output the result to the summary file.

  my $SNPs_in_range_number=0;

  if($sequence_type_code == 1){
    while(<$loc_in_range_file_fh>){
      chomp;
      my($fragment_name,$pos1,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        for my $SNP (@{$SNPs{$scaffold_name}}){
          if( ($SNP>=$pos1 && $SNP<=($pos1+$se_sequence_length-1)) || ( $SNP>= ($pos2-$se_sequence_length+1) && $SNP<=$pos2 ) ){
          $SNPs_in_range_number++;
          }
        }
      }
    }
    $SNPs_in_range_number=int($SNPs_in_range_number/2);
    $summary_digestion_fh->print("Expected SNPs at fragments in range via $sequence_type:\t$SNPs_in_range_number\n");
  }
  elsif($sequence_type_code ==2){
    while(<$loc_in_range_file_fh>){
      chomp;
      my($fragment_name,$pos1,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        for my $SNP (@{$SNPs{$scaffold_name}}){
          if( ($SNP>=$pos1 && $SNP<=($pos1+$pe_sequence_length-1)) || ( $SNP>= ($pos2-$pe_sequence_length+1) && $SNP<=$pos2 ) ){
          $SNPs_in_range_number++;
          $SNPs_at_frags_in_range_fh->print("$scaffold_name\t$SNP\n");
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected SNPs at fragments in range via $sequence_type:\t$SNPs_in_range_number\n");
  }
  print "Done!\n";
}

=head2 add_gff

    $digest->add_gff(-ref=>'Full path of the gff file');
        Adds the absolute path of gff file to RestrictionDigest variable. 
        If User want to calculate the coverage ratio of different genomic strutures by the digested 
        fragments, this method is obligatory.
        
=cut

sub add_gff { 
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-gff$/){
      $self->{gff}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters of the method <add_gff>, please perldoc RestrictionDigest for example or read README for more help\n";
    }
  }
}

=head2 all_frags_coverage_ratio

     $digest->all_frags_coverage_ratio();
          Calculate the coverage ratio of different parts of the genome reference, 
          including mRNA, CDS, UTR and intergenic region, by all digested fragments.
          This process will take a very long time, please be patient.  
=cut

sub all_frags_coverage_ratio {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $enzyme=$self->{enzyme};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $loc_file_all_frags="$output_dir/position_frags_${basename_ref}_by_${enzyme}";

  $self->genome_structure_coverage_ratio($loc_file_all_frags);
}

=head2 frags_in_range_coverage_ratio

     $digest->frags_in_range_coverage_ratio();
          Calculate the coverage ratio of different parts of the genome reference, 
          including mRNA, CDS, UTR and intergenic region, by the fragments in range.
          This process will take a very long time, please be patient. 
=cut

sub frags_in_range_coverage_ratio {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $enzyme=$self->{enzyme};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $loc_file_all_frags="$output_dir/position_frags_in_range_${basename_ref}_by_${enzyme}";

  $self->genome_structure_coverage_ratio($loc_file_all_frags);
}

=head2  genome_structure_coverage_ratio

        A subroutine invoked by the 'all_frags_coverage_ratio' and 
        the 'frags_in_range_coverage_ratio' methods. 
        User do not use this subroutine.

=cut

 
sub genome_structure_coverage_ratio {
  my $self=shift;
  my $loc_frags_file=shift;

  # Get the reference name, enzyme names from hash.
  my $ref=$self->{ref};
  my $ref_fh=IO::File->new($ref,'r');

  my $basename_ref=basename($ref);
  my $enzyme=$self->{enzyme};
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  # Get the full path and name of gff file and make a filehandle to it.
  my $gff=$self->{gff};
  
  if($gff){
    my $gff_basename=basename($gff);
    print "The GFF file provided is\t", $gff_basename, ".\n";
  }
  else{
    die"No GFF file provided! Please perldoc RestrictionDigest for example.\n";
  }

   
  # Make  output filehandles to generic region and intergenic region coverage ratio files. 
  my $coverage_ratio_fh;
  
  my $basename_loc_frags_file=basename($loc_frags_file);

  if($basename_loc_frags_file =~/in_range/){
    $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_ratio_in_range_${basename_ref}_by_${enzyme}");
  }
  else{
    $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_ratio_${basename_ref}_by_${enzyme}");
  }

  $coverage_ratio_fh->print("ScaffoldName\tIntergenicLength\tIntergenicMapLength\tIntergenicMapRatio\tGenesLength\tGenesMapLength\tGenesMapRatio\tExonsLength\tExonsMapLength\tExonsMapRatio\tIntronsLength\tIntronsMapLength\tIntronsMapRatio\n");


  # Set variables used in this program.
  my($all_intergenic_length,$all_intergenic_map_length,$all_gene_length,$all_gene_map_length,$all_exon_length,$all_exon_map_length,$all_intron_length,$all_intron_map_length)=qw(0 0 0 0 0 0 0 0);

  # Fetch the lengths of different parts of this scaffold from the GFF file.
      
      my $gff_fh=IO::File->new($gff,'r');
      my %gff; my %exon_check; my %gene_check;
      my ($gff_gene_count,$gff_exon_count)=qw(0 0);
      print "Scanning the GFF file...\t";
      while(<$gff_fh>){
        chomp;
        unless($_=~/^#/){
          my @gff=split /\t/, $_;
          if(@gff == 9){
           my ($gff_scfd_name,$gff_type,$gff_start,$gff_stop)=@gff[0,2,3,4];
              if($gff_type=~/CDS|exon/){
                $gff_type='exon';
                unless($exon_check{$gff_scfd_name}{start}{$gff_start}){
                  $gff_exon_count++;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_exon_count}{start}=$gff_start;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_exon_count}{stop}=$gff_stop;
                  $exon_check{$gff_scfd_name}{start}{$gff_start}=1;
                }
              }
              elsif($gff_type=~/mRNA|gene/){
                $gff_type='gene';
                unless($gene_check{$gff_scfd_name}{start}{$gff_start}){
                  $gff_gene_count++;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_gene_count}{start}=$gff_start;
                  $gff{$gff_scfd_name}{$gff_type}{$gff_gene_count}{stop}=$gff_stop;
                  $gene_check{$gff_scfd_name}{start}{$gff_start}=1;
                }
              }
          }
        }
      }
      $gff_fh->close;
      print"Done!\n";

      # Scan the loc file and fetch the information.
      
      my $loc_frags_file_fh=IO::File->new($loc_frags_file,'r');
      my %all_locs;
      print "Scanning the positions file...\t";
      while(<$loc_frags_file_fh>){
        chomp;
        my($frgt_name,$start_coor,$stop_coor,$length)=split /\t/,$_;
        my $scfd_frgt;
        if($frgt_name=~/^>(.+)-\d+$/){
          $scfd_frgt=$1;
          push @{$all_locs{$scfd_frgt}}, ($start_coor, $stop_coor);
        }
      }
      $loc_frags_file_fh->close;
      print"Done!\n";



  # Handle every scaffold one by one.
  my($scfd_name);
  while(<$ref_fh>){
    chomp;
    if($_=~/^>(\S+)$/){
    $scfd_name=$1;
    }
    else{
      my $scaffold_length=length $_;
      print "CoverageRatio: Processing $scfd_name...\t";


      # Calculate the lengths of different parts of this scaffold.
      my ($exon_length,$intron_length,$gene_length,$intergenic_length)=qw(0 0 0 0);
      
      if($gff{$scfd_name}{gene}){
        for my $gene_order(keys %{$gff{$scfd_name}{gene}}){
          my $front=$gff{$scfd_name}{gene}{$gene_order}{start};
          my $behind=$gff{$scfd_name}{gene}{$gene_order}{stop};
          $gene_length+=$behind-$front+1;
        }
        $intergenic_length=$scaffold_length-$gene_length;
      }
      else{
        $intergenic_length=$scaffold_length;
      }
        
      if($gff{$scfd_name}{exon}){
        for my $exon_order(keys %{$gff{$scfd_name}{exon}}){
          my $front=$gff{$scfd_name}{exon}{$exon_order}{start};
          my $behind=$gff{$scfd_name}{exon}{$exon_order}{stop};
          $exon_length+=$behind-$front+1;
        }
        $intron_length=$gene_length-$exon_length;
      }

      $all_intergenic_length+=$intergenic_length;
      $all_gene_length+=$gene_length;
      $all_exon_length+=$exon_length;
      $all_intron_length+=$intron_length;     



      # Scan the locs file and store information about this scaffold.

      my $gene_map_length=0;
      my $exon_map_length=0;
      my $intron_map_length=0;
      my $intergenic_map_length=0; 

      if($gff{$scfd_name}{gene}){
        # Sort orders of 'genes' according to the 'start' value;
        my (@sorted_gene_orders,@unsorted_gene_orders,@sorted_exon_orders,@unsorted_exon_orders);
        for my $gene_order(keys %{$gff{$scfd_name}{gene}}){
          push @unsorted_gene_orders, $gene_order;
        }
        for my $exon_order(keys %{$gff{$scfd_name}{exon}}){
          push @unsorted_exon_orders, $exon_order;
        }
        @sorted_gene_orders =sort {$gff{$scfd_name}{gene}{$a}{start} <=> $gff{$scfd_name}{gene}{$b}{start} } @unsorted_gene_orders;
        @sorted_exon_orders =sort {$gff{$scfd_name}{exon}{$a}{start} <=> $gff{$scfd_name}{exon}{$b}{start} } @unsorted_exon_orders;


        # Extract fragments' locs for %all_locs;
        my $scfd_loc_length=0;

        if($all_locs{$scfd_name}){
          my @tmp_locs=@{$all_locs{$scfd_name}};
          while(@tmp_locs){
            my $loc_front=shift @tmp_locs;
            my $loc_behind=shift @tmp_locs;
            my $length_loc=$loc_behind-$loc_front+1;
            $scfd_loc_length+=$length_loc;

            while (@sorted_gene_orders){
              my $gene_order=$sorted_gene_orders[0];
              my $gff_front=$gff{$scfd_name}{gene}{$gene_order}{start};
              my $gff_behind=$gff{$scfd_name}{gene}{$gene_order}{stop};
              my $length_gff=$gff_behind-$gff_front+1;
              if($length_gff == $length_loc){
                if($gff_behind <= $loc_behind){
                  shift @sorted_gene_orders; 
                }
                if($gff_front == $loc_front ){
                   $gene_map_length+=$length_gff;
                }elsif($loc_behind >= $gff_front && $gff_front > $loc_front ){
                          $gene_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front < $loc_front  && $loc_front <= $gff_behind){
                          $gene_map_length+=$gff_behind-$loc_front+1;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
              elsif($length_gff < $length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_gene_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $gene_map_length+=$length_gff;
                }elsif($gff_front< $loc_front && $loc_front<=$gff_behind ){
                    $gene_map_length+=$gff_behind-$loc_behind+1;
                }elsif($loc_front< $gff_front && $gff_front <=$loc_behind && $gff_behind > $loc_behind){
                    $gene_map_length+=$loc_behind-$gff_front+1;
                }elsif($loc_front<$gff_front && $gff_behind<$loc_behind){
                    $gene_map_length+=$length_gff;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
              elsif($length_gff>$length_loc){
                if($gff_behind <= $loc_behind){
                  shift @sorted_gene_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $gene_map_length+=$length_loc;
                }elsif($loc_front<$gff_front && $gff_front<=$loc_behind){
                    $gene_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front<$loc_front && $loc_front <=$gff_behind && $loc_behind> $gff_behind){
                    $gene_map_length+=$gff_behind-$loc_front+1;
                }elsif($loc_front>$gff_front && $loc_behind<$gff_behind){
                    $gene_map_length+=$length_loc;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
            }

            while (@sorted_exon_orders){
              my $exon_order=$sorted_exon_orders[0];
              my $gff_front=$gff{$scfd_name}{exon}{$exon_order}{start};
              my $gff_behind=$gff{$scfd_name}{exon}{$exon_order}{stop};
              my $length_gff=$gff_behind-$gff_front+1;      
              if($length_gff == $length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_exon_orders;
                }
                if($gff_front == $loc_front ){
                    $exon_map_length+=$length_gff;
                }elsif($loc_behind >= $gff_front && $gff_front > $loc_front ){
                    $exon_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front < $loc_front && $loc_front <= $gff_behind){
                    $exon_map_length+=$gff_behind-$loc_front+1;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }
              elsif($length_gff < $length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_exon_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $exon_map_length+=$length_gff;
                }elsif($gff_front< $loc_front && $loc_front <=$gff_behind ){
                    $exon_map_length+=$gff_behind-$loc_behind+1;
                }elsif($loc_front< $gff_front && $gff_front<=$loc_behind && $gff_behind > $loc_behind){
                    $exon_map_length+=$loc_behind-$gff_front+1;
                }elsif($loc_front<$gff_front && $gff_behind<$loc_behind){
                    $exon_map_length+=$length_gff;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }  
              elsif($length_gff>$length_loc){
                if($gff_behind <= $loc_behind){
                    shift @sorted_exon_orders;
                }
                if($gff_front==$loc_front || $gff_behind==$loc_behind){
                    $exon_map_length+=$length_loc;
                }elsif($loc_front<$gff_front && $gff_front<=$loc_behind){
                    $exon_map_length+=$loc_behind-$gff_front+1;
                }elsif($gff_front<$loc_front && $loc_front <=$gff_behind && $loc_behind> $gff_behind){
                    $exon_map_length+=$gff_behind-$loc_front+1;
                }elsif($loc_front>$gff_front && $loc_behind<$gff_behind){
                    $exon_map_length+=$length_loc;
                }
                if($loc_behind< $gff_behind ){
                  last;
                }
              }                       
            }      
          }
        }
        $intergenic_map_length=$scfd_loc_length-$gene_map_length;
        $intron_map_length=$gene_map_length-$exon_map_length;
      }

      # Output part!!!
      my($intergenic_map_ratio,$gene_map_ratio,$exon_map_ratio,$intron_map_ratio)=qw(0 0 0 0);
      unless($intergenic_length==0){
        $intergenic_map_ratio=$intergenic_map_length/$intergenic_length;
      }
      else{$intergenic_map_ratio="N/A";
      }

      unless($gene_length==0){
        $gene_map_ratio=$gene_map_length/$gene_length;
      }
      else{$gene_map_ratio="N/A";
      }
         
      unless($exon_length==0){
        $exon_map_ratio=$exon_map_length/$exon_length;
      }
      else{$exon_map_ratio="N/A";
      }

      unless($intron_length==0){
        $intron_map_ratio=$intron_map_length/$intron_length;
      }
      else{$intron_map_ratio="N/A";
      }
      $all_intergenic_map_length+=$intergenic_map_length;
      $all_gene_map_length+=$gene_map_length;
      $all_exon_map_length+=$exon_map_length;
      $all_intron_map_length+=$intron_map_length;


      $coverage_ratio_fh->print("$scfd_name\t$intergenic_length\t$intergenic_map_length\t$intergenic_map_ratio\t$gene_length\t$gene_map_length\t$gene_map_ratio\t$exon_length\t$exon_map_length\t$exon_map_ratio\t$intron_length\t$intron_map_length\t$intron_map_ratio\n");
      print "Done!\n";
    }
  }
  my ($all_intergenic_map_ratio,$all_gene_map_ratio,$all_exon_map_ratio,$all_intron_map_ratio)=qw(0 0 0 0);
  if($all_gene_length==0){$all_gene_map_ratio='N/A';}
  else{
  $all_gene_map_ratio=$all_gene_map_length/$all_gene_length;}
  if($all_exon_length==0){$all_exon_map_ratio='N/A';}
  else{
  $all_exon_map_ratio=$all_exon_map_length/$all_exon_length;}
  if($all_intron_length==0){$all_intron_map_ratio='N/A';}
  else{
  $all_intron_map_ratio=$all_intron_map_length/$all_intron_length;}
  if($all_intergenic_length==0){$all_intergenic_map_ratio='N/A';}
  else{
  $all_intergenic_map_ratio=$all_intergenic_map_length/$all_intergenic_length;}
  
  $coverage_ratio_fh->print("Intergenic region map ratio is\t$all_intergenic_map_ratio\nGenes region map ratio is\t$all_gene_map_ratio\nExon region map ratio is\t$all_exon_map_ratio\nIntron region map ratio is\t$all_intron_map_ratio\n");
}



=head1 AUTHOR

Jinpeng Wang, Li Li, Guofan Zhang, C<< <RestrictionDigest at 163.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-RestrictionDigest at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=RestrictionDigest>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc RestrictionDigest


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=RestrictionDigest>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/RestrictionDigest>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/RestrictionDigest>

=item * Search CPAN

L<http://search.cpan.org/dist/RestrictionDigest/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Jinpeng Wang, Li Li, Guofan Zhang.

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


=cut

1; # End of RestrictionDigest
