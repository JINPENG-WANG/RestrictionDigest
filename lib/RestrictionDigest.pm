
package RestrictionDigest;

use 5.8.0;
use strict;
use warnings FATAL => 'all';



=head1 NAME

RestrictionDigest - A powerful Perl module for simulating genomic restriction digests!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 DESCRIPTION

Restriction site Associated DNA sequence (RAD), Genotyping by Sequencing (GBS) and their variations are popular methods
in reduced-representation sequencing. A key factor of these methods is to determine the enzyme or enzyme-combination to
use. Prior to practical experiments, simulating the digestions is useful. RestrictionDigest can simulate these digestions
and generate useful information of the evaluation which would help users to determine the usability of the enzyme evaluated
in reduced-representation library construction.

=cut





package RestrictionDigest::SingleItem::Double;

use 5.8.0;
use strict;
use warnings FATAL => 'all';

=head1 NAME

RestrictionDigest::SingleItem::Double

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';



=head1 SYNOPSIS

RestrictionDigest::SingleItem::Double is used for the simulation of double-enzyme digestion of one reference genome. 


    use RestrictionDigest;

    $digest = RestrictionDigest::SingleItem::Double->new();
    
    $digest->add_ref(-reference=>'Full path of the reference file');

    $digest->new_enzyme(-enzyme_name=>'Ncol',-recognition_site=>'C|CATGG');

    $digest->add_enzyme_pair(-front_enzyme=>'EcoRI',-behind_enzyme=>'HinfI');

    $digest->add_output_dir(-output_dir=>'Full path of the output directory');

    $digest->change_range(-start => 301, -end => 500);

    $digest->change_lengths_distribution_parameters(-front => 200, -behind => 800, -step => 25);

    $digest->double_digest();

    $digest->add_SNPs(-SNPs => 'full path to the SNPs coordinate file');

    $digest->count_SNPs_at_fragments(-sequence_type=>'125SE', -sequence_end=>'front_enzyme');

    $digest->add_gff(-gff=>'full path to the GFF file');

    $digest->frags_in_range_coverage();


=cut

use IO::File;
use File::Basename;

use vars qw(%fields %enzyme_list);

# Define the fields to be used as identifiers.
%fields = (
  ref => undef,
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
	  MslI => 'CAYNN|NNRTG',



);

=head1 SUBROUTINES/METHODS

=head2 new
     
     $digest=RestrictionDigest::SingleItem::Double->new(); 
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
          Add the reference file to the object. 
          
=cut

sub add_ref {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-reference$/){
      $self->{ref}=$parameters{$_};
    }
    else{
    die "Unacceptable parameters in the method <add_ref>, please perldoc RestrictionDigest for example or read
    README for more help!\n";
    }
  }
}

=head2 add_enzyme_pair
     
     $digest->add_enzyme_pair(-front_enzyme=>'EcoRI', -behind_enzyme=>'HinfI');
          Add two enzymes to the object. The names of the two enzymes are case-insensitive.

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
      die"Unacceptable parameters in the method <add_enzyme_pair>, please perldoc RestrictionDigest for example
      or read README for more help!\n";
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
    die "The front restrict endonuclease  User provides does not exist.\nHowever User can add this enzyme and 
    ts recognition sites to enzyme-container via the method 'new_enzyme' .\n";
  }
  if($enzyme2_exists >0 ){
    
  }else{
    die "The behind restrict endonuclease  User provides does not exist.\nHowever User can add this enzyme and
    its recognition sites to enzyme-container via the method 'new_enzyme' .\n";
  }
}

=head2 new_enzyme

     $digest->new_enzyme(-enzyme_name=>'EcoRI', -recognition_site=>'G|AATTC');
     
     If the enzyme User wants to use does not exists in the enzyme resevior of the module, the enzyme can be added
     temporarily to the enzyme resevior.

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
        die "The 'recognition_site' must contain a cut flag '|'. Please perldoc RestrictionDigest for example
        or read README for moer help\n";
      }
    }
    else{
      die"Unacceptable parameters in the method <new_enzyme>, please perldoc RestrictionDigest for example or
      read README for more help\n";
    }
  }
  $enzyme_list{$new_enzyme_name}=$new_enzyme_site;
}


=head2 change_range

      $digest->change_range(-start=>201, -end=>500);
       This function is used to change the length range corresponding to size selection process. The default range is from
       201 bp to 500 bp.


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
      die"Unacceptable parameters in the method <change_range>,please perldoc RestrictionDigest for example or
      read README for more help.\n";
    }
  }
  if($self->{range_start} < $self->{range_end} ) {
    
  }else{
    die "The first parameter of the method <change_range> must be smaller than the second parameter\n";
  }
}

=head2 change_lengths_distribution_parameters
    
     $digest->change_lengths_distribution_parameters(-front=>100,-behind=>1000,-step=>50);
      This function is used to change the resolution of fragments lengths distribution. This function has three
      parameters: front and behind define the two boundary length values, and step defines the step length.
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

=cut


sub add_output_dir {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-output_dir$/){
      $self->{output_dir}=$parameters{$_};
    }
    else{
      die "Unacceptable parameter in the method <add_output_dir>, please perldoc RestrictionDigest for
      example or read README for more help!\n";
    }
  }
  if(-d $self->{output_dir} ) {
  }else{
    die "The output directory User provides does not exist.\n";
  }
}

=head2 double_digest

     $digest->double_digest();
      Execute the digestion process. This process will produce several result files which will be 
      ocated in the output directory.

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
  my $all_FB_BF_seq_file_fh=IO::File->new(">>$output_dir/seq_FB_BF_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $all_FF_seq_file_fh=IO::File->new(">>$output_dir/seq_FF_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $all_BB_seq_file_fh=IO::File->new(">>$output_dir/seq_BB_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $all_FB_BF_loc_file_fh=IO::File->new(">>$output_dir/position_FB_BF_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $all_FF_loc_file_fh=IO::File->new(">>$output_dir/position_FF_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $all_BB_loc_file_fh=IO::File->new(">>$output_dir/position_BB_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $seq_FB_BF_in_range_file_fh=IO::File->new(">>$output_dir/seq_FB_BF_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $seq_FF_in_range_file_fh=IO::File->new(">>$output_dir/seq_FF_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $seq_BB_in_range_file_fh=IO::File->new(">>$output_dir/seq_BB_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $loc_FB_BF_in_range_file_fh=IO::File->new(">>$output_dir/position_FB_BF_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $loc_FF_in_range_file_fh=IO::File->new(">>$output_dir/position_FF_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  my $loc_BB_in_range_file_fh=IO::File->new(">>$output_dir/position_BB_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");

  
  my $reduced_ratio_every_chromosome_file_fh=IO::File->new(">>$output_dir/coverage_every_chromosome_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  $reduced_ratio_every_chromosome_file_fh->print("chromosome_name\tall_type_overall_rate\tall_type_in_range_rate\t");
  $reduced_ratio_every_chromosome_file_fh->print("FB_BF_overall_rate\tFB_BF_in_range_rate\tFF_overall_rate\tFF_in_range_rate\tBB_overall_rate\tBB_in_range_rate\n");
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
  my @lengths_of_FB_BF_fragments=();
  my @lengths_of_FF_fragments=();
  my @lengths_of_BB_fragments=();

  my @lengths_fragments_in_range=();
  my @lengths_FB_BF_fragments_in_range=();
  my @lengths_FF_fragments_in_range=();
  my @lengths_BB_fragments_in_range=();

  my $overall_fragments_length=0;
  my $overall_FB_BF_fragments_length=0;
  my $overall_FF_fragments_length=0;
  my $overall_BB_fragments_length=0;

  my $overall_fragments_in_range_length=0;
  my $overall_FB_BF_fragments_in_range_length=0;
  my $overall_FF_fragments_in_range_length=0;
  my $overall_BB_fragments_in_range_length=0;

  my $overall_length_of_scfd=0;

  my $total_front_enzyme_loci=0;
  my $total_later_enzyme_loci=0;
  ########## CIRCLES START HERE!###########

  # Process one scaffold one circle;
  my $seq_tmp="";
  my $seq_name_tmp="";
  my $scaffold_name; my $scaffold_seq;
  while(<$ref_fh>){
    chomp;
    my $line=$_;
    if($_=~/^(>\S+)/){
      my $seq_name=$1;
      my $seq_name_length=length $seq_name_tmp;
      if($seq_name_length !=0){
        $scaffold_name=$seq_name_tmp;
        $scaffold_seq=$seq_tmp;
		$scaffold_seq = uc $scaffold_seq;
        print"Digesting $scaffold_name...\t";
        $reduced_ratio_every_chromosome_file_fh->print("$scaffold_name\t");
        my $length_of_scafd=length $scaffold_seq;
        
        # Set the arrays which contain cut locations of front enzyme and later enzyme;
        my @front_enzyme_locs=();
        
        my @later_enzyme_locs=();
        
        my $seq_for_search= $scaffold_seq;
        # Find locs of recognition of front enzyme in this scaffold;the locs are meaningful in this array; not in the human context;
          my $string=$seq_string_front_enzyme;
		  
		  ## for Enzyme Cut Site with no Degenerate Site.
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
		  ## for Enzyme Cut Site with one Degenerate Site.
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
			 ## for Enzyme Cut Site with two Degenerate Sites.
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
			## for Enzyme Cut Site with two Degenerate Sites.
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
		  ## for Enzyme Cut Site with three Degenerate Sites.
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
		  ## for Enzyme Cut Site with two Degenerate Sites.
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

		  ##  for Enzyme Cut Site with Six Degenerate Sites: MSlI: CAYNNNNRTG.
		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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

		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }




        
        # Count all enzymes loci of the genome.
        $total_front_enzyme_loci+=@front_enzyme_locs;
        $total_later_enzyme_loci+=@later_enzyme_locs;

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

        # Three types of fragments: FB_BF, FF, BB (F for Front_enzyme and B for Behind_enzyme)
        my ($length_FB_BF_fragments,$length_FF_fragments,$length_BB_fragments,$length_FB_BF_fragments_in_range)=qw(0 0 0 0);
        my ($length_FF_fragments_in_range, $length_BB_fragments_in_range)=qw(0 0);
        my ($count_all,$count_in_range,$count_FB_BF,$count_FF,$count_BB,$count_FB_BF_in_range,$count_FF_in_range,$count_BB_in_range)=qw(0 0 0 0 0 0 0 0);
        my ($before_cut_coorr,$after_cut_coorr,$bf_ct_crr_hm_rd,$af_ct_crr_hm_rd);

        my (@FB_BF_locs_before, @FB_BF_enzyme_before, @FB_BF_locs_after, @FB_BF_enzyme_after);
        my (@FF_locs_before,  @FF_locs_after,@BB_locs_before, @BB_locs_after);
        my ($after_loc,$enzyme_before, $enzyme_after);
        my $before_loc=shift @locs_sorted;

        while(@locs_sorted){
          $after_loc=shift @locs_sorted;
          $enzyme_before=$hash{$before_loc};
          $enzyme_after=$hash{$after_loc};
          unless($enzyme_before eq $enzyme_after){
            push @FB_BF_locs_before, $before_loc;
            push @FB_BF_enzyme_before, $enzyme_before;
            push @FB_BF_locs_after,$after_loc;
            push @FB_BF_enzyme_after, $enzyme_after;
            my $length_of_FB_BF_fragment=0; my $seq_FB_BF;
            $count_all++; $count_FB_BF++;
            if($front_enzyme=~/$enzyme_before/i){
              $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
              $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
              $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
              $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);   
            }
            elsif($later_enzyme=~/$enzyme_before/i){
                $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
                $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
                $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
                $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);
            }

            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;

            # Output the sequence and fragment positions of all_length_FB_BF fragments to files;
            $all_FB_BF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
            my $strand;
            if($enzyme_before=~/$front_enzyme/i){$strand="+";}
            elsif($enzyme_before=~/$later_enzyme/i){$strand="-";}
            $all_FB_BF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
            push @lengths_of_fragments, $length_of_FB_BF_fragment;
            push @lengths_of_FB_BF_fragments, $length_of_FB_BF_fragment; 

            # Process fragments in range. 
            if($length_of_FB_BF_fragment >= $range1 && $length_of_FB_BF_fragment <= $range2){
              $count_in_range++;
              $count_FB_BF_in_range++;
              $seq_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
              $loc_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FB_BF_fragment;
              push @lengths_FB_BF_fragments_in_range, $length_of_FB_BF_fragment;
        
              $length_FB_BF_fragments_in_range+=$length_of_FB_BF_fragment;              
            }

            # Sum all fragments length no matter in range or not!!!
            
            $length_FB_BF_fragments+=$length_of_FB_BF_fragment;
          }elsif($enzyme_before eq $front_enzyme){
            push @FF_locs_before, $before_loc;
            push @FF_locs_after, $after_loc;
            my $length_of_FF_fragment=0; my $seq_FF;
            $count_all++; $count_FF++;
            $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
            $length_of_FF_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_FF=substr($scaffold_seq, $before_cut_coorr, $length_of_FF_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            $all_FF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
            my $strand="+";
            $all_FF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
            push @lengths_of_fragments, $length_of_FF_fragment;
            push @lengths_of_FF_fragments, $length_of_FF_fragment;

            # Process fragments in range.
            if($length_of_FF_fragment >= $range1 && $length_of_FF_fragment <= $range2){
              $count_in_range++;
              $count_FF_in_range++;
              $seq_FF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
              $loc_FF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FF_fragment;
              push @lengths_FF_fragments_in_range, $length_of_FF_fragment;
             
              $length_FF_fragments_in_range+=$length_of_FF_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
            
            $length_FF_fragments+=$length_of_FF_fragment;

          }elsif($enzyme_before eq $later_enzyme){
            push @BB_locs_before, $before_loc;
            push @BB_locs_after, $after_loc;
            my $length_of_BB_fragment=0; my $seq_BB;
            $count_all++; $count_BB++;
            $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
            $length_of_BB_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_BB=substr($scaffold_seq, $before_cut_coorr, $length_of_BB_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            $all_BB_seq_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
            my $strand="+";
            $all_BB_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
            push @lengths_of_fragments, $length_of_BB_fragment;
            push @lengths_of_BB_fragments, $length_of_BB_fragment;

            # Process fragments in range.
            if($length_of_BB_fragment >= $range1 && $length_of_BB_fragment <= $range2){
              $count_in_range++;
              $count_BB_in_range++;
              $seq_BB_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
              $loc_BB_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
              push @lengths_fragments_in_range, $length_of_BB_fragment;
              push @lengths_BB_fragments_in_range, $length_of_BB_fragment;
              
              $length_BB_fragments_in_range+=$length_of_BB_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
          
            $length_BB_fragments+=$length_of_BB_fragment;

          }
          $before_loc=$after_loc;
        }

        # Calculate reduced rate of three types of fragments on this chromosome.
        my $all_type_length_overall=$length_FB_BF_fragments+$length_FF_fragments+$length_BB_fragments;
        my $all_type_length_in_range=$length_FB_BF_fragments_in_range+$length_FF_fragments_in_range+$length_BB_fragments_in_range;
        my $all_type_overall_rate=$all_type_length_overall/$length_of_scafd;
        my $all_type_in_range_rate=$all_type_length_in_range/$length_of_scafd;
        my $FB_BF_overall_rate=$length_FB_BF_fragments/$length_of_scafd;
        my $FB_BF_in_range_rate=$length_FB_BF_fragments_in_range/$length_of_scafd;
        my $FF_overall_rate=$length_FF_fragments/$length_of_scafd;
        my $FF_in_range_rate=$length_FF_fragments_in_range/$length_of_scafd;
        my $BB_overall_rate=$length_BB_fragments/$length_of_scafd;
        my $BB_in_range_rate=$length_BB_fragments_in_range/$length_of_scafd;

        $overall_fragments_length+=$all_type_length_overall;
        $overall_fragments_in_range_length+=$all_type_length_in_range;
        $overall_length_of_scfd+=$length_of_scafd;

        $overall_FB_BF_fragments_length+=$length_FB_BF_fragments;
        $overall_FF_fragments_length+=$length_FF_fragments;
        $overall_BB_fragments_length+=$length_BB_fragments;

        $overall_FB_BF_fragments_in_range_length+=$length_FB_BF_fragments_in_range;
        $overall_FF_fragments_in_range_length+=$length_FF_fragments_in_range;
        $overall_BB_fragments_in_range_length+=$length_BB_fragments_in_range;

        $reduced_ratio_every_chromosome_file_fh->print("$all_type_overall_rate\t$all_type_in_range_rate\t$FB_BF_overall_rate\t$FB_BF_in_range_rate\t");
        $reduced_ratio_every_chromosome_file_fh->print("$FF_overall_rate\t$FF_in_range_rate\t$BB_overall_rate\t$BB_in_range_rate\n");
        print"Done!\n";
      }
      $seq_name_tmp=$seq_name;
      $seq_tmp="";
    }else{
      $seq_tmp=$seq_tmp.$line;
    }
    if(eof($ref_fh)){
      $scaffold_name=$seq_name_tmp;
      $scaffold_seq=$seq_tmp;
	  $scaffold_seq = uc $scaffold_seq;
        print"Digesting $scaffold_name...\t";
        $reduced_ratio_every_chromosome_file_fh->print("$scaffold_name\t");
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


		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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

		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }









        # Count all enzymes loci of the genome.
        $total_front_enzyme_loci+=@front_enzyme_locs;
        $total_later_enzyme_loci+=@later_enzyme_locs;
    
     
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
        

        # Three types of fragments: FB_BF, FF, BB (F for Front_enzyme and B for Behind_enzyme)
        my ($length_FB_BF_fragments,$length_FF_fragments,$length_BB_fragments,$length_FB_BF_fragments_in_range)=qw(0 0 0 0);
        my ($length_FF_fragments_in_range, $length_BB_fragments_in_range)=qw(0 0);
        my ($count_all,$count_in_range,$count_FB_BF,$count_FF,$count_BB,$count_FB_BF_in_range,$count_FF_in_range,$count_BB_in_range)=qw(0 0 0 0 0 0 0 0);
        my ($before_cut_coorr,$after_cut_coorr,$bf_ct_crr_hm_rd,$af_ct_crr_hm_rd);

        my (@FB_BF_locs_before, @FB_BF_enzyme_before, @FB_BF_locs_after, @FB_BF_enzyme_after);
        my (@FF_locs_before,  @FF_locs_after,@BB_locs_before, @BB_locs_after);
        my ($after_loc,$enzyme_before, $enzyme_after);
        my $before_loc=shift @locs_sorted;

        while(@locs_sorted){
          $after_loc=shift @locs_sorted;
          $enzyme_before=$hash{$before_loc};
          $enzyme_after=$hash{$after_loc};
          unless($enzyme_before eq $enzyme_after){
            push @FB_BF_locs_before, $before_loc;
            push @FB_BF_enzyme_before, $enzyme_before;
            push @FB_BF_locs_after,$after_loc;
            push @FB_BF_enzyme_after, $enzyme_after;
            my $length_of_FB_BF_fragment=0; my $seq_FB_BF;
            $count_all++; $count_FB_BF++;
            if($front_enzyme=~/$enzyme_before/i){
              $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
              $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
              $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
              $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);   
            }
            elsif($later_enzyme=~/$enzyme_before/i){
                $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
                $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
                $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
                $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);
            }

            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;

            # Output the sequence and fragment positions of all_length_FB_BF fragments to files;
            $all_FB_BF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
            my $strand;
            if($enzyme_before=~/$front_enzyme/i){$strand="+";}
            elsif($enzyme_before=~/$later_enzyme/i){$strand="-";}
            $all_FB_BF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
            push @lengths_of_fragments, $length_of_FB_BF_fragment;
            push @lengths_of_FB_BF_fragments, $length_of_FB_BF_fragment; 

            # Process fragments in range. 
            if($length_of_FB_BF_fragment >= $range1 && $length_of_FB_BF_fragment <= $range2){
              $count_in_range++;
              $count_FB_BF_in_range++;
              $seq_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
              $loc_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FB_BF_fragment;
              push @lengths_FB_BF_fragments_in_range, $length_of_FB_BF_fragment;
              
              $length_FB_BF_fragments_in_range+=$length_of_FB_BF_fragment;              
            }

            # Sum all fragments length no matter in range or not!!!
            
            $length_FB_BF_fragments+=$length_of_FB_BF_fragment;
          }elsif($enzyme_before eq $front_enzyme){
            push @FF_locs_before, $before_loc;
            push @FF_locs_after, $after_loc;
            my $length_of_FF_fragment=0; my $seq_FF;
            $count_all++; $count_FF++;
            $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
            $length_of_FF_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_FF=substr($scaffold_seq, $before_cut_coorr, $length_of_FF_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            $all_FF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
            my $strand="+";
            $all_FF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
            push @lengths_of_fragments, $length_of_FF_fragment;
            push @lengths_of_FF_fragments, $length_of_FF_fragment;

            # Process fragments in range.
            if($length_of_FF_fragment >= $range1 && $length_of_FF_fragment <= $range2){
              $count_in_range++;
              $count_FF_in_range++;
              $seq_FF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
              $loc_FF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FF_fragment;
              push @lengths_FF_fragments_in_range, $length_of_FF_fragment;
           
              $length_FF_fragments_in_range+=$length_of_FF_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
            
            $length_FF_fragments+=$length_of_FF_fragment;

          }elsif($enzyme_before eq $later_enzyme){
            push @BB_locs_before, $before_loc;
            push @BB_locs_after, $after_loc;
            my $length_of_BB_fragment=0; my $seq_BB;
            $count_all++; $count_BB++;
            $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
            $length_of_BB_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_BB=substr($scaffold_seq, $before_cut_coorr, $length_of_BB_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            $all_BB_seq_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
            my $strand="+";
            $all_BB_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
            push @lengths_of_fragments, $length_of_BB_fragment;
            push @lengths_of_BB_fragments, $length_of_BB_fragment;

            # Process fragments in range.
            if($length_of_BB_fragment >= $range1 && $length_of_BB_fragment <= $range2){
              $count_in_range++;
              $count_BB_in_range++;
              $seq_BB_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
              $loc_BB_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
              push @lengths_fragments_in_range, $length_of_BB_fragment;
              push @lengths_BB_fragments_in_range, $length_of_BB_fragment;
             
              $length_BB_fragments_in_range+=$length_of_BB_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
           
            $length_BB_fragments+=$length_of_BB_fragment;

          }
          $before_loc=$after_loc;
        }

        # Calculate reduced rate of three types of fragments on this chromosome.
        my $all_type_length_overall=$length_FB_BF_fragments+$length_FF_fragments+$length_BB_fragments;
        my $all_type_length_in_range=$length_FB_BF_fragments_in_range+$length_FF_fragments_in_range+$length_BB_fragments_in_range;
        my $all_type_overall_rate=$all_type_length_overall/$length_of_scafd;
        my $all_type_in_range_rate=$all_type_length_in_range/$length_of_scafd;
        my $FB_BF_overall_rate=$length_FB_BF_fragments/$length_of_scafd;
        my $FB_BF_in_range_rate=$length_FB_BF_fragments_in_range/$length_of_scafd;
        my $FF_overall_rate=$length_FF_fragments/$length_of_scafd;
        my $FF_in_range_rate=$length_FF_fragments_in_range/$length_of_scafd;
        my $BB_overall_rate=$length_BB_fragments/$length_of_scafd;
        my $BB_in_range_rate=$length_BB_fragments_in_range/$length_of_scafd;

        $overall_fragments_length+=$all_type_length_overall;
        $overall_fragments_in_range_length+=$all_type_length_in_range;
        $overall_length_of_scfd+=$length_of_scafd;

        $overall_FB_BF_fragments_length+=$length_FB_BF_fragments;
        $overall_FF_fragments_length+=$length_FF_fragments;
        $overall_BB_fragments_length+=$length_BB_fragments;

        $overall_FB_BF_fragments_in_range_length+=$length_FB_BF_fragments_in_range;
        $overall_FF_fragments_in_range_length+=$length_FF_fragments_in_range;
        $overall_BB_fragments_in_range_length+=$length_BB_fragments_in_range;

        $reduced_ratio_every_chromosome_file_fh->print("$all_type_overall_rate\t$all_type_in_range_rate\t$FB_BF_overall_rate\t$FB_BF_in_range_rate\t");
        $reduced_ratio_every_chromosome_file_fh->print("$FF_overall_rate\t$FF_in_range_rate\t$BB_overall_rate\t$BB_in_range_rate\n");
        print"Done!\n";
    }
  }

  #Summary the whole reference genome.
  my $overall_all_type_coverage=$overall_fragments_length/$overall_length_of_scfd;
  my $all_type_in_range_coverage=$overall_fragments_in_range_length/$overall_length_of_scfd;
  my $overall_FB_BF_coverage=$overall_FB_BF_fragments_length/$overall_length_of_scfd;
  my $FB_BF_in_range_coverage=$overall_FB_BF_fragments_in_range_length/$overall_length_of_scfd;
  my $overall_FF_coverage=$overall_FF_fragments_length/$overall_length_of_scfd;
  my $FF_in_range_coverage=$overall_FF_fragments_in_range_length/$overall_length_of_scfd;
  my $overall_BB_coverage=$overall_BB_fragments_length/$overall_length_of_scfd;
  my $BB_in_range_coverage=$overall_BB_fragments_in_range_length/$overall_length_of_scfd;

  my $count_all=@lengths_of_fragments;
  my $count_FB_BF=@lengths_of_FB_BF_fragments;
  my $count_FF=@lengths_of_FF_fragments;
  my $count_BB=@lengths_of_BB_fragments;
  my $count_in_range=@lengths_fragments_in_range;
  my $count_FB_BF_in_range=@lengths_FB_BF_fragments_in_range;
  my $count_FF_in_range=@lengths_FF_fragments_in_range;
  my $count_BB_in_range=@lengths_BB_fragments_in_range;

  $summary_digestion_fh->print("loci_number_of_$front_enzyme\tloci_number_of_$later_enzyme\n");
  $summary_digestion_fh->print("$total_front_enzyme_loci\t$total_later_enzyme_loci\n");
  # Output the summary of digestion.
  $summary_digestion_fh->print("all_3types_fragments_number\tall_3types_fragments_coverage\tall_FB_BF_fragments_number\t");
  $summary_digestion_fh->print("all_FB_BF_fragments_coverage\tall_FF_fragments_number\tall_FF_fragments_coverage\t");
  $summary_digestion_fh->print("all_BB_fragments_number\tall_BB_fragments_coverage\n");

  $summary_digestion_fh->print("$count_all\t$overall_all_type_coverage\t$count_FB_BF\t$overall_FB_BF_coverage\t");
  $summary_digestion_fh->print("$count_FF\t$overall_FF_coverage\t$count_BB\t$overall_BB_coverage\n");


  $summary_digestion_fh->print("3types_fragments_in_range_number\t3types_fragments_in_range_coverage\t");
  $summary_digestion_fh->print("FB_BF_fragments_in_range_number\tFB_BF_fragments_in_range_coverage\t");
  $summary_digestion_fh->print("FF_fragments_in_range_number\tFF_fragments_in_range_coverage\t");
  $summary_digestion_fh->print("BB_fragments_in_range_number\tBB_fragments_in_range_coverage\n");

  $summary_digestion_fh->print("$count_in_range\t$all_type_in_range_coverage\t$count_FB_BF_in_range\t$FB_BF_in_range_coverage");
  $summary_digestion_fh->print("$count_FF_in_range\t$FF_in_range_coverage\t$count_BB_in_range\t$BB_in_range_coverage\n");


  my $lengths_distribution_start=$self->{lengths_distribution_start};
  my $lengths_distribution_end=$self->{lengths_distribution_end};
  my $lengths_distribution_step=$self->{lengths_distribution_step};
  my $lengths_distribution_tmp;
  my (@starts,@ends,@counts,%three_types,%FB_BF,%FF,%BB);

  my $num_pair=0;

  for($lengths_distribution_tmp= $lengths_distribution_start;$lengths_distribution_tmp<$lengths_distribution_end+$lengths_distribution_step;$lengths_distribution_tmp=$lengths_distribution_tmp+$lengths_distribution_step){
    if($lengths_distribution_tmp== $lengths_distribution_start){
      $num_pair++;
      push @starts, '0';
      push @ends , $lengths_distribution_tmp;
      $three_types{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
      $FB_BF{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
      $FF{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
      $BB{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;

      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      $three_types{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      $FB_BF{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      $FF{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      $BB{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      
    }
    elsif($lengths_distribution_tmp  < $lengths_distribution_end ){
      if($lengths_distribution_tmp+$lengths_distribution_step <$lengths_distribution_end ){
        $num_pair++;
        my $tmp_start=$lengths_distribution_tmp+1;
        my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
        push @starts, $tmp_start;
        push @ends, $tmp_end;
        $three_types{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        $FB_BF{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        $FF{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        $BB{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        
    }
    elsif($lengths_distribution_tmp +$lengths_distribution_step >= $lengths_distribution_end){
      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_end;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      $three_types{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      $FB_BF{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      $FF{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      $BB{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;

      $num_pair++;
      $tmp_start=$tmp_end+1;
      
      push @starts, $tmp_start;
      push @ends, 'bigger';
      $three_types{$num_pair}{"${tmp_start}bp-longer"}=0;
      $FB_BF{$num_pair}{"${tmp_start}bp-longer"}=0;
      $FF{$num_pair}{"${tmp_start}bp-longer"}=0;
      $BB{$num_pair}{"${tmp_start}bp-longer"}=0;
      
      }
    }
  }
  for (@lengths_of_FB_BF_fragments){
    my $length=$_;
    my $num_pair=1;
    my $num_last_pair=@starts;
    for($num_pair=1;$num_pair<=$num_last_pair;$num_pair++){
      
      my $array_index=$num_pair-1;
      my $start=$starts[$array_index];
      my $end=$ends[$array_index];
      if($num_pair<$num_last_pair){
        if($length>=$start && $length<=$end){
          for my $description (keys %{$FB_BF{$num_pair}}){
            $FB_BF{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
      elsif($num_pair==$num_last_pair){
        if($length>=$start){
          for my $description (keys %{$FB_BF{$num_pair}}){
            $FB_BF{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
    }
  }
  for (@lengths_of_FF_fragments){
    my $length=$_;
    my $num_pair=1;
    my $num_last_pair=@starts;
    for($num_pair=1;$num_pair<=$num_last_pair;$num_pair++){
      
      my $array_index=$num_pair-1;
      my $start=$starts[$array_index];
      my $end=$ends[$array_index];
      if($num_pair<$num_last_pair){
        if($length>=$start && $length<=$end){
          for my $description (keys %{$FF{$num_pair}}){
            $FF{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
      elsif($num_pair==$num_last_pair){
        if($length>=$start){
          for my $description (keys %{$FF{$num_pair}}){
            $FF{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
    }
  }
  for (@lengths_of_BB_fragments){
    my $length=$_;
    my $num_pair=1;
    my $num_last_pair=@starts;
    for($num_pair=1;$num_pair<=$num_last_pair;$num_pair++){
      
      my $array_index=$num_pair-1;
      my $start=$starts[$array_index];
      my $end=$ends[$array_index];
      if($num_pair<$num_last_pair){
        if($length>=$start && $length<=$end){
          for my $description (keys %{$BB{$num_pair}}){
            $BB{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
      elsif($num_pair==$num_last_pair){
        if($length>=$start){
          for my $description (keys %{$BB{$num_pair}}){
            $BB{$num_pair}{$description}++;
            last;
          }
          last;
        }
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
          for my $description (keys %{$three_types{$num_pair}}){
            $three_types{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
      elsif($num_pair==$num_last_pair){
        if($length>=$start){
          for my $description (keys %{$three_types{$num_pair}}){
            $three_types{$num_pair}{$description}++;
            last;
          }
          last;
        }
      }
    }
  }
  
  $summary_digestion_fh->print("\nLength_range\tthree_types_number\tthree_types_ratio\tFB_BF_number\tFB_BF_ratio\tFF_number\tFF_ratio\tBB_number\tBB_ratio\n");



  
  for(1..@starts){
    my $num_pair=$_;

    for my $description(keys %{$three_types{$num_pair}}){

      my $three_types_fragments_num=$three_types{$num_pair}{$description};
      my $FB_BF_fragments_num=$FB_BF{$num_pair}{$description};
      my $FF_fragments_num=$FF{$num_pair}{$description};
      my $BB_fragments_num=$BB{$num_pair}{$description};

      my ($three_types_ratio,$FB_BF_ratio,$FF_ratio,$BB_ratio);

      if($count_all){
        $three_types_ratio=$three_types_fragments_num/$count_all;
      }else{
        $three_types_ratio="N/A";
      }

      if($count_FB_BF){
        $FB_BF_ratio=$FB_BF_fragments_num/$count_FB_BF;
      }else{
        $FB_BF_ratio="N/A";
      } 

      if($count_FF){
        $FF_ratio=$FF_fragments_num/$count_FF;
      }else{
        $FF_ratio="N/A";
      }

      if($count_BB){
        $BB_ratio=$BB_fragments_num/$count_BB;
      }else{
        $BB_ratio="N/A";
      }

      $summary_digestion_fh->print("$description\t$three_types_fragments_num\t$three_types_ratio\t");
      $summary_digestion_fh->print("$FB_BF_fragments_num\t$FB_BF_ratio\t");
      $summary_digestion_fh->print("$FF_fragments_num\t$FF_ratio\t");
      $summary_digestion_fh->print("$BB_fragments_num\t$BB_ratio\n");
    }
  }
  $summary_digestion_fh->print("\n");
}

=head2 add_SNPs
  
    $digest->add_SNPs(-SNPs=>'Full path of the SNPs file');
      Adds the SNPs file to the object.
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
        Count the expected SNPs located within the framgents.
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
        die "Unacceptable parameters in the method <count_SNPs_at_fragments>, please perldoc RestrictionDigest
        for example or read README for more help\n";
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
      die"Unacceptable sequence_type when range's front border is $range1 bp. Please perldoc RestrictionDigest
      for example or read README for more help!\n";
    }
    else{
      $pe_sequence_length=$sequence_length;
    }
  }
  elsif($sequence_type=~/^(\d+)se$/i){
    $sequence_type_code=1;
    my $sequence_length=$1;
    unless($sequence_length<=$range1){
      die "Unacceptable sequence_type when range's front border is $range1 bp. Please perldoc RestrictionDigest
      for example or read README for more help!\n";
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

  my $FB_BF_loc_fh=IO::File->new("$output_dir/position_FB_BF_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  my $FB_BF_loc_in_range_fh=IO::File->new("$output_dir/position_FB_BF_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  my $FF_loc_fh=IO::File->new("$output_dir/position_FF_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  my $FF_loc_in_range_fh=IO::File->new("$output_dir/position_FF_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  my $BB_loc_fh=IO::File->new("$output_dir/position_BB_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  my $BB_loc_in_range_fh=IO::File->new("$output_dir/position_BB_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}",'r');
  
  # Hold SNPs into hash.
  my $SNPs_fh=IO::File->new("$SNPs",'r');
  my %SNPs;
  while(<$SNPs_fh>){
    chomp;
    my ($SNPs_scaffold_name,$SNPs_pos,$type)=split /\s+/, $_;
    $SNPs{$SNPs_scaffold_name}{$SNPs_pos}=1;
    #push @{$SNPs{$SNPs_scaffold_name}}, $SNPs_pos;
  }
  print "Counting SNPs at three types of fragments...\t";
  
  # Create the SNPs files to hold SNPs of all fragments and fragments in range.
  #my $SNPs_at_FB_BF_all_frags_fh=IO::File->new(">>$output_dir/SNPs_at_FB_BF_all_frags_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");
  #my $SNPs_at_FB_BF_frags_in_range_fh=IO::File->new(">>$output_dir/SNPs_at_FB_BF_frags_in_range_${name_reference}_by_${front_enzyme}_and_${later_enzyme}");

  # Scan all locs file to count SNPs and output the result to the summary file.
  my ($FB_BF_SNPs_number, $FF_SNPs_number, $BB_SNPs_number)=qw(0 0 0);

  if($sequence_type_code == 1){
    if($se_sequence_end_code==1){
      while(<$FB_BF_loc_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)$/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          if($front_enzyme=~/$enzyme1/i){
            my $frag_start=$pos1;
            my $frag_end=$pos1+$se_sequence_length-1;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_number++;
              }
            }
          }
          elsif($front_enzyme=~/$enzyme2/i){
            my $frag_start=$pos2-$se_sequence_length+1;
            my $frag_end=$pos2;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_number++;
              }
            }
          }
        }
      }
    }
    elsif($se_sequence_end_code==2){
      while(<$FB_BF_loc_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)$/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          if($later_enzyme=~/$enzyme1/i){
            my $frag_start=$pos1;
            my $frag_end=$pos1+$se_sequence_length-1;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_number++;
              }
            }
          }
          elsif($later_enzyme=~/$enzyme2/i){
            my $frag_start=$pos2-$se_sequence_length+1;
            my $frag_end=$pos2;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_number++;
              }
            }
          }
        }
      }
    }
    while(<$FF_loc_fh>){
      chomp;
      my ($fragment_name,$strand, $enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$se_sequence_length-1;
        my $frag_start2=$pos2-$se_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_number++;
          }
        }
      }
    }
    $FF_SNPs_number=int($FF_SNPs_number/2);
    while(<$BB_loc_fh>){
        chomp;
        my ($fragment_name,$strand, $enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)$/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$se_sequence_length-1;
        my $frag_start2=$pos2-$se_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_number++;
          }
        }
      }
    }
    $BB_SNPs_number=int($BB_SNPs_number/2);
    $summary_digestion_fh->print("Expected_SNPs_at_FB_BF_all_fragments_via_$sequence_type:\t$FB_BF_SNPs_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_FF_all_fragments_via_$sequence_type:\t$FF_SNPs_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_BB_all_fragments_via_$sequence_type:\t$BB_SNPs_number\n");
  }
  elsif($sequence_type_code ==2){
    while(<$FB_BF_loc_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FB_BF_SNPs_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FB_BF_SNPs_number++;
          }
        }
      }
    }
    while(<$FF_loc_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_number++;
          }
        }
      }
    }
    while(<$BB_loc_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_number++;
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected_SNPs_at_FB_BF_all_fragments_via_$sequence_type:\t$FB_BF_SNPs_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_FF_all_fragments_via_$sequence_type:\t$FF_SNPs_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_BB_all_fragments_via_$sequence_type:\t$BB_SNPs_number\n");
  }

  
  # Scan locs in range file to count SNPs and output the result to the summary file.
  my ($FB_BF_SNPs_in_range_number, $FF_SNPs_in_range_number, $BB_SNPs_in_range_number)=qw(0 0 0);

  # Single end sequencing.
  if($sequence_type_code == 1){
    # Process FB_BF file.
    if($se_sequence_end_code==1){
      while(<$FB_BF_loc_in_range_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)$/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          if($front_enzyme=~/$enzyme1/i){
            my $frag_start=$pos1;
            my $frag_end=$pos1+$se_sequence_length-1;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_in_range_number++;
              }
            }
          }
          elsif($front_enzyme=~/$enzyme2/i){
            my $frag_start=$pos2-$se_sequence_length+1;
            my $frag_end=$pos2;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_in_range_number++;
              }
            }
          }
        }
      }
    }
    elsif($se_sequence_end_code==2){
      while(<$FB_BF_loc_in_range_fh>){
        chomp;
        my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
        my $scaffold_name;
        if($fragment_name=~/^>(.+)-(\d+)$/){
          $scaffold_name=$1;
        }
        if($SNPs{$scaffold_name}){
          if($later_enzyme=~/$enzyme1/i){
            my $frag_start=$pos1;
            my $frag_end=$pos1+$se_sequence_length-1;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_in_range_number++;
              }
            }
          }
          elsif($later_enzyme=~/$enzyme2/i){
            my $frag_start=$pos2-$se_sequence_length+1;
            my $frag_end=$pos2;
            for my $frag_pos ($frag_start..$frag_end){
              if($SNPs{$scaffold_name}{$frag_pos}){
                $FB_BF_SNPs_in_range_number++;
              }
            }
          }
        }
      }
    }
    # Process FF file.
    while(<$FF_loc_in_range_fh>){
      chomp;
      my ($fragment_name,$strand, $enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$se_sequence_length-1;
        my $frag_start2=$pos2-$se_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_in_range_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_in_range_number++;
          }
        }
      }
    }
    $FF_SNPs_in_range_number=int($FF_SNPs_in_range_number/2);
    # Process BB file.
    while(<$BB_loc_in_range_fh>){
      chomp;
      my ($fragment_name,$strand, $enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$se_sequence_length-1;
        my $frag_start2=$pos2-$se_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_in_range_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_in_range_number++;
          }
        }
      }
    }
    $BB_SNPs_in_range_number=int($BB_SNPs_in_range_number/2);
    # Output SNPs numbers.
    $summary_digestion_fh->print("Expected_SNPs_at_FB_BF_fragments_in_range_via_$sequence_type:\t$FB_BF_SNPs_in_range_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_FF_fragments_in_range_via_$sequence_type:\t$FF_SNPs_in_range_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_BB_fragments_in_range_via_$sequence_type:\t$BB_SNPs_in_range_number\n");
  }
  # Pair end sequencing.
  elsif($sequence_type_code ==2){
    while(<$FB_BF_loc_in_range_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FB_BF_SNPs_in_range_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FB_BF_SNPs_in_range_number++;
          }
        }
      }
    }
    while(<$FF_loc_in_range_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_in_range_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $FF_SNPs_in_range_number++;
          }
        }
      }
    }
    while(<$BB_loc_in_range_fh>){
      chomp;
      my($fragment_name,$strand,$enzyme1,$pos1,$enzyme2,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_in_range_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $BB_SNPs_in_range_number++;
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected_SNPs_at_FB_BF_fragments_in_range_via_$sequence_type:\t$FB_BF_SNPs_in_range_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_FF_fragments_in_range_via_$sequence_type:\t$FF_SNPs_in_range_number\n");
    $summary_digestion_fh->print("Expected_SNPs_at_BB_fragments_in_range_via_$sequence_type:\t$BB_SNPs_in_range_number\n");
  }
  print "Done!\n";
}


=head2 add_gff

    $digest->add_gff(-ref=>'Full path of the gff file');
        Adds the Gff file to RestrictionDigest the object. 

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

=head2 all_frags_coverage

     $digest->all_frags_coverage();
          Calculate the coverage of all restriction fragments on different genome regions.

=cut

sub all_frags_coverage {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $e1=$self->{enzyme1};
  my $e2=$self->{enzyme2};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $FB_BF_loc_all_frags="$output_dir/position_FB_BF_frags_${basename_ref}_by_${e1}_and_${e2}";
  my $FF_loc_all_frags="$output_dir/position_FF_frags_${basename_ref}_by_${e1}_and_${e2}";
  my $BB_loc_all_frags="$output_dir/position_BB_frags_${basename_ref}_by_${e1}_and_${e2}";
  $self->genome_structure_coverage_ratio($FB_BF_loc_all_frags);
  $self->genome_structure_coverage_ratio($FF_loc_all_frags);
  $self->genome_structure_coverage_ratio($BB_loc_all_frags);
}

=head2 frags_in_range_coverage

     $digest->frags_in_range_coverage();
      Calculate the coverage of restriction fragments on different genome regions.
 
=cut

sub frags_in_range_coverage {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $e1=$self->{enzyme1};
  my $e2=$self->{enzyme2};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $FB_BF_loc_frags_in_range="$output_dir/position_FB_BF_frags_in_range_${basename_ref}_by_${e1}_and_${e2}";
  my $FF_loc_frags_in_range="$output_dir/position_FF_frags_in_range_${basename_ref}_by_${e1}_and_${e2}";
  my $BB_loc_frags_in_range="$output_dir/position_BB_frags_in_range_${basename_ref}_by_${e1}_and_${e2}";

  $self->genome_structure_coverage($FB_BF_loc_frags_in_range);
  $self->genome_structure_coverage($FF_loc_frags_in_range);
  $self->genome_structure_coverage($BB_loc_frags_in_range);
}

=head2  genome_structure_coverage

        A subroutine invoked by the 'all_frags_coverage' and the 'frags_in_range_coverage' methods. 
        User do not use this subroutine.

=cut

 
sub genome_structure_coverage {
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
  
  # Make  output filehandles to generic region and intergenic region coverage  files. 
  my $coverage_ratio_fh;
  
  my $basename_loc_frags_file=basename($loc_frags_file);

  if($basename_loc_frags_file =~/frags_in_range/){
    if($basename_loc_frags_file=~/^position_FB_BF_frags_in_range/){
      $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_of_FB_BF_frags_in_range_${basename_ref}_by_${e1}_and_${e2}");
    }
    elsif($basename_loc_frags_file=~/^position_FF_frags_in_range/){
      $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_of_FF_frags_in_range_${basename_ref}_by_${e1}_and_${e2}");
    }
    elsif($basename_loc_frags_file=~/^position_BB_frags_in_range/){
      $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_of_BB_frags_in_range_${basename_ref}_by_${e1}_and_${e2}");
    }
  }
  else{
    if($basename_loc_frags_file=~/^position_FB_BF_frags/){
      $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_of_FB_BF_frags_${basename_ref}_by_${e1}_and_${e2}");
    }
    elsif($basename_loc_frags_file=~/^position_FF_frags/){
      $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_of_FF_frags_${basename_ref}_by_${e1}_and_${e2}");
    }
    elsif($basename_loc_frags_file=~/^position_BB_frags/){
      $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_of_BB_frags_${basename_ref}_by_${e1}_and_${e2}");
    }
  }

  $coverage_ratio_fh->print("ScaffoldName\tIntergenicLength\tIntergenicMapLength\tIntergenicMapRatio\tGenesLength\tGenesMapLength\tGenesMapRatio\tExonsLength\tExonsMapLength\tExonsMapRatio\tIntronsLength\tIntronsMapLength\tIntronsMapRatio\n");


  # Set variables used in this program.
  my($all_intergenic_length,$all_intergenic_map_length,$all_gene_length,$all_gene_map_length,$all_exon_length,$all_exon_map_length,$all_intron_length,$all_intron_map_length)=qw(0 0 0 0 0 0 0 0);
  


  # Fetch the lengths of different parts of this scaffold from the GFF file.
      my $gff_fh=IO::File->new($gff,'r');
      # Create two hashes to check if start or stop position exists.
      my %gff; my %exon_check; my %gene_check;
      print "Reading the GFF file...\t";

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

      # Read the loc file and fetch the information.
      
      my $loc_frags_file_fh=IO::File->new($loc_frags_file,'r');
      my %all_locs;
      print "Reading the positions file...\t";
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


  # Process every scaffold one by one.
  my $seq_tmp="";
  my $seq_name_tmp="";
  my $scfd_name;
  my $scfd_seq;

  while(<$ref_fh>){
    chomp;
    my $line=$_;
    if($_=~/^>(\S+)/){
      my $seq_name=$1;
      my $seq_name_length=length $seq_name_tmp;
      if($seq_name_length !=0){
        $scfd_name=$seq_name_tmp;
        $scfd_seq=$seq_tmp;
      
        my $scaffold_length=length $scfd_seq;
        print "Coverage: Processing $scfd_name...\t";


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
      $seq_name_tmp=$seq_name;
      $seq_tmp="";
    }else{
      $seq_tmp=$seq_tmp.$line;
    }
    if(eof($ref_fh)){
      $scfd_name=$seq_name_tmp;
      $scfd_seq=$seq_tmp;
      my $scaffold_length=length $scfd_seq;
      print "Coverage: Processing $scfd_name...\t";


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
  
  my $all_map_length=$all_gene_map_length+$all_intergenic_map_length;
  my $frag_length_on_intergenic_rate=$all_intergenic_map_length/$all_map_length;
  my $frag_length_on_gene_rate=$all_gene_map_length/$all_map_length;
  my $frag_length_on_exon_rate=$all_exon_map_length/$all_map_length;
  my $frag_length_on_intron_rate=$all_intron_map_length/$all_map_length;

  $coverage_ratio_fh->print("Intergenic region coverage is\t$all_intergenic_map_ratio\nGenes region coverage is\t$all_gene_map_ratio\nExon region coverage  is\t$all_exon_map_ratio\nIntron region coverage is\t$all_intron_map_ratio\n");
  $coverage_ratio_fh->print("Total length of fragments:\t$all_map_length\n");
  $coverage_ratio_fh->print("Total bases mapped on intergenic regions:\t$all_intergenic_map_length\tLength rate on all fragments bases:\t$frag_length_on_intergenic_rate\n");
  $coverage_ratio_fh->print("Total bases mapped on gene regions:\t$all_gene_map_length\tLength rate on all fragments bases:\t$frag_length_on_gene_rate\n");
  $coverage_ratio_fh->print("Total bases mapped on exon regions:\t$all_exon_map_length\tLength rate on all fragments bases:\t$frag_length_on_exon_rate\n");
  $coverage_ratio_fh->print("Total bases mapped on intron regions:\t$all_intron_map_length\tLength rate on all fragments bases:\t$frag_length_on_intron_rate\n");
}




package RestrictionDigest::SingleItem::Single;

use 5.8.0;
use strict;
use warnings FATAL => 'all';

=head1 NAME


RestrictionDigest::SingleItem::Single

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';



=head1 SYNOPSIS

RestrictionDigest::SingleItem::Single is used for the simulation of single-enzyme digestion of one reference genome.


    use RestrictionDigest;

    $digest=RestrictionDigest::SingleItem::Single->new();
    
    $digest->add_ref(-reference=>'Full path of the reference file');

    $digest->new_enzyme(-enzyme_name=>'Ncol',-recognition_site=>'C|CATGG');

    $digest->add_single_enzyme(-enzyme=>'EcoRI');

    $digest->add_output_dir(-output_dir=>'Full path of the output directory');

    $digest->change_range(-start => 301, -end => 500);

    $digest->change_lengths_distribution_parameters(-front => 200, -behind => 800, -step => 25);

    $digest->single_digest();

    $digest->add_SNPs(-SNPs => 'full path to the SNPs coordinate file');

    $digest->count_SNPs_at_fragments(-sequence_type=>'125SE', -sequence_end=>'front_enzyme');

<<<<<<< HEAD
    $digest->add_gff(-gff=>'full path to the gff file');
=======
    $digest->add_gff(-gff=>'full path to the GFF file');
>>>>>>> origin/master

    $digest->frags_in_range_coverage();


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
	  MslI => 'CAYNN|NNRTG',
);

=head1 SUBROUTINES/METHODS

=head2 new
     
     $digest=RestrictionDigest::SingleItem::Single->new(); 
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
          Add the reference file to the object.

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
}

=head2 add_single_enzyme
     
     $digest->add_single_enzyme(-enzyme=>'EcoRI');
          Add the single enzyme name used to digest the reference genome.
          
=cut


sub add_single_enzyme {
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-enzyme$/){
      $self->{enzyme}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters in the method <add_singel_enzyme>, please perldoc RestrictionDigest
      for example or read README for more help!\n";
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
      If the enzyme User wants to use does not exists in the enzyme resevior of the module, the enzyme can be added
     temporarily to the enzyme resevior.

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

     $digest->change_lengths_distribution_parameters(-front=>100,-behind=>1000,-step=>50);
      This function is used to change the resolution of fragments lengths distribution. This function has three
      parameters: front and behind define the two boundary length values, and step defines the step length.

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
      This function is used to change the resolution of fragments lengths distribution. This function has three
      parameters: front and behind define the two boundary length values, and step defines the step length.
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
      die"Unacceptable parameters in the method <change_lengths_distribution_parameters>, please
      perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
}



=head2 add_output_dir
    
     $digest->add_output_dir(-output_dir=>'Full path of the output directory');
          Adds the directory into which the result files are put by the program. 

=cut

sub add_output_dir {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-output_dir$/){
      $self->{output_dir}=$parameters{$_};
    }
    else{
      die "Unacceptable parameter in the method <add_output_dir>, please perldoc RestrictionDigest
      for example or read README for more help!\n";
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
          This process will produce several result files which will be located in the output directory a
          dded through the add_output_dir method. 

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
  my $reduced_ratio_every_scaffold_file_fh=IO::File->new(">>$output_dir/coverage_every_chromosome_${name_reference}_by_${enzyme}");
  my $summary_digestion_fh=IO::File->new(">>$output_dir/digestion_summary_${name_reference}_by_${enzyme}");
  $reduced_ratio_every_scaffold_file_fh->print("chromosome_name\toverall_rate\tin_range_rate\n");

  
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
  my $total_enzyme_loci=0;
  ########## CIRCLES START HERE!###########

  # Process one scaffold one circle;
  my $seq_tmp="";
  my $seq_name_tmp="";
  my ($scaffold_name, $scaffold_seq);
  while(<$ref_fh>){
    chomp;
    my $line=$_;
    if($_=~/^>/){
      my $seq_name=$line;
      my $seq_name_length=length $seq_name_tmp;
      if($seq_name_length !=0){
        $scaffold_name=$seq_name_tmp;
        $scaffold_seq=$seq_tmp;
		$scaffold_seq = uc $scaffold_seq;
        
        print"Digesting $scaffold_name...\t";
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
        
        elsif($enzyme_seq=~/CCNNNNNNN/i){
			  my $loc_in_array;
			  $loc_in_array=index($seq_for_search,'CC');
			  unless($loc_in_array==-1){
				  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				  if(length $seq_after_CC >19){
					push @enzyme_locs, $loc_in_array;
				  }
			  }
			  while($loc_in_array!=-1){
				  $loc_in_array=index($seq_for_search, 'CC', $loc_in_array+1);
				  unless($loc_in_array==-1){
					  my $previous_loc=$enzyme_locs[$#enzyme_locs];
					  if(($loc_in_array- $previous_loc) > 13){
						  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				          if(length $seq_after_CC >19){
							  push @enzyme_locs, $loc_in_array;
						  }
					  }
				  }
			  }
		 }
		 elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/i){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }





        
        
        # Sort the array.
        $total_enzyme_loci+=@enzyme_locs;
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
      $seq_name_tmp=$seq_name;
      $seq_tmp="";
    }
    else{
      $seq_tmp=$seq_tmp.$line;
    }
    if(eof($ref_fh)){
      $scaffold_name=$seq_name_tmp;
      $scaffold_seq=$seq_tmp;
	  $scaffold_seq = uc $scaffold_seq;
      
        print"Digesting $scaffold_name...\t";
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
        elsif($enzyme_seq=~/CCNNNNNNN/i){
			  my $loc_in_array;
			  $loc_in_array=index($seq_for_search,'CC');
			  unless($loc_in_array==-1){
				  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				  if(length $seq_after_CC >19){
					push @enzyme_locs, $loc_in_array;
				  }
			  }
			  while($loc_in_array!=-1){
				  $loc_in_array=index($seq_for_search, 'CC', $loc_in_array+1);
				  unless($loc_in_array==-1){
					  my $previous_loc=$enzyme_locs[$#enzyme_locs];
					  if(($loc_in_array- $previous_loc) > 13){
						  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				          if(length $seq_after_CC >19){
							  push @enzyme_locs, $loc_in_array;
						  }
					  }
				  }
			  }
		 }

		 elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/i){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }



        
        
        # Sort the array.
        $total_enzyme_loci+=@enzyme_locs;
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

  $summary_digestion_fh->print("loci_number_of_$enzyme\n$total_enzyme_loci\n");
  # Output the summary of digestion.
  $summary_digestion_fh->print("all_fragments_number\tall_fragments_coverage\t");
  $summary_digestion_fh->print("fragments_in_range_number\tfragments_in_range_coverage\n");
  $summary_digestion_fh->print("$num_all_frags\t$length_ratio_all_frags\t");
  $summary_digestion_fh->print("$num_frags_in_range\t$length_ratio_frags_in_range\n");
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
      $lengths_distribution{$num_pair}{"${tmp_start}bp-longer"}=0;
      
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
  
  $summary_digestion_fh->print("\nLength_range\tfragments_number\tfragments_ratio\n");
  
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
        Adds the SNPs coordinate file to the object.
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
        Count the expected SNPs located within the framgents.
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
        die "Unacceptable parameters in the method <count_SNPs_at_fragments>, please perldoc
        RestrictionDigest for example or read README for more help\n";
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
  
  # Store SNPs into hash.
  my $SNPs_fh=IO::File->new("$SNPs",'r');
  my %SNPs;
  while(<$SNPs_fh>){
    chomp;
    my ($SNPs_scaffold_name,$SNPs_pos,$type)=split /\s+/, $_;
    $SNPs{$SNPs_scaffold_name}{$SNPs_pos}=1;
    #push @{$SNPs{$SNPs_scaffold_name}}, $SNPs_pos;
  }

  # Create the SNPs file to hold SNPs of all fragments and fragments in range when sequencing through pair-end.
  #my($SNPs_at_all_frags_fh,$SNPs_at_frags_in_range_fh);
  #if($sequence_type_code==2){
    #$SNPs_at_all_frags_fh=IO::File->new(">>$output_dir/SNPs_at_all_frags_${name_reference}_by_${enzyme}");
    #$SNPs_at_frags_in_range_fh=IO::File->new(">>$output_dir/SNPs_at_frags_in_range_${name_reference}_by_${enzyme}");
  #}

  print "Counting SNPs at the output fragments...\t";
  # Scan all locs file to count SNPs and output the result to the summary file.
  my $all_SNPs_number=0;
  if($sequence_type_code == 1){
    while(<$all_loc_file_fh>){
      chomp;
      my($fragment_name,$pos1,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$se_sequence_length-1;
        my $frag_start2=$pos2-$se_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $all_SNPs_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $all_SNPs_number++;
          }
        }
      }
    }
    $all_SNPs_number=int($all_SNPs_number/2);
    $summary_digestion_fh->print("Expected_SNPs_at_all_fragments_via_$sequence_type:\t$all_SNPs_number\n");
  }
  elsif($sequence_type_code ==2){
    while(<$all_loc_file_fh>){
      chomp;
      my($fragment_name,$pos1,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $all_SNPs_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $all_SNPs_number++;
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected_SNPs_at_all_fragments_via_$sequence_type:\t$all_SNPs_number\n");
  }
  # Scan locs in range file to count SNPs and output the result to the summary file.

  my $SNPs_in_range_number=0;

  if($sequence_type_code == 1){
    while(<$loc_in_range_file_fh>){
      chomp;
      my($fragment_name,$pos1,$pos2,$length)=split /\t/, $_;
      my $scaffold_name;
      if($fragment_name=~/^>(.+)-(\d+)$/){
        $scaffold_name=$1;
      }
      if($SNPs{$scaffold_name}){
        my $frag_start1=$pos1; my $frag_end1=$pos1+$se_sequence_length-1;
        my $frag_start2=$pos2-$se_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $SNPs_in_range_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $SNPs_in_range_number++;
          }
        }
      }
    }
    $SNPs_in_range_number=int($SNPs_in_range_number/2);
    $summary_digestion_fh->print("Expected_SNPs_at_fragments_in_range_via_$sequence_type:\t$SNPs_in_range_number\n");
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
        my $frag_start1=$pos1; my $frag_end1=$pos1+$pe_sequence_length-1;
        my $frag_start2=$pos2-$pe_sequence_length+1; my $frag_end2=$pos2;
        for my $frag_pos ($frag_start1..$frag_end1){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $SNPs_in_range_number++;
          }
        }
        for my $frag_pos ($frag_start2..$frag_end2){
          if($SNPs{$scaffold_name}{$frag_pos}){
            $SNPs_in_range_number++;
          }
        }
      }
    }
    $summary_digestion_fh->print("Expected_SNPs_at_fragments_in_range_via_$sequence_type:\t$SNPs_in_range_number\n");
  }
  print "Done!\n";
}

=head2 add_gff

    $digest->add_gff(-ref=>'Full path of the gff file');
        Adds the Gff file to the object. 
 
=cut

sub add_gff { 
  my $self=shift;
  my %parameters=@_;
  for (keys %parameters){
    if($_=~/^-gff$/){
      $self->{gff}=$parameters{$_};
    }
    else{
      die"Unacceptable parameters of the method <add_gff>, please perldoc RestrictionDigest
      for example or read README for more help\n";
    }
  }
}

=head2 all_frags_coverage

     $digest->all_frags_coverage();
       Calculate the coverage of all restriction fragments on different genome regions.
   
=cut

sub all_frags_coverage {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $enzyme=$self->{enzyme};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $loc_file_all_frags="$output_dir/position_frags_${basename_ref}_by_${enzyme}";

  $self->genome_structure_coverage($loc_file_all_frags);
}

=head2 frags_in_range_coverage

     $digest->frags_in_range_coverage();
      Calculate the coverage of restriction fragments in range on different genome regions.
=cut

sub frags_in_range_coverage {
  my $self=shift;
  
  my $ref=$self->{ref};
  my $basename_ref=basename($ref);
  
  my $enzyme=$self->{enzyme};
  
  my $output_dir=$self->{output_dir};
  $output_dir=~s/^(.+)\/?$/$1/;

  my $loc_file_all_frags="$output_dir/position_frags_in_range_${basename_ref}_by_${enzyme}";

  $self->genome_structure_coverage($loc_file_all_frags);
}

=head2  genome_structure_coverage

        A subroutine invoked by the 'all_frags_coverage' and the 'frags_in_range_coverage' methods. 
        User do not use this subroutine.

=cut

 
sub genome_structure_coverage {
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

   
  # Make  output filehandles to generic region and intergenic region coverage files. 
  my $coverage_ratio_fh;
  
  my $basename_loc_frags_file=basename($loc_frags_file);

  if($basename_loc_frags_file =~/in_range/){
    $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_in_range_${basename_ref}_by_${enzyme}");
  }
  else{
    $coverage_ratio_fh=IO::File->new(">>$output_dir/genome_coverage_${basename_ref}_by_${enzyme}");
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
  my $seq_tmp="";
  my $seq_name_tmp="";
  my $scfd_name;
  my $scfd_seq;
  while(<$ref_fh>){
    chomp;
    my $line=$_;

    if($_=~/^>(.+)/){
      my $seq_name=$1;
      my $seq_name_length=length $seq_name_tmp;
      if($seq_name_length !=0){
        $scfd_name=$seq_name_tmp;
        $scfd_seq=$seq_tmp;
        
      my $scaffold_length=length $scfd_seq;
      print "Coverage: Processing $scfd_name...\t";


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
      $seq_name_tmp=$seq_name;
      $seq_tmp="";
    }
    else{
      $seq_tmp=$seq_tmp.$line;
    }
    if(eof($ref_fh)){
      $scfd_name=$seq_name_tmp;
      $scfd_seq=$seq_tmp;
      
      my $scaffold_length=length $scfd_seq;
      print "Coverage: Processing $scfd_name...\t";


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
  
  my $all_map_length=$all_gene_map_length+$all_intergenic_map_length;
  my $frag_length_on_intergenic_rate=$all_intergenic_map_length/$all_map_length;
  my $frag_length_on_gene_rate=$all_gene_map_length/$all_map_length;
  my $frag_length_on_exon_rate=$all_exon_map_length/$all_map_length;
  my $frag_length_on_intron_rate=$all_intron_map_length/$all_map_length;

  $coverage_ratio_fh->print("Intergenic region coverage is\t$all_intergenic_map_ratio\nGenes region coverage is\t$all_gene_map_ratio\nExon region coverage  is\t$all_exon_map_ratio\nIntron region coverage is\t$all_intron_map_ratio\n");
  $coverage_ratio_fh->print("Total length of fragments:\t$all_map_length\n");
  $coverage_ratio_fh->print("Total bases mapped on intergenic regions:\t$all_intergenic_map_length\tLength rate on all fragments bases:\t$frag_length_on_intergenic_rate\n");
  $coverage_ratio_fh->print("Total bases mapped on gene regions:\t$all_gene_map_length\tLength rate on all fragments bases:\t$frag_length_on_gene_rate\n");
  $coverage_ratio_fh->print("Total bases mapped on exon regions:\t$all_exon_map_length\tLength rate on all fragments bases:\t$frag_length_on_exon_rate\n");
  $coverage_ratio_fh->print("Total bases mapped on intron regions:\t$all_intron_map_length\tLength rate on all fragments bases:\t$frag_length_on_intron_rate\n");

  }



package RestrictionDigest::MultipleItems::Double;

use 5.8.0;
use strict;
use warnings FATAL => 'all';


=head1 NAME

RestrictionDigest::MultipleItems::Double

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';



=head1 SYNOPSIS

RestrictionDigest::MultipleItems::Double is used for the simulation of double-enzyme digestion of multiple genomes
at one run. 


    use RestrictionDigest;

    $digest=RestrictionDigest::MultipleItems::Double->new();
    
    $digest->add_refs(-ref1=>'Full path of reference 1 file', -ref2=>'Full path to reference 2 file',...,-refX=>'Full path to reference X file');

    $digest->new_enzyme(-enzyme_name=>'Ncol',-recognition_site=>'C|CATGG');

    $digest->add_enzyme_pair(-front_enzyme=>'EcoRI',-behind_enzyme=>'HinfI');

    $digest->add_output_dir(-output_dir=>'Full path of the output directory');

    $digest->change_lengths_distribution_parameters(-front=>200,-behind=>800,-step=>25);

    $digest->digests_and_compare();

=cut

use IO::File;
use File::Basename;

use vars qw(%fields %enzyme_list);

# Define the fields to be used as identifiers.
%fields = (
  ref1 => undef,
  ref2 => undef,
  ref3 => undef,
  ref4 => undef,
  ref5 => undef,
  ref6 => undef,
  ref7 => undef,
  ref8 => undef,
  ref9 => undef,
  ref10 => undef,
  ref11 => undef,
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
	  MslI => 'CAYNN|NNRTG',
);

=head1 SUBROUTINES/METHODS

=head2 new
     
     $digest=RestrictionDigest::SingleItem::Double->new(); 
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

=head2 add_refs
     
     $digest->add_refs(-ref1=>'Full path of reference 1 file', -ref2=>'Full path to reference 2 file',...,-refX=>'Full path to reference X file');
      Add the reference files to the object. 
          
=cut

sub add_refs {

  my $self=shift;
  my %parameters=@_;
  my $refs_count=0;
  for my $key (keys %parameters){
    if($key=~/^-(ref\d+)$/){
      my $ref_code=$1;
      unless($parameters{$key}){
        die"No file provided for $key!\n";
      }
      $refs_count++;
      $self->{$ref_code}=$parameters{$key};
    }
    else{
    die "Unacceptable parameters in the method <add_ref>, please perldoc RestrictionDigest for
    example or read README for more help!\n";
    }
  }
  unless($refs_count>=2){
    die "At least 2 references are needed for RestrictionDigest::MultipleItems\n";
  } 
}

=head2 add_enzyme_pair
     
     $digest->add_enzyme_pair(-front_enzyme=>'EcoRI', -behind_enzyme=>'HinfI');
      Add  two enzymes to the object. The names of the two enzymes are case-insensitive.

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
      die"Unacceptable parameters in the method <add_enzyme_pair>, please perldoc RestrictionDigest
      for example or read README for more help!\n";
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
    die "The front restrict endonuclease  User provides does not exist.\nHowever User can add this
    enzyme and its recognition sites to enzyme-container via the method 'new_enzyme' .\n";
  }
  if($enzyme2_exists >0 ){
    
  }else{
    die "The behind restrict endonuclease  User provides does not exist.\nHowever User can add this
    enzyme and its recognition sites to enzyme-container via the method 'new_enzyme' .\n";
  }
}

=head2 new_enzyme

     $digest->new_enzyme(-enzyme_name=>'EcoRI', -recognition_site=>'G|AATTC');
     
     If the enzyme User wants to use does not exists in the enzyme resevior of the module, the enzyme can be added
     temporarily to the enzyme resevior.

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
      This function is used to change the resolution of fragments lengths distribution. This function has three
      parameters: front and behind define the two boundary length values, and step defines the step length.


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
         This function is used to change the resolution of fragments lengths distribution. This function has three
      parameters: front and behind define the two boundary length values, and step defines the step length.
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
      die"Unacceptable parameters in the method <change_lengths_distribution_parameters>,
      please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
}


=head2 add_output_dir
    
     $digest->add_output_dir(-output_dir=>'Full path of the output directory');
          Adds the directory into which the result files are put by the program. 

=cut


sub add_output_dir {
  my $self=shift;
  my %parameters=@_;
  for(keys %parameters){
    if($_=~/^-output_dir$/){
      $self->{output_dir}=$parameters{$_};
    }
    else{
      die "Unacceptable parameter in the method <add_output_dir>, please perldoc RestrictionDigest
      for example or read README for more help!\n";
    }
  }
  if(-d $self->{output_dir} ) {
  }else{
    die "The output directory User provides does not exist.\n";
  }
}

=head2 digests_and_compare

     $digest->digests_and_compare();
          Exexute the digestions precess.
          This process will digest all genomes/sequences inputed. Loci numbers and fragments lengths 
          distributions will be generated.

=cut

sub digests_and_compare {

  my $self=shift;
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

  my %loci_lengths;
  my @reference_names;
  my %order;

  for my $key (keys %{$self}){
    if($key=~/^ref(\d+)$/){
      my $ref_no=$1;
      my $ref_genome=$self->{$key};
      if($ref_genome){
        $order{$ref_no}=$key;
        my $ref_basename=basename($ref_genome);
        print "The ${ref_no}th reference file is:\t$ref_basename\n";
        # Process digests and store results into %loci_lengths;
        my @refs_to_data=$self->multiple_genomes_double_digest($ref_genome);
        $loci_lengths{$key}{loci}=shift @refs_to_data;
        $loci_lengths{$key}{three_types}=shift @refs_to_data;
        $loci_lengths{$key}{FB_BF}=shift @refs_to_data;
        $loci_lengths{$key}{FF}=shift @refs_to_data;
        $loci_lengths{$key}{BB}=shift @refs_to_data;
      }
    }
  }
  my $summary_digestion_fh=IO::File->new(">>$output_dir/digestion_summary_multiple_genomes_by_${front_enzyme}_and_${later_enzyme}");
  print "Compare fragments lengths distributions...\t";

  # Output header line which is the names of reference genomes.
  $summary_digestion_fh->print("Number_of_enzyme_loci");
  for my $ref_no (sort {$a<=>$b} keys %order){
    my $key=$order{$ref_no};
    my $ref_genome=$self->{$key};
    my $ref_basename=basename($ref_genome);
    $summary_digestion_fh->print("\t$ref_basename");
  }
  # Output front enzyme loci number.
  $summary_digestion_fh->print("\nLoci_number_of_$front_enzyme");
  for my $ref_no(sort {$a<=>$b} keys %order){
    my $key=$order{$ref_no};
    my $ref_to_loci=$loci_lengths{$key}{loci};
    my @two_enzymes_loci_no=@{$ref_to_loci};
    my $front_enzyme_loci=$two_enzymes_loci_no[0];
    $summary_digestion_fh->print("\t$front_enzyme_loci");
    
  }
  # Output later enzyme loci number.
  $summary_digestion_fh->print("\nLoci_number_of_$later_enzyme");
  for my $ref_no (sort {$a<=>$b} keys %order){
    my $key=$order{$ref_no};
    my $ref_to_loci=$loci_lengths{$key}{loci};
    my @two_enzymes_loci_no=@{$ref_to_loci};
    my $later_enzyme_loci=$two_enzymes_loci_no[1];
    $summary_digestion_fh->print("\t$later_enzyme_loci");
  }

  # Output header of lengths distributions.
  $summary_digestion_fh->print("\nLength_range");
  for my $ref_no (sort {$a<=>$b} keys %order){
    my $key=$order{$ref_no};
    my $ref_genome=$self->{$key};
    my $ref_basename=basename($ref_genome);
    $summary_digestion_fh->print("\tthree_types_number_$ref_basename\tthree_types_ratio_$ref_basename");
    $summary_digestion_fh->print("\tFB_BF_number_$ref_basename\tFB_BF_ratio_$ref_basename");
    $summary_digestion_fh->print("\tFF_number_$ref_basename\tFF_ratio_$ref_basename");
    $summary_digestion_fh->print("\tBB_number_$ref_basename\tBB_ratio_$ref_basename");
  }
  $summary_digestion_fh->print("\n");

  # Parsing distribution parameters and initialize all values.
  my $lengths_distribution_start=$self->{lengths_distribution_start};
  my $lengths_distribution_end=$self->{lengths_distribution_end};
  my $lengths_distribution_step=$self->{lengths_distribution_step};
  my $lengths_distribution_tmp;
  my (@starts,@ends,@counts,%ref_genomes_depot,%description_copy);

  my $num_pair=0;
  for($lengths_distribution_tmp= $lengths_distribution_start;$lengths_distribution_tmp<$lengths_distribution_end+$lengths_distribution_step;$lengths_distribution_tmp=$lengths_distribution_tmp+$lengths_distribution_step){
    if($lengths_distribution_tmp== $lengths_distribution_start){
      $num_pair++;
      push @starts, '0';
      push @ends , $lengths_distribution_tmp;

      for my $ref_code (keys %loci_lengths){
        $ref_genomes_depot{$ref_code}{three_types}{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
        $ref_genomes_depot{$ref_code}{FB_BF}{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
        $ref_genomes_depot{$ref_code}{FF}{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
        $ref_genomes_depot{$ref_code}{BB}{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
      }
      $description_copy{$num_pair}="0bp-${lengths_distribution_tmp}bp";

      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      for my $ref_code (keys %loci_lengths){
        $ref_genomes_depot{$ref_code}{three_types}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        $ref_genomes_depot{$ref_code}{FB_BF}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        $ref_genomes_depot{$ref_code}{FF}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        $ref_genomes_depot{$ref_code}{BB}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      }
      $description_copy{$num_pair}="${tmp_start}bp-${tmp_end}bp";
    }
    elsif($lengths_distribution_tmp  < $lengths_distribution_end){
      if($lengths_distribution_tmp+$lengths_distribution_step <$lengths_distribution_end ){
        $num_pair++;
        my $tmp_start=$lengths_distribution_tmp+1;
        my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
        push @starts, $tmp_start;
        push @ends, $tmp_end;
        for my $ref_code (keys %loci_lengths){
          $ref_genomes_depot{$ref_code}{three_types}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
          $ref_genomes_depot{$ref_code}{FB_BF}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
          $ref_genomes_depot{$ref_code}{FF}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
          $ref_genomes_depot{$ref_code}{BB}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        }
        $description_copy{$num_pair}="${tmp_start}bp-${tmp_end}bp";
      }
      elsif($lengths_distribution_tmp +$lengths_distribution_step >= $lengths_distribution_end){
        $num_pair++;
        my $tmp_start=$lengths_distribution_tmp+1;
        my $tmp_end=$lengths_distribution_end;

        push @starts, $tmp_start;
        push @ends, $tmp_end;

        for my $ref_code (keys %loci_lengths){
          $ref_genomes_depot{$ref_code}{three_types}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
          $ref_genomes_depot{$ref_code}{FB_BF}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
          $ref_genomes_depot{$ref_code}{FF}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
          $ref_genomes_depot{$ref_code}{BB}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        }
        $description_copy{$num_pair}="${tmp_start}bp-${tmp_end}bp";
      
        $num_pair++;
        $tmp_start=$tmp_end+1;

        push @starts, $tmp_start;
        push @ends, 'bigger';
        for my $ref_code (keys %loci_lengths){
          $ref_genomes_depot{$ref_code}{three_types}{$num_pair}{"${tmp_start}bp-longer"}=0;
          $ref_genomes_depot{$ref_code}{FB_BF}{$num_pair}{"${tmp_start}bp-longer"}=0;
          $ref_genomes_depot{$ref_code}{FF}{$num_pair}{"${tmp_start}bp-longer"}=0;
          $ref_genomes_depot{$ref_code}{BB}{$num_pair}{"${tmp_start}bp-longer"}=0;
        }
        $description_copy{$num_pair}="${tmp_start}bp-longer";
      }
    }
  }

  # Analyze the distributions of lengths of all reference genomes.
  for my $ref_code (keys %loci_lengths){
    for my $frag_type (keys %{$loci_lengths{$ref_code}}){
      unless($frag_type eq "loci"){
        my @lengths_depot=@{$loci_lengths{$ref_code}{$frag_type}};
        for (@lengths_depot){
          my $length=$_;
          my $num_pair=1;
          my $num_last_pair=@starts;
          for($num_pair=1;$num_pair<=$num_last_pair;$num_pair++){
            my $array_index=$num_pair-1;
            my $start=$starts[$array_index];
            my $end=$ends[$array_index];
            if($num_pair<$num_last_pair){
              if($length>=$start && $length<=$end){
                for my $description (keys %{$ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}}){
                  $ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}{$description}++;
                  last;
                }
                last;
              }
            }
            elsif($num_pair==$num_last_pair){
              if($length>=$start){
                for my $description (keys %{$ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}}){
                  $ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}{$description}++;
                  last;
                }
                last;
              }
            }
          }
        }
      }
    }
  }

  # Output details of lengths distributions.
  for (1..@starts){
    my $num_pair=$_;
    my $description_copy=$description_copy{$num_pair};
    $summary_digestion_fh->print("$description_copy");
    for my $ref_no(sort {$a<=>$b} keys %order){
      my $key=$order{$ref_no};
      for my $description (keys %{$ref_genomes_depot{$key}{three_types}{$num_pair}}){
        my $three_types_fragments_num=$ref_genomes_depot{$key}{three_types}{$num_pair}{$description};
        my $FB_BF_fragments_num=$ref_genomes_depot{$key}{FB_BF}{$num_pair}{$description};
        my $FF_fragments_num=$ref_genomes_depot{$key}{FF}{$num_pair}{$description};
        my $BB_fragments_num=$ref_genomes_depot{$key}{BB}{$num_pair}{$description};
        my ($three_types_ratio,$FB_BF_ratio,$FF_ratio,$BB_ratio);
        my $three_types_count=@{$loci_lengths{$key}{three_types}};
        my $FB_BF_count=@{$loci_lengths{$key}{FB_BF}};
        my $FF_count=@{$loci_lengths{$key}{FF}};
        my $BB_count=@{$loci_lengths{$key}{BB}};
        if($three_types_count){
          $three_types_ratio=$three_types_fragments_num/$three_types_count;
        }else{
          $three_types_ratio="N/A";
        }
        if($FB_BF_count){
          $FB_BF_ratio=$FB_BF_fragments_num/$FB_BF_count;
        }else{
          $FB_BF_ratio="N/A";
        } 
        if($FF_count){
          $FF_ratio=$FF_fragments_num/$FF_count;
        }else{
          $FF_ratio="N/A";
        }
        if($BB_count){
          $BB_ratio=$BB_fragments_num/$BB_count;
        }else{
          $BB_ratio="N/A";
        }
        $summary_digestion_fh->print("\t$three_types_fragments_num\t$three_types_ratio");
        $summary_digestion_fh->print("\t$FB_BF_fragments_num\t$FB_BF_ratio");
        $summary_digestion_fh->print("\t$FF_fragments_num\t$FF_ratio");
        $summary_digestion_fh->print("\t$BB_fragments_num\t$BB_ratio");
      }
    }
    $summary_digestion_fh->print("\n");
  }
  print "Done!\n";
  

}          

=head2 multiple_genomes_double_digest

     $digest->multiple_genomes_double_digest();
          A subroutine invoked by RestrictionDigest::MultipleItems::Double.

=cut

sub multiple_genomes_double_digest {
  my $self=shift;
  my $ref=shift;
  my $front_enzyme=$self->{enzyme1};
  my $later_enzyme=$self->{enzyme2};
  my $range1=$self->{range_start};
  my $range2=$self->{range_end};
  my $output_dir=$self->{output_dir};
  
  $output_dir=~s/^(.+)\/?$/$1/;
  
  # Make the file handle to the reference file.
  my $ref_fh=IO::File->new("$ref",'r');
  
  # Get the basename of reference,not include the path.
  my $name_reference=basename($ref);
  
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
  my @lengths_of_FB_BF_fragments=();
  my @lengths_of_FF_fragments=();
  my @lengths_of_BB_fragments=();

  my @lengths_fragments_in_range=();
  my @lengths_FB_BF_fragments_in_range=();
  my @lengths_FF_fragments_in_range=();
  my @lengths_BB_fragments_in_range=();

  my $overall_fragments_length=0;
  my $overall_FB_BF_fragments_length=0;
  my $overall_FF_fragments_length=0;
  my $overall_BB_fragments_length=0;

  my $overall_fragments_in_range_length=0;
  my $overall_FB_BF_fragments_in_range_length=0;
  my $overall_FF_fragments_in_range_length=0;
  my $overall_BB_fragments_in_range_length=0;

  my $overall_length_of_scfd=0;

  my $total_front_enzyme_loci;
  my $total_later_enzyme_loci;

  ########## CIRCLES START HERE!###########

  # Process one scaffold one circle;
  my $seq_tmp="";
  my $seq_name_tmp="";
  my $scaffold_name; my $scaffold_seq;
  print "Digesting $name_reference...\t";
  while(<$ref_fh>){
    chomp;
    my $line=$_;
    if($_=~/^>/){
      my $seq_name=$line;
      my $seq_name_length=length $seq_name_tmp;
      if($seq_name_length !=0){
        $scaffold_name=$seq_name_tmp;
        $scaffold_seq=$seq_tmp;
		$scaffold_seq = uc $scaffold_seq;
        #print"Digesting $name_reference:$scaffold_name...\t";
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

		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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


		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }
     
        # Count all enzymes loci of the genome.
        $total_front_enzyme_loci+=@front_enzyme_locs;
        $total_later_enzyme_loci+=@later_enzyme_locs;

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

        # Three types of fragments: FB_BF, FF, BB (F for Front_enzyme and B for Behind_enzyme)
        my ($length_FB_BF_fragments,$length_FF_fragments,$length_BB_fragments,$length_FB_BF_fragments_in_range)=qw(0 0 0 0);
        my ($length_FF_fragments_in_range, $length_BB_fragments_in_range)=qw(0 0);
        my ($count_all,$count_in_range,$count_FB_BF,$count_FF,$count_BB,$count_FB_BF_in_range,$count_FF_in_range,$count_BB_in_range)=qw(0 0 0 0 0 0 0 0);
        my ($before_cut_coorr,$after_cut_coorr,$bf_ct_crr_hm_rd,$af_ct_crr_hm_rd);

        my (@FB_BF_locs_before, @FB_BF_enzyme_before, @FB_BF_locs_after, @FB_BF_enzyme_after);
        my (@FF_locs_before,  @FF_locs_after,@BB_locs_before, @BB_locs_after);
        my ($after_loc,$enzyme_before, $enzyme_after);
        my $before_loc=shift @locs_sorted;

        while(@locs_sorted){
          $after_loc=shift @locs_sorted;
          $enzyme_before=$hash{$before_loc};
          $enzyme_after=$hash{$after_loc};
          unless($enzyme_before eq $enzyme_after){
            push @FB_BF_locs_before, $before_loc;
            push @FB_BF_enzyme_before, $enzyme_before;
            push @FB_BF_locs_after,$after_loc;
            push @FB_BF_enzyme_after, $enzyme_after;
            my $length_of_FB_BF_fragment=0; my $seq_FB_BF;
            $count_all++; $count_FB_BF++;
            if($front_enzyme=~/$enzyme_before/i){
              $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
              $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
              $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
              $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);   
            }
            elsif($later_enzyme=~/$enzyme_before/i){
                $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
                $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
                $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
                $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);
            }

            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;

            # Output the sequence and fragment positions of all_length_FB_BF fragments to files;
            #$all_FB_BF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
            my $strand;
            if($enzyme_before=~/$front_enzyme/i){$strand="+";}
            elsif($enzyme_before=~/$later_enzyme/i){$strand="-";}
            #$all_FB_BF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
            push @lengths_of_fragments, $length_of_FB_BF_fragment;
            push @lengths_of_FB_BF_fragments, $length_of_FB_BF_fragment; 

            # Process fragments in range. 
            if($length_of_FB_BF_fragment >= $range1 && $length_of_FB_BF_fragment <= $range2){
              $count_in_range++;
              $count_FB_BF_in_range++;
              #$seq_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
              #$loc_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FB_BF_fragment;
              push @lengths_FB_BF_fragments_in_range, $length_of_FB_BF_fragment;
        
              $length_FB_BF_fragments_in_range+=$length_of_FB_BF_fragment;              
            }

            # Sum all fragments length no matter in range or not!!!
            
            $length_FB_BF_fragments+=$length_of_FB_BF_fragment;
          }elsif($enzyme_before eq $front_enzyme){
            push @FF_locs_before, $before_loc;
            push @FF_locs_after, $after_loc;
            my $length_of_FF_fragment=0; my $seq_FF;
            $count_all++; $count_FF++;
            $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
            $length_of_FF_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_FF=substr($scaffold_seq, $before_cut_coorr, $length_of_FF_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            #$all_FF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
            my $strand="+";
            #$all_FF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
            push @lengths_of_fragments, $length_of_FF_fragment;
            push @lengths_of_FF_fragments, $length_of_FF_fragment;

            # Process fragments in range.
            if($length_of_FF_fragment >= $range1 && $length_of_FF_fragment <= $range2){
              $count_in_range++;
              $count_FF_in_range++;
              #$seq_FF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
              #$loc_FF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FF_fragment;
              push @lengths_FF_fragments_in_range, $length_of_FF_fragment;
             
              $length_FF_fragments_in_range+=$length_of_FF_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
            
            $length_FF_fragments+=$length_of_FF_fragment;

          }elsif($enzyme_before eq $later_enzyme){
            push @BB_locs_before, $before_loc;
            push @BB_locs_after, $after_loc;
            my $length_of_BB_fragment=0; my $seq_BB;
            $count_all++; $count_BB++;
            $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
            $length_of_BB_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_BB=substr($scaffold_seq, $before_cut_coorr, $length_of_BB_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            #$all_BB_seq_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
            my $strand="+";
            #$all_BB_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
            push @lengths_of_fragments, $length_of_BB_fragment;
            push @lengths_of_BB_fragments, $length_of_BB_fragment;

            # Process fragments in range.
            if($length_of_BB_fragment >= $range1 && $length_of_BB_fragment <= $range2){
              $count_in_range++;
              $count_BB_in_range++;
              #$seq_BB_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
              #$loc_BB_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
              push @lengths_fragments_in_range, $length_of_BB_fragment;
              push @lengths_BB_fragments_in_range, $length_of_BB_fragment;
              
              $length_BB_fragments_in_range+=$length_of_BB_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
          
            $length_BB_fragments+=$length_of_BB_fragment;

          }
          $before_loc=$after_loc;
        }

        # Calculate reduced rate of three types of fragments on this chromosome.
        my $all_type_length_overall=$length_FB_BF_fragments+$length_FF_fragments+$length_BB_fragments;
        my $all_type_length_in_range=$length_FB_BF_fragments_in_range+$length_FF_fragments_in_range+$length_BB_fragments_in_range;
        my $all_type_overall_rate=$all_type_length_overall/$length_of_scafd;
        my $all_type_in_range_rate=$all_type_length_in_range/$length_of_scafd;
        my $FB_BF_overall_rate=$length_FB_BF_fragments/$length_of_scafd;
        my $FB_BF_in_range_rate=$length_FB_BF_fragments_in_range/$length_of_scafd;
        my $FF_overall_rate=$length_FF_fragments/$length_of_scafd;
        my $FF_in_range_rate=$length_FF_fragments_in_range/$length_of_scafd;
        my $BB_overall_rate=$length_BB_fragments/$length_of_scafd;
        my $BB_in_range_rate=$length_BB_fragments_in_range/$length_of_scafd;

        $overall_fragments_length+=$all_type_length_overall;
        $overall_fragments_in_range_length+=$all_type_length_in_range;
        $overall_length_of_scfd+=$length_of_scafd;

        $overall_FB_BF_fragments_length+=$length_FB_BF_fragments;
        $overall_FF_fragments_length+=$length_FF_fragments;
        $overall_BB_fragments_length+=$length_BB_fragments;

        $overall_FB_BF_fragments_in_range_length+=$length_FB_BF_fragments_in_range;
        $overall_FF_fragments_in_range_length+=$length_FF_fragments_in_range;
        $overall_BB_fragments_in_range_length+=$length_BB_fragments_in_range;

        #$reduced_ratio_every_chromosome_file_fh->print("$all_type_overall_rate\t$all_type_in_range_rate\t$FB_BF_overall_rate\t$FB_BF_in_range_rate\t");
        #$reduced_ratio_every_chromosome_file_fh->print("$FF_overall_rate\t$FF_in_range_rate\t$BB_overall_rate\t$BB_in_range_rate\n");
        #print"Done!\n";
      }
      $seq_name_tmp=$seq_name;
      $seq_tmp="";
    }else{
      $seq_tmp=$seq_tmp.$line;
    }
    if(eof($ref_fh)){
      $scaffold_name=$seq_name_tmp;
      $scaffold_seq=$seq_tmp;
	  $scaffold_seq = uc $scaffold_seq;
        #print"Digesting $name_reference:$scaffold_name...\t";
        #$reduced_ratio_every_chromosome_file_fh->print("$scaffold_name\t");
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

		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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

		  elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }



        
        # Count all enzymes loci of the genome.
        $total_front_enzyme_loci+=@front_enzyme_locs;
        $total_later_enzyme_loci+=@later_enzyme_locs;
     
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
        

        # Three types of fragments: FB_BF, FF, BB (F for Front_enzyme and B for Behind_enzyme)
        my ($length_FB_BF_fragments,$length_FF_fragments,$length_BB_fragments,$length_FB_BF_fragments_in_range)=qw(0 0 0 0);
        my ($length_FF_fragments_in_range, $length_BB_fragments_in_range)=qw(0 0);
        my ($count_all,$count_in_range,$count_FB_BF,$count_FF,$count_BB,$count_FB_BF_in_range,$count_FF_in_range,$count_BB_in_range)=qw(0 0 0 0 0 0 0 0);
        my ($before_cut_coorr,$after_cut_coorr,$bf_ct_crr_hm_rd,$af_ct_crr_hm_rd);

        my (@FB_BF_locs_before, @FB_BF_enzyme_before, @FB_BF_locs_after, @FB_BF_enzyme_after);
        my (@FF_locs_before,  @FF_locs_after,@BB_locs_before, @BB_locs_after);
        my ($after_loc,$enzyme_before, $enzyme_after);
        my $before_loc=shift @locs_sorted;

        while(@locs_sorted){
          $after_loc=shift @locs_sorted;
          $enzyme_before=$hash{$before_loc};
          $enzyme_after=$hash{$after_loc};
          unless($enzyme_before eq $enzyme_after){
            push @FB_BF_locs_before, $before_loc;
            push @FB_BF_enzyme_before, $enzyme_before;
            push @FB_BF_locs_after,$after_loc;
            push @FB_BF_enzyme_after, $enzyme_after;
            my $length_of_FB_BF_fragment=0; my $seq_FB_BF;
            $count_all++; $count_FB_BF++;
            if($front_enzyme=~/$enzyme_before/i){
              $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
              $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
              $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
              $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);   
            }
            elsif($later_enzyme=~/$enzyme_before/i){
                $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
                $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
                $length_of_FB_BF_fragment=$after_cut_coorr-$before_cut_coorr+1;
                $seq_FB_BF=substr($scaffold_seq, $before_cut_coorr, $length_of_FB_BF_fragment);
            }

            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;

            # Output the sequence and fragment positions of all_length_FB_BF fragments to files;
            #$all_FB_BF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
            my $strand;
            if($enzyme_before=~/$front_enzyme/i){$strand="+";}
            elsif($enzyme_before=~/$later_enzyme/i){$strand="-";}
            #$all_FB_BF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
            push @lengths_of_fragments, $length_of_FB_BF_fragment;
            push @lengths_of_FB_BF_fragments, $length_of_FB_BF_fragment; 

            # Process fragments in range. 
            if($length_of_FB_BF_fragment >= $range1 && $length_of_FB_BF_fragment <= $range2){
              $count_in_range++;
              $count_FB_BF_in_range++;
              #$seq_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FB_BF\n");
              #$loc_FB_BF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$enzyme_before\t$bf_ct_crr\t$enzyme_after\t$af_ct_crr\t$length_of_FB_BF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FB_BF_fragment;
              push @lengths_FB_BF_fragments_in_range, $length_of_FB_BF_fragment;
              
              $length_FB_BF_fragments_in_range+=$length_of_FB_BF_fragment;              
            }

            # Sum all fragments length no matter in range or not!!!
            
            $length_FB_BF_fragments+=$length_of_FB_BF_fragment;
          }elsif($enzyme_before eq $front_enzyme){
            push @FF_locs_before, $before_loc;
            push @FF_locs_after, $after_loc;
            my $length_of_FF_fragment=0; my $seq_FF;
            $count_all++; $count_FF++;
            $before_cut_coorr=$before_loc+$front_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$front_enzyme_cut_loc-2;
            $length_of_FF_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_FF=substr($scaffold_seq, $before_cut_coorr, $length_of_FF_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            #$all_FF_seq_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
            my $strand="+";
            #$all_FF_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
            push @lengths_of_fragments, $length_of_FF_fragment;
            push @lengths_of_FF_fragments, $length_of_FF_fragment;

            # Process fragments in range.
            if($length_of_FF_fragment >= $range1 && $length_of_FF_fragment <= $range2){
              $count_in_range++;
              $count_FF_in_range++;
              #$seq_FF_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_FF\n");
              #$loc_FF_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$front_enzyme\t$bf_ct_crr\t$front_enzyme\t$af_ct_crr\t$length_of_FF_fragment\n");
              push @lengths_fragments_in_range, $length_of_FF_fragment;
              push @lengths_FF_fragments_in_range, $length_of_FF_fragment;
           
              $length_FF_fragments_in_range+=$length_of_FF_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
            
            $length_FF_fragments+=$length_of_FF_fragment;

          }elsif($enzyme_before eq $later_enzyme){
            push @BB_locs_before, $before_loc;
            push @BB_locs_after, $after_loc;
            my $length_of_BB_fragment=0; my $seq_BB;
            $count_all++; $count_BB++;
            $before_cut_coorr=$before_loc+$later_enzyme_cut_loc-1;
            $after_cut_coorr=$after_loc+$later_enzyme_cut_loc-2;
            $length_of_BB_fragment=$after_cut_coorr-$before_cut_coorr+1;
            $seq_BB=substr($scaffold_seq, $before_cut_coorr, $length_of_BB_fragment);
            my $bf_ct_crr=$before_cut_coorr+1; my $af_ct_crr=$after_cut_coorr+1;
            # Output the sequence and fragment positions of all_length_FF fragments to files;
            #$all_BB_seq_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
            my $strand="+";
            #$all_BB_loc_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
            push @lengths_of_fragments, $length_of_BB_fragment;
            push @lengths_of_BB_fragments, $length_of_BB_fragment;

            # Process fragments in range.
            if($length_of_BB_fragment >= $range1 && $length_of_BB_fragment <= $range2){
              $count_in_range++;
              $count_BB_in_range++;
              #$seq_BB_in_range_file_fh->print("$scaffold_name-$count_all\n$seq_BB\n");
              #$loc_BB_in_range_file_fh->print("$scaffold_name-$count_all\t$strand\t$later_enzyme\t$bf_ct_crr\t$later_enzyme\t$af_ct_crr\t$length_of_BB_fragment\n");
              push @lengths_fragments_in_range, $length_of_BB_fragment;
              push @lengths_BB_fragments_in_range, $length_of_BB_fragment;
             
              $length_BB_fragments_in_range+=$length_of_BB_fragment;              
            }
            # Sum all fragments length no matter in range or not!!!
           
            $length_BB_fragments+=$length_of_BB_fragment;

          }
          $before_loc=$after_loc;
        }
        #print"Done!\n";
    }
  }
  print "Done!\n";
  my @two_enzymes_loci_no=($total_front_enzyme_loci,$total_later_enzyme_loci);
  return (\@two_enzymes_loci_no, \@lengths_of_fragments, \@lengths_of_FB_BF_fragments, \@lengths_of_FF_fragments,\@lengths_of_BB_fragments);
}



package RestrictionDigest::MultipleItems::Single;

use 5.8.0;
use strict;
use warnings FATAL => 'all';

=head1 NAME


RestrictionDigest::MultipleItems::Single
=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';



=head1 SYNOPSIS

RestrictionDigest::MultipleItems::Single is used for the simulation of single-enzyme digestion of multiple genomes.

    use RestrictionDigest;

    $digest=RestrictionDigest::MultipleItems::Single->new();
    
    $digest->add_refs(-ref1=>'Full path of reference 1 file', -ref2=>'Full path to reference 2 file',...,-refX=>'Full path to reference X file');

    $digest->new_enzyme(-enzyme_name=>'Ncol',-recognition_site=>'C|CATGG');

    $digest->add_enzyme_pair(-front_enzyme=>'EcoRI',-behind_enzyme=>'HinfI');

    $digest->add_output_dir(-output_dir=>'Full path of the output directory');

    $digest->change_lengths_distribution_parameters(-front=>200,-behind=>800,-step=>25);

    $digest->digests_and_compare();

=cut

use IO::File;
use File::Basename;

use vars qw(%fields %enzyme_list);

# Define the fields to be used as identifiers.
%fields = (
  ref1 => undef,
  ref2 => undef,
  ref3 => undef,
  ref4 => undef,
  ref5 => undef,
  ref6 => undef,
  ref7 => undef,
  ref8 => undef,
  ref9 => undef,
  ref10 => undef,
  ref11 => undef,
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
	  MslI => 'CAYNN|NNRTG',
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
     
     $digest->add_refs(-ref1=>'Full path of reference 1 file', -ref2=>'Full path to reference 2 file',...,-refX=>'Full path to reference X file');
      Add the reference files to the object.

=cut

sub add_refs {

  my $self=shift;
  my %parameters=@_;
  my $refs_count=0;
  for my $key (keys %parameters){
    if($key=~/^-(ref\d+)$/){
      my $ref_code=$1;
      unless($parameters{$key}){
        die"No file provided for $key!\n";
      }
      $refs_count++;
      $self->{$ref_code}=$parameters{$key};
    }
    else{
    die "Unacceptable parameters in the method <add_ref>, please perldoc RestrictionDigest for example or read README for more help!\n";
    }
  }
  unless($refs_count>=2){
    die "At least 2 references are needed for RestrictionDigest::MultipleItems\n";
  } 
}

=head2 add_single_enzyme
     
     $digest->add_single_enzyme(-enzyme=>'EcoRI');
          Add the single enzyme to the object.
          
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
     
     If the enzyme User wants to use does not exists in the enzyme resevior of the module, the enzyme can be added
     temporarily to the enzyme resevior.


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
        This function is used to change the resolution of fragments lengths distribution. This function has three
      parameters: front and behind define the two boundary length values, and step defines the step length.
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
        This function is used to change the resolution of fragments lengths distribution. This function has three
      parameters: front and behind define the two boundary length values, and step defines the step length.
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

=head2 digests_and_compare

     $digest->digests_and_compare();
          Exexute the digestions precess.
          This process will digest all genomes/sequences inputed. Loci numbers and fragments lengths 
          distributions will be generated.

=cut

sub digests_and_compare {

  my $self=shift;
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

  my %loci_lengths;
  my @reference_names;
  my %order;

  for my $key (keys %{$self}){
    if($key=~/^ref(\d+)$/){
      my $ref_no=$1;
      my $ref_genome=$self->{$key};
      if($ref_genome){
        $order{$ref_no}=$key;
        my $ref_basename=basename($ref_genome);
        print "The ${ref_no}th reference file is:\t$ref_basename\n";
        # Process digests and store results into %loci_lengths;

        my @refs_to_data=$self->multiple_genomes_single_digest($ref_genome);

        $loci_lengths{$key}{loci}=shift @refs_to_data;
        $loci_lengths{$key}{lengths}=shift @refs_to_data;
      }
    }
  }

  my $summary_digestion_fh=IO::File->new(">>$output_dir/digestion_summary_multiple_genomes_by_$enzyme");
  print "Compare fragments lengths distributions...\t";

  # Output header line which is the names of reference genomes.
  $summary_digestion_fh->print("Number_of_enzyme_loci");
  for my $ref_no (sort {$a<=>$b} keys %order){
    my $key=$order{$ref_no};
    my $ref_genome=$self->{$key};
    my $ref_basename=basename($ref_genome);
    $summary_digestion_fh->print("\t$ref_basename");
  }
  # Output enzyme loci number.
  $summary_digestion_fh->print("\nLoci_number_of_$enzyme");
  for my $ref_no(sort {$a<=>$b} keys %order){
    my $key=$order{$ref_no};
    my $enzyme_loci=$loci_lengths{$key}{loci};
    $summary_digestion_fh->print("\t$enzyme_loci");
  }
  
  # Output header of lengths distributions.
  $summary_digestion_fh->print("\nLength_range");
  for my $ref_no (sort {$a<=>$b} keys %order){
    my $key=$order{$ref_no};
    my $ref_genome=$self->{$key};
    my $ref_basename=basename($ref_genome);
    $summary_digestion_fh->print("\tfragments_number_$ref_basename\tfragments_ratio_$ref_basename");
  }
  $summary_digestion_fh->print("\n");

  # Parsing distribution parameters and initialize all values.
  my $lengths_distribution_start=$self->{lengths_distribution_start};
  my $lengths_distribution_end=$self->{lengths_distribution_end};
  my $lengths_distribution_step=$self->{lengths_distribution_step};
  my $lengths_distribution_tmp;
  my (@starts,@ends,@counts,%ref_genomes_depot,%description_copy);
  my $num_pair=0;
  for($lengths_distribution_tmp= $lengths_distribution_start;$lengths_distribution_tmp<$lengths_distribution_end+$lengths_distribution_step;$lengths_distribution_tmp=$lengths_distribution_tmp+$lengths_distribution_step){
    if($lengths_distribution_tmp== $lengths_distribution_start){
      $num_pair++;
      push @starts, '0';
      push @ends , $lengths_distribution_tmp;
      for my $ref_code (keys %loci_lengths){
        $ref_genomes_depot{$ref_code}{lengths}{$num_pair}{"0bp-${lengths_distribution_tmp}bp"}=0;
      }
      $description_copy{$num_pair}="0bp-${lengths_distribution_tmp}bp";

      $num_pair++;
      my $tmp_start=$lengths_distribution_tmp+1;
      my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
      push @starts, $tmp_start;
      push @ends, $tmp_end;
      for my $ref_code (keys %loci_lengths){
        $ref_genomes_depot{$ref_code}{lengths}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
      }
      $description_copy{$num_pair}="${tmp_start}bp-${tmp_end}bp";
    }
    elsif($lengths_distribution_tmp  < $lengths_distribution_end){
      if($lengths_distribution_tmp+$lengths_distribution_step <$lengths_distribution_end ){
        $num_pair++;
        my $tmp_start=$lengths_distribution_tmp+1;
        my $tmp_end=$lengths_distribution_tmp+$lengths_distribution_step;
        push @starts, $tmp_start;
        push @ends, $tmp_end;
        for my $ref_code (keys %loci_lengths){
          $ref_genomes_depot{$ref_code}{lengths}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        }
        $description_copy{$num_pair}="${tmp_start}bp-${tmp_end}bp";
      }
      elsif($lengths_distribution_tmp +$lengths_distribution_step >= $lengths_distribution_end){
        $num_pair++;
        my $tmp_start=$lengths_distribution_tmp+1;
        my $tmp_end=$lengths_distribution_end;
        push @starts, $tmp_start;
        push @ends, $tmp_end;
        for my $ref_code (keys %loci_lengths){
          $ref_genomes_depot{$ref_code}{lengths}{$num_pair}{"${tmp_start}bp-${tmp_end}bp"}=0;
        }
        $description_copy{$num_pair}="${tmp_start}bp-${tmp_end}bp";
        $num_pair++;
        $tmp_start=$tmp_end+1;
        push @starts, $tmp_start;
        push @ends, 'bigger';
        for my $ref_code (keys %loci_lengths){
          $ref_genomes_depot{$ref_code}{lengths}{$num_pair}{"${tmp_start}bp-longer"}=0;
        }
        $description_copy{$num_pair}="${tmp_start}bp-longer"; 
      }
    }
  }
  # Analyze the distributions of lengths of all reference genomes.
  for my $ref_code (keys %loci_lengths){
    for my $frag_type (keys %{$loci_lengths{$ref_code}}){
      unless($frag_type eq "loci"){
        my @lengths_depot=@{$loci_lengths{$ref_code}{$frag_type}};
        for (@lengths_depot){
          my $length=$_;
          my $num_pair=1;
          my $num_last_pair=@starts;
          for($num_pair=1;$num_pair<=$num_last_pair;$num_pair++){
            my $array_index=$num_pair-1;
            my $start=$starts[$array_index];
            my $end=$ends[$array_index];
            if($num_pair<$num_last_pair){
              if($length>=$start && $length<=$end){
                for my $description (keys %{$ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}}){
                  $ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}{$description}++;
                  last;
                }
                last;
              }
            }
            elsif($num_pair==$num_last_pair){
              if($length>=$start){
                for my $description (keys %{$ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}}){
                  $ref_genomes_depot{$ref_code}{$frag_type}{$num_pair}{$description}++;
                  last;
                }
                last;
              }
            }
          }
        }
      }
    }
  }

  # Output details of lengths distributions.
  for (1..@starts){
    my $num_pair=$_;
    my $description_copy=$description_copy{$num_pair};
    $summary_digestion_fh->print("$description_copy");
    for my $ref_no(sort {$a<=>$b} keys %order){
      my $key=$order{$ref_no};
      for my $description (keys %{$ref_genomes_depot{$key}{lengths}{$num_pair}}){
        my $fragments_num=$ref_genomes_depot{$key}{lengths}{$num_pair}{$description};
        my ($fragments_ratio);
        my $fragments_count=@{$loci_lengths{$key}{lengths}};
        if($fragments_count){
          $fragments_ratio=$fragments_num/$fragments_count;
        }
        else{
          $fragments_ratio="N/A";
        }
        $summary_digestion_fh->print("\t$fragments_num\t$fragments_ratio");
      }
    }
    $summary_digestion_fh->print("\n");
  }
  print "Done!\n";
}

  

=head2 multiple_genomes_single_digest

     A subroutine invoked by RestricitonDigest::MultipleItems::Single;

=cut

sub multiple_genomes_single_digest {

  my $self=shift;
  my $ref=shift;
  my $enzyme=$self->{enzyme};
  my $range1=$self->{range_start};
  my $range2=$self->{range_end};
  my $output_dir=$self->{output_dir};

  $output_dir=~s/^(.+)\/?$/$1/;
  
  # Make the file handle to the reference file.
  my $ref_fh=IO::File->new("$ref",'r');
  
  # Get the basename of reference,not include the path.
  my $name_reference=basename($ref);
  
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

  my $total_enzyme_loci=0;

  ########## CIRCLES START HERE!###########

  # Process one scaffold one circle;
  my $seq_tmp="";
  my $seq_name_tmp="";
  my ($scaffold_name, $scaffold_seq);
  print"Digesting $name_reference...\t";
  
  while(<$ref_fh>){
    chomp;
    my $line=$_;

    if($_=~/^>/){
      my $seq_name=$line;
      my $seq_name_length=length $seq_name_tmp;
      
      if($seq_name_length !=0){
        $scaffold_name=$seq_name_tmp;
        $scaffold_seq=$seq_tmp;
		$scaffold_seq = uc $scaffold_seq;

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
        elsif($enzyme_seq=~/CCNNNNNNN/i){
			  my $loc_in_array;
			  $loc_in_array=index($seq_for_search,'CC');
			  unless($loc_in_array==-1){
				  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				  if(length $seq_after_CC >19){
					push @enzyme_locs, $loc_in_array;
				  }
			  }
			  while($loc_in_array!=-1){
				  $loc_in_array=index($seq_for_search, 'CC', $loc_in_array+1);
				  unless($loc_in_array==-1){
					  my $previous_loc=$enzyme_locs[$#enzyme_locs];
					  if(($loc_in_array- $previous_loc) > 13){
						  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				          if(length $seq_after_CC >19){
							  push @enzyme_locs, $loc_in_array;
						  }
					  }
				  }
			  }
		 }

		 elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }




        
        $total_enzyme_loci+=@enzyme_locs;
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
       
          my $front_loc_on_scaffold=$front_end_loc+1;
          my $behind_loc_on_scaffold=$behind_end_loc+1;
        
          push @lengths_of_fragments, $length_of_fragment;
                    
          # We treat the fragments in range in the same way of all fragments as above;
          if($length_of_fragment>=$range1){
            if($length_of_fragment<=$range2){
              $count_in_range++;
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
      }
      $seq_name_tmp=$seq_name;
      $seq_tmp="";
    }

    else{
      $seq_tmp=$seq_tmp.$line;
    }
    
    if( eof($ref_fh) ) {
      $scaffold_name=$seq_name_tmp;
      $scaffold_seq=$seq_tmp;
	  $scaffold_seq = uc $scaffold_seq;
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
      elsif($enzyme_seq=~/CCNNNNNNN/i){
			  my $loc_in_array;
			  $loc_in_array=index($seq_for_search,'CC');
			  unless($loc_in_array==-1){
				  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				  if(length $seq_after_CC >19){
					push @enzyme_locs, $loc_in_array;
				  }
			  }
			  while($loc_in_array!=-1){
				  $loc_in_array=index($seq_for_search, 'CC', $loc_in_array+1);
				  unless($loc_in_array==-1){
					  my $previous_loc=$enzyme_locs[$#enzyme_locs];
					  if(($loc_in_array- $previous_loc) > 13){
						  my $seq_after_CC=substr($seq_for_search,$loc_in_array);
				          if(length $seq_after_CC >19){
							  push @enzyme_locs, $loc_in_array;
						  }
					  }
				  }
			  }
		}


		elsif($string=~/(\w+)\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\]\[(\w+)\](\w+)/){
			  my $breast_bases=$1;
			  my $g1_bases=$2;
			  my $g2_bases=$3;
			  my $g3_bases=$4;
			  my $g4_bases=$5;
			  my $g5_bases=$6;
			  my $g6_bases=$7;
			  my $back_bases=$8;
			  my @g1_bases=split //, $g1_bases;
			  my @g2_bases=split //, $g2_bases;
			  my @g3_bases=split //, $g3_bases;
			  my @g4_bases=split //, $g4_bases;
			  my @g5_bases=split //, $g5_bases;
			  my @g6_bases=split //, $g6_bases;
			  for my $one_base(@g1_bases) {
				  for my $two_base(@g2_bases) {
					  for my $three_base(@g3_bases) {
						  for my $four_base (@g4_bases) {
							  for my $five_base( @g5_bases) {
								  for my $six_base(@g6_bases) {
									  $string=$breast_bases.$one_base.$two_base.$three_base.$four_base.$five_base.$six_base.$back_bases;
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
				  }
			  }
		  }
      
      
      
      $total_enzyme_loci+= @enzyme_locs; 
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
        my $front_loc_on_scaffold=$front_end_loc+1;
        my $behind_loc_on_scaffold=$behind_end_loc+1;
        push @lengths_of_fragments, $length_of_fragment;          
        # We treat the fragments in range in the same way of all fragments as above;
        if($length_of_fragment>=$range1){
          if($length_of_fragment<=$range2){
            $count_in_range++;
            push @lengths_fragments_in_range,$length_of_fragment;
            $length_of_fragments_in_range+=$length_of_fragment;
          }
        }
      }
      
    }
  }
  print "Done!\n";
  return ($total_enzyme_loci, \@lengths_of_fragments);
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
