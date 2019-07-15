#!/usr/bin/perl

use strict;
use warnings;
use lib "/wrk/hyvonen1/lib";
use MatrixReal; 

   # my ($input) = $ARGV[0];
    my $input = "";
    my $number_of_atoms;
    my $modes=1;
    my @zero_point_geo = ();
    my @eigen_modes = ();
    my @frequencies = ();
    my @ir_intensities = ();
    my @sorted_frequencies = ();
    my @translation_modes = ();
    my @unstable_modes = ();
    my $unstable_frequ = 0 ;
    my @distorting_modes = (); # For later: if user wants to put it by hand
    my @prefactor_mode = ();    # -> array lateron 
    my @distorted_geometry = () ;  
    my $project_geometry_file="";
    my $do_project = 0;
    my $xyz2molden = 0;
 

# TODO: lateron we might also change the overall norm for 
#       a global distortion


######## THESES VARIABLE SHOULD BE CAPTURED BY  A READING ROUTINE ############
    #TODO: romve this!
    my $chose_mode = 1 ;
##############################################################################

for (my $arg=0; $arg<@ARGV; $arg ++){
   if    ($ARGV[$arg] eq "-mode")   {push @distorting_modes, $ARGV[$arg+1]; print "Distortion along mode:       " , $ARGV[$arg+1],  " \n"}
   elsif ($ARGV[$arg] eq "-norm")   {push @prefactor_mode, $ARGV[$arg+1]; print "Norm of distorting mode:       " , $ARGV[$arg+1]," Angstroms \n"}
   elsif ($ARGV[$arg] eq "-xyz2molden"){$xyz2molden=1}
   elsif ($ARGV[$arg] eq "-project"){     $project_geometry_file=$ARGV[$arg+1]; $do_project = 1; 
				     print "Decomposition of " , $project_geometry_file," with respect to modes.\n"}
   elsif ($ARGV[$arg] =~ /.xyz/    ){$input=$ARGV[$arg]; print "# Input file: ",  $ARGV[$arg], "  \n"}
}


   # (1) Read vib.xyz 
   & read_vib ($input) ; 
   #& test_read_vib () ; 

  if($xyz2molden ==1) {
     & xyz2molden ();
  }

  elsif ($do_project==1) {
     print "\n...Projecting geometries \n\n";
     my @project_geo = () ;
     @project_geo = & read_geometry ($project_geometry_file) ;
 
     # Check geomtries for consistency and substract zero-point energy
     my @delta_geometry = & substract_geo (\@zero_point_geo, \@project_geo) ;

     # Project and write result
     & project(\@delta_geometry);
     
   }
 
   else {  # no projection -> distortion

   # (2) Check of arrays are consistent
   if($#distorting_modes != $#prefactor_mode ){
      print "\n 
      *** Number of chosen modes (",$#distorting_modes+1, ") 
          not equal to number of requested norms (", $#prefactor_mode+1 ,").\n\n" ;
      die;
   }

   # (3) Analyze frequencies --> enumerate with respect to mode number!
   & analyze_modes();

   # (4) )Decide along which modes we want to distort and print
   
   # No choice made -> distort all unstable ones
   if ($#distorting_modes == -1){   
      @distorting_modes=@unstable_modes; # All unstable ones
      # For now assign the prefactors to be 1
      for (my $k=0;$k < @distorting_modes ;$k++){
         $prefactor_mode[$k]=1.0;
      }

      if ($unstable_frequ != 0){
         & distort_and_print();
      }
   }

   else {    # Use input modes
         & distort_and_print();
   }
   } # No projection


############################################################################################################################
# SUBROUTIONES
############################################################################################################################
sub  test_read_vib{
   
   print "Number of atoms	", $number_of_atoms, "\n"; 
   print "\n";

   # TEST Zero point geometries 
   print  "TEST zero point energy \n";
   for (my $j=0;$j<@zero_point_geo;$j++){
   printf "%10.5f\t%10.5f\t%10.5f\t%s\n",
   $zero_point_geo[$j][0], $zero_point_geo[$j][1], $zero_point_geo[$j][2], $zero_point_geo[$j][3];
   }
   print "\n\n";

   # TEST Distorted Geometries
   for (my $k=0; $k < $modes-1 ; $k++){
      print "\n\n Mode  ",$k, "     Frequency:   ",$frequencies[$k][0],   "\n";
      for (my $j=0+($k*$number_of_atoms); $j<$number_of_atoms+($k*$number_of_atoms) ;$j++){
      printf "%10.5f\t%10.5f\t%10.5f\n",
      $eigen_modes[$j][0], $eigen_modes[$j][1], $eigen_modes[$j][2];
      }
   }

}


sub read_vib {

    my ($file_name) = @_; 
    my $frequency;
    my $ir_intensity;
    my $i;

    open  (INPUT, $file_name) or die "Need xyz file as input \n" ;
    my @line = <INPUT> ;

    # Get number of atoms
    my @split_line = split " ", $line[0] ;
    $number_of_atoms = $split_line[0];

    # Get zero-point goemetry
    for ($i=2; $i<=$number_of_atoms+1; $i++) {
       my @split_line = split " ", $line[$i] ;
       push @zero_point_geo , ([$split_line[1], $split_line[2], $split_line[3], $split_line[0]]);
    }


     for ($i=0; $i<=$#line; $i ++) {
        if ($line[$i] =~ /frequency/) {
        # (1) Get frequency
        my @split_line = split " ", $line[$i] ;
        for (my $k=0; $k<=$#split_line; $k ++) {
           if ($split_line[$k] =~ /at/){
           $frequency=$split_line[$k+1];
           $ir_intensity=$split_line[$k+6];
           push @frequencies ,  ([$frequency, $modes]);
           push @ir_intensities , $ir_intensity;
           #print "Mode:   ", $modes, "\t", "Frequency:   ", $frequency, "\n"; # TEST
           }
        }
        # (2) Get the geometry
        for (my $j=$i+1 ; $j<=$number_of_atoms+$i; $j++) {
           my @split_line = split " ", $line[$j] ;
           push @eigen_modes , ([$split_line[4], $split_line[5], $split_line[6]]);
        }

        $modes ++ ;
        }  
     }

    close (INPUT) ;
}


sub analyze_modes{
     #==================TRANSLATIONS=======================================

     # get absolute values in frequencies 
     for (my $i = 0; $i<@frequencies; $i++){
        $sorted_frequencies[$i][0] = abs ($frequencies[$i][0]);
        $sorted_frequencies[$i][1] =     ($frequencies[$i][1]);
     }

     # Sort  w.r.t. abs frequencies 
     @sorted_frequencies =  sort {$a->[0] <=> $b->[0]} @sorted_frequencies;

     # Print the 6 lowest lying frequencies
     print " Analyzing file: ", $input, "\n \n";
     print " TRANSLATIONS/ROTATIONS		# Mode	   ", "Frequency   (The six frequencies which are closest to zero)	 \n" ;
     print "----------------------------------------------------------------------------------------------\n" ;

     for (my $i=0; $i<6; $i++){
        my $transl_mode = $sorted_frequencies[$i][1] ;
        push @translation_modes, $transl_mode;
        printf "\t\t\t%10.0f\t%12.5f\n", $transl_mode, $frequencies[$transl_mode-1][0];
     }
    
     # Error message if first non-trivial frequency too close to translations:
     if (abs ($sorted_frequencies[5][0]-$sorted_frequencies[6][0]) < 10){
        print "*** Warning: there seem to be more frequencies close to zero than 6! Check your stuff ... \n";
     }

	#     for (my $i=0;$i<@sorted_frequencies;$i++){
	#     printf "%10.5f\t%2.0f\n",  $sorted_frequencies[$i][0], $sorted_frequencies[$i][1];
	#     }
      
    #==================NEGATIVE FREQUENCIES==================================

     print "\n UNSTABLE MODES			# Mode	   ", "Frequency	(Negative modes that are not unitary) \n" ;
     print "----------------------------------------------------------------------------------------------\n" ;
     for (my $i=0; $i<@frequencies; $i++){
        my $matches = grep {/$frequencies[$i][1]/} @translation_modes; # how if does it match?
        if (($frequencies[$i][0] < 0) &&  ($matches==0)){
           printf "\t\t\t%10.0f\t%12.5f\n", $frequencies[$i][1], $frequencies[$i][0];
           push @unstable_modes,$frequencies[$i][1] ; 
           $unstable_frequ ++;
        }
     }
     if ($unstable_frequ == 0){
        print "\n";
        print "There are no unstable frequencies -> structure seems to be in a local minimum! \n";
        print "\n \n";
        }

} # sub


sub distort_and_print{
   
     # counter
     my $atom = 0; 
     my $j    = 0; 
     my $k    = 0; 

     # Get distorted geometry
     for ($atom=0; $atom < $number_of_atoms; $atom ++){
        for ($j=0; $j<=3; $j++){
           $distorted_geometry[$atom][$j] =    $zero_point_geo [$atom][$j];
	}
     }
     # For every mode ...
     for ($k=0;$k<@distorting_modes;$k++){
        my $index = $distorting_modes[$k]-1; # which mode in the array eigen_modes?
        # ...distort all atoms along this mode
        for ($atom=0; $atom < $number_of_atoms; $atom ++){
           for ($j=0; $j<=2; $j++){
              $distorted_geometry[$atom][$j] +=    
              $prefactor_mode[$k]  * 
              $eigen_modes[$atom+($number_of_atoms*$index)][$j];
           }
        }
     }
     

     # Print distorted geometry
     print "\n";
     print "\n";
     print "# Distorted geometry obtained by distorting along the modes: \n";
     for (my $j=0;$j<@distorting_modes;$j++){
        print "# Mode number:  ", $distorting_modes[$j],
              "\t Norm:  ",$prefactor_mode[$j], " Angstroms \n" ;
     }
     print "\n";
     for ($atom=0; $atom < $number_of_atoms; $atom ++){
        printf "%s\t%12.6f\t%12.6f\t%12.6f\t%s\n",
        "atom", $distorted_geometry[$atom][0], $distorted_geometry[$atom][1], $distorted_geometry[$atom][2], $distorted_geometry[$atom][3]; 
     }
}


sub read_geometry {
   my $input = $_[0];
   my @geo   = (); 

   open  (INPUT, $input) or die "Need geometry file as input \n" ;
      while ($_ = <INPUT>){
      	  # last input fits regular expression ????
          if (/atom/){
      	  # create an array with all the entries in a line, starting index is 0
          my @line = split " ", $_ ; 
          if ( $line[0] eq "atom" ){
             # add coord triple to the end of array coords (higher dimensional array)
       	     push @geo, [($line[1], $line[2], $line[3],$line[4])] ;
             # $n_atoms++ ; -> for consistency check lateron
             }
         } 
   }
   close (INPUT) ;
   return (@geo) ;
}

sub substract_geo{
   my @geo_1 = @{$_[0]};
   my @geo_2 = @{$_[1]};
   my @diff_geo = () ;

   #------------------------------------------------------------------------------------------------------------------------------------------
   # Debug 
   #
   #         my  $j=@geo_2;
   #         print "\n ohne raute", $j,"mit raute", $#geo_2, " \n" ;   
   #
   #    		for (my $k=0;$k<@geo_1;$k++){
   #    		printf "%s\t%12.9f\t%12.9f\t%12.9f\t%s\n", "atom", $geo_1[$k][0],  $geo_1[$k][1],  $geo_1[$k][2],  $geo_1[$k][3]  ; 
   #             	}
   #
   #    		for (my $k=0;$k<@geo_2;$k++){
   #   		printf "%s\t%12.9f\t%12.9f\t%12.9f\t%s\n", "atom", $geo_2[$k][0],  $geo_2[$k][1],  $geo_2[$k][2],  $geo_2[$k][3] ; 
   #    		 }
   # Debug 
   #------------------------------------------------------------------------------------------------------------------------------------------


   # Have geometries the same dimension ?
   if ($#geo_1 ne  $#geo_2){
      print " *** Error in when substracting geometries: arrays have different dimensions:" , $#geo_1," vs " , $#geo_2,
            "  Please check that input structures have the same dimension !"; die ;}
  

   for (my $k=0;$k<@geo_1;$k++){

      # Check if species are the same
      if ($geo_1[$k][3] ne $geo_2[$k][3]){ 
         print " *** Species are not the same!", $k+1, " (uncommented) input atomic specification.\n\n"; die ; 
      }
      push @diff_geo, [($geo_2[$k][0]-$geo_1[$k][0],
                        $geo_2[$k][1]-$geo_1[$k][1],
                        $geo_2[$k][2]-$geo_1[$k][2],
                                      $geo_1[$k][3])] ;
   }


   return (@diff_geo);
}

sub project_v0 {
   my @geo = @{$_[0]};

  # For all modes ...
  for (my $m=0; $m<$modes-1; $m++){   
     my $index = ($number_of_atoms)*$m;
     # ... Project on every atom
     my $alpha = 0;   # scalar product for a given mode
     for (my $k=0;$k<@geo;$k++){

           $alpha += $geo[$k][0]*$eigen_modes[$index+$k][0] +
                     $geo[$k][1]*$eigen_modes[$index+$k][1] +
                     $geo[$k][2]*$eigen_modes[$index+$k][2] ;

     }
    printf "\t%s%5.0f%s\t%12.9f\n", " Projection on Mode ",$m+1, ":",$alpha; 
  }
  print "\n";
}

sub project {

   my @geo = @{$_[0]};
   my $geometry = new Math::MatrixReal(3*$number_of_atoms,1); 
      for(my $k = 0; $k <  $number_of_atoms ; $k++){
      $geometry->assign(3*$k+1,1,$geo[$k][0]); 
      $geometry->assign(3*$k+2,1,$geo[$k][1]); 
      $geometry->assign(3*$k+3,1,$geo[$k][2]);  
      }
  
#   print $geometry;
# DEBUG begin   compare with original array:
#   print "\n Compare to original array:\n"; 
#      for(my $k = 0; $k <  $number_of_atoms ; $k++){
#      printf "%12.9f\t%12.9f\t%12.9f\n", $geo[$k][0],  $geo[$k][1],  $geo[$k][2]
#      }
# DEBUG  end  

   my $eigen_mode_matrix = new Math::MatrixReal(3*$number_of_atoms, 3*$number_of_atoms); 
      for(my $l = 0; $l <  3*$number_of_atoms ; $l++){     # modes
         my $index = $l*$number_of_atoms;
         for(my $k = 0; $k <  $number_of_atoms ; $k++){  # components of each mode
            $eigen_mode_matrix->assign(3*$k+1,$l+1,$eigen_modes[$k+$index][0]);  # x-component
            $eigen_mode_matrix->assign(3*$k+2,$l+1,$eigen_modes[$k+$index][1]);  # y-component
            $eigen_mode_matrix->assign(3*$k+3,$l+1,$eigen_modes[$k+$index][2]);  # z-component
         }
      }
#   print $eigen_mode_matrix;
   
   # Invert eigen_mode_matrix
   my $inverse_modes = $eigen_mode_matrix->inverse();

   # Get prefactors of decomposition by multiplying M^-1 X, where X is geoemtry and M is the matrix of eigenmodes
   my $alpha = $inverse_modes->multiply($geometry);
#   print $alpha;    

   for(my $l = 0; $l <  3*$number_of_atoms ; $l++){     # modes
      printf "\t%s%5.0f%s\t%12.9f\n", " Projection on Mode ",$l+1, ":",$alpha->element($l+1,1); 
   }

}

sub  xyz2molden {
my $bohr=0.52918;
print "[Molden Format] \n";
print "[GEOMETRIES] XYZ \n";
print "   ",$number_of_atoms, "\n\n";

for (my $j=0;$j<@zero_point_geo;$j++){
   printf "%s\t%10.5f\t%10.5f\t%10.5f\n",
   $zero_point_geo[$j][3], $zero_point_geo[$j][0], $zero_point_geo[$j][1], $zero_point_geo[$j][2];
}

print  "[FREQ]\n";
for (my $k=0; $k<=$#frequencies; $k ++) {
           print $frequencies[$k][0] , "\n" ;
}

print  "[INT]\n";
for(my $k=0; $k<=$#frequencies; $k++){
      print $ir_intensities[$k], "\n";
}

print  "[FR-COORD]\n";
for (my $j=0;$j<@zero_point_geo;$j++){
   printf "%s\t%10.5f\t%10.5f\t%10.5f\n",
   $zero_point_geo[$j][3], $zero_point_geo[$j][0]/$bohr, $zero_point_geo[$j][1]/$bohr, $zero_point_geo[$j][2]/$bohr;
}


print  "[FR-NORM-COORD]\n";
for (my $m=0; $m<$modes-1; $m++){   
  print "vibration\t",$m+1, "\n";
  my $index = ($number_of_atoms)*$m;
  for (my $k=0;$k<$number_of_atoms;$k++){
         printf "\t%12.6f\t%12.6f\t%12.6f\n",   $eigen_modes[$index+$k][0]/$bohr, $eigen_modes[$index+$k][1]/$bohr, $eigen_modes[$index+$k][2]/$bohr;

     }
  }
} 


