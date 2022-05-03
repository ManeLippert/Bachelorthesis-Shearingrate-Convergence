#!/bin/sh
#  -*-Perl-*-  (for Emacs)    vim:set filetype=perl:  (for vim)
#======================================================================#
# Run the right perl version:
if [ -x /usr/local/bin/perl ]; then
  perl=/usr/local/bin/perl
elif [ -x /usr/bin/perl ]; then
  perl=/usr/bin/perl
else
  perl=`which perl| sed 's/.*aliased to *//'`
fi  
#   
exec $perl -x -S $0 "$@";     # -x: start from the following line
#======================================================================#
#! /Good_Path/perl -w 
# line 17
#
#  Update the GKW input files to work with a the code when some switches have
#  been moved from the CONTROL namelist to DIAGNOSTIC. Takes the local file
#  `input.dat' and outputs the update as `input.dat.new'.
#  

use strict;

my $control_nml  = '^\s*\&control\s*';
my $diagnostic_nml  = '^\s*\&diagnostic\s*';
my $nml_end = '^\s*\/\s*';
my @lines_to_move;
my @moved_vars = qw(
       output3d
       lwrite_output1
       lphi_diagnostics
       screen_output
       lfinal_output
       lcalc_fluxes
       lcalc_freq
   );
my @deleted_vars = qw(
       ldiagnostic_namelist
   );
my $input_file = 'input.dat';
my $new_file_suffix = '.new';


# first find the lines in control to move
find_lines_to_move($input_file);

# create a new file with the removed/moved lines
create_new_file($input_file);



#----------------------------------------------------------------------------
# find the lines of &CONTROL to move to &DIAGNOSTIC
# 
sub find_lines_to_move{
    my $file = shift;
    my $check_line = 0;
	
    open(INPUT,"<$file") || die "can not open file $file\n";
	
	# read through the whole file a line at a time
    while( my $line = <INPUT> ) {
		
        chop($line);
		# check if we are at the control nml
        $check_line = 1 if ($line =~ m/$control_nml/i );
	    $check_line = 0 if ($line =~ m/$nml_end/i);
		
		# make a list of the items to move from control
	    if ($check_line) {
            foreach my $item (@moved_vars) {
	            if ($line =~ m/^\s*$item/i ) {
		            push @lines_to_move,$line;
				    #   print "WILL MOVE:\n";
				    #   print "   $line\n";
		        }
	        }
	    }
    }
    close(INPUT);
}

#----------------------------------------------------------------------------
# write the modified namelist to a new file
# 
sub create_new_file{
	
    my $oldfile = shift;
    my $newfile = $oldfile.$new_file_suffix;

    # check if the diagnostic namelist exists already
    my $have_diagnostic = check_for_diagnostic($oldfile);
    my $check_line =0;
    my $write_line =1;
    
	open(INPUT,"<$oldfile") || die "can not open file $oldfile\n";
    open(OUTPUT,">$newfile") || die "can not open file $newfile for output\n";

	while( my $line = <INPUT> ) {
    
		my $line_out = $line;
        chop($line);
        $check_line = 1 if ($line =~ m/$control_nml/i );
  	    $check_line = 0 if ($line =~ m/$nml_end/i);
  	    $write_line = 1;
  	    
		# do not write moved or deleted lines back into the CONTROL namelist
		if ($check_line) {
            foreach my $var (@moved_vars) {
  	            $write_line = 0 if ($line =~ m/$var/i);
            }
  	        foreach my $var (@deleted_vars) {
  	            $write_line = 0 if ($line =~ m/$var/i);
  	        }
  	    }
        
		# write the moved lines at the start of the existing diagnostic
		# namelist if it exists
		if ($line =~ m/$diagnostic_nml/i ) {
           
  	        print OUTPUT "$line_out";
			
  	        foreach my $new_line (@lines_to_move) {
  	            print OUTPUT "$new_line\n";
  	        }
			
  	    } else {
	    # just output the original line if we are not in the diagnostic part
  	        print OUTPUT "$line_out" if ($write_line);
  	    }
		
    }

    # add the diagnostic namelist at the end if it did not exist already
    if ( ! $have_diagnostic) {
		
        print OUTPUT "\&DIAGNOSTIC\n";
		
        foreach my $line (@lines_to_move) {
            print OUTPUT "$line\n";
        }
		
        print OUTPUT "\/\n";
    }
	
    close(INPUT);
    close(OUTPUT);
}

#----------------------------------------------------------------------------
# check if the diagnostic namelist exists anywhere in the input file
# 
sub check_for_diagnostic {

	my $file = shift;
    my $diagnostic_present = 0;

	open(INPUT,"<$file") || die "can not open file $file\n";

	while( my $line = <INPUT> ) {
        chop($line);
	    $diagnostic_present = 1 if ($line =~ m/$diagnostic_nml/i );
    }
	
    close(INPUT);
    
	return $diagnostic_present;
}

#----------------------------------------------------------------------------
