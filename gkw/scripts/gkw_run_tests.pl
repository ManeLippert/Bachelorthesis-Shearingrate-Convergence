#!/bin/sh
#  -*-Perl-*-  (for Emacs)    vim:set filetype=perl:  (for vim)
#============================================================================#
# Run the right perl version:
if [ -x /usr/local/bin/perl ]; then
    perl=/usr/local/bin/perl
elif [ -x /usr/bin/perl ]; then
    perl=/usr/bin/perl
else
    perl=`which perl | sed 's/.*aliased to *//'`
fi

exec $perl -x -S $0 "$@"     # -x: start from the following line
#============================================================================#
#! /Good_Path/perl -w

use strict;
use File::Compare;
use Getopt::Long;
use POSIX;
use File::Copy;
use File::Path qw(mkpath);
use File::Basename;
use File::Temp qw(tempfile);

# renice itself
my $nice_output = `renice -n +15 -p $$`;

my $help = 0;
my $barrier =
    "+----+------------------------------------------------------------\n";
my $FORCE_NONE = 3;
my $FORCE_SERIAL = 1;
my $FORCE_PARALLEL = 2;
my $force_serialparallel = 0;
my $force_serialruns_w_mpirun = 0;
my $compare_to_serial = 0;
my $just_list = 0;

# the command to run mpi executables
my $mpiruncmd      = "mpirun";
# get the right mpirun command (or mpiexec) to use later
get_mpiruncmd();

# the argument to make the tests run with the newest executable, but
# ignore newer compilations while the tests run.
my $fix_executable = "fix";

# default: the input filename for GKW
my $input_filename = "input.dat";
my $cwd = getcwd();
chomp($cwd);

my @testdirs;
my %cap;

# information about the tests
my @failed_tests_running;    # somewhere to put the output
my @failed_tests_output;
my $nfailed_running  = 0;
my $nfailed_output   = 0;
my $ntests_completed = 0;
my $ntests_failed    = 0;

my (%opts);                  # variables written by GetOption
(my $cmdname = $0) =~ s{.*/}{};

# How do we determine if the code ran succesfully?
my $successful_run_line =
  "Run successfully completed";    # find this in the output

# get input from various switches
Getopt::Long::config("bundling");
GetOptions(
    \%opts,
    qw( -h     --help
        -l     --list
        -s     --serial
        -p     --parallel
        -m     --mpirun
               --compare-to-serial
               --clean
        -d=s   --code_home=s
        -e=s   --executable=s
        -i=s   --inputfile=s
        -I=s   --ignore=s
               --diff-command=s
               --gkw-check-cleanup=s
               --gkw-generate-parallelisations=s
      )
    ) or $help = 1;

# get and check the GKW_HOME environment variable
my $code_home = ($opts{'d'} || $opts{'code_home'} || "$ENV{GKW_HOME}" || "none");
if ("$code_home" eq "none")
{
    die
      "You need to set the environment variable GKW_HOME to where the GKW code is located.\n";
}


# the standard directory to look for testcases
my $default_testdir = "$code_home/tests/standard";
my @testcases = ("$default_testdir");

my $filename_regexp_not_to_diff_default = "^perform.*|^perfloop.*|^resource_usage";
my $filename_regexp_not_to_diff = ($opts{'I'} || $opts{'ignore'} || "$filename_regexp_not_to_diff_default");
print "Files matching the pattern $filename_regexp_not_to_diff are ignored for diffs.\n";

my $STANDARD_DIFF = "diff";
my $diff_command = ($opts{'diff-command'} || $STANDARD_DIFF);

# where to look for executables (relative to the GKW_HOME directory for now)
my $execdir = "$code_home/run";


# ARGV holds now the command line arguments which were not parsed.
# Any remaining command line arguments should be
#  * input files or
#  * directories
if (@ARGV) {
    @ARGV = map { s{^(?!/)}{$cwd/}; $_ } @ARGV;
    for my $d (@ARGV) {
        #print "$d\n";
    }
    @testcases = @ARGV ;
    chomp(@testcases);
}

if ($opts{'h'} || $opts{'help'} || $help == 1) { usage(); }
if ($opts{'s'} || $opts{'serial'})
{
    $force_serialparallel = $force_serialparallel | $FORCE_SERIAL;
}
if ($opts{'p'} || $opts{'parallel'})
{
    $force_serialparallel = $force_serialparallel | $FORCE_PARALLEL;
}
if($force_serialparallel == 0)
{
    $force_serialparallel = $FORCE_NONE;
}

if ($opts{'m'} || $opts{'mpirun'})
{
    $force_serialruns_w_mpirun = 1;
}

if ($opts{'compare-to-serial'})
{
    $compare_to_serial = 1;
}

my $gkw_check_cleanup = ($opts{'gkw-check-cleanup'} || "none");

if ($opts{'l'} || $opts{'list'})
{
    $just_list = 1;
}

my $gkw_generate_parallelisations = ($opts{'gkw-generate-parallelisations'} || "none");

# if the executable is specified as a command line argument, read it.
my $exe = ($opts{'e'} || $opts{'executable'} || "none");
if (("$exe" eq "none" || "$exe" eq "$fix_executable") && $just_list == 0)
{
    # default: find all the executables in the GKW run directory
    # and sort
    my @executables = find_gkw_testing_executables("$execdir");
    if(@executables > 0){
        if ("$exe" eq "$fix_executable")
        {
            my $tmphandle;
            my $tmpfilename;
            ($tmphandle, $tmpfilename) = tempfile();
            copy ($executables[0],"$tmpfilename");
            chmod(0744, $tmpfilename);
            print "Using a fixed executable for runs.\n";
            $exe = $tmpfilename;
        }
        else
        {
            $exe = $executables[0];
        }
    }else{
        $exe = "_none_";
    }
}
if ($exe eq "_none_")
{
    print "\nERROR: NO executables found!\n";
    print "(maybe you need to compile the code)\n";
    die;
}
print "Tests are performed with the executable \n\n";
print " * $exe\n\n";

# if the executable is specified as a command line argument, read it.
$input_filename = ($opts{'i'} || $opts{'inputfile'} || $input_filename);

# Try to parse the capabilities of the executable.
# At the moment this is not used.
#get_executable_capabilities($default_testdir, $exe);


print "Check if h5diff command is available...";
my $h5diff_avail = (system("h5diff --version") == 0);
if(!$h5diff_avail)
{
    print "-> differences in HDF5 files will not be checked.\n";
}else{
    print "-> differences in HDF5 files can be checked.\n";
}

print "Check if /usr/bin/time command is available...";
# FIXME ? can you make this quiet:
my $time_avail = (system("/usr/bin/time --version") == 0);
if(!$time_avail)
{
    print "-> resource usage will not be measured.\n";
}else{
    print "-> resource usage will be measured in a very simple way.\n";
}

# do the tests
perform_tests();

die "\n\n" . " DONE " . "\n";

#=============================================================================

sub usage
{
    print STDERR << "EOF";

Usage:  $cmdname [options] [dir1 [dir2 [..]]]

Run GKW testcases or show information about them.

If one or more directories are given, they will be searched
recursively for testcases.
If no directories are given, then all testcases below the folder
    $default_testdir
will be run.

In addition, it is checked if the source code, consisting of all
Fortran files (filenames match *[Ff]90) below the code_home,
contains non-ascii characters or tabulators. If it does, a message is
displayed as testcase no. -1

For now, the script looks for any executable files gkw*.x in the
directory
    $execdir
The default is to test the executable with the most recent
modification time.  You can set a previously compiled executable to be
the one to test by simply doing "touch <your-old-gkw-executable> ".

This script evaluates the env variable \$GKW_MPIRUNCMD to determine the
correct run command for parallel runs, e.g. mpirun (with full path if
necessary).

Options:
  -h,  --help              Display this usage message
  -l,  --list              List the test input files, then exit
  -d,  --code_home=<dir>   Use <dir>, instead of the directory specified
                           in the environment variable \$GKW_HOME
  -e,  --executable=<file> Use <file> as the executable, instead of looking for
                           the most recent GKW executable
  -e=fix, --executable=fix Use a temporary copy of the most recent GKW executable.
                           This allows to recompile the code without affecting
                           testcases which have not yet started.
  -i,  --inputfile=<name>  Use <name> as name for input
                           files, instead of 'input.dat'
  -I,  --ignore=<regexp>   Ignore filenames matching the pattern <regexp> when
                           it comes to finding diffs. The inputfile and .md5
                           files are always ignored regardless of this pattern.
                           Default value is $filename_regexp_not_to_diff_default .
  -s,  --serial            Consider only the serial tests
  -p,  --parallel          Consider only the parallel tests
  -m,  --mpirun            Do the serial runs with $mpiruncmd, too.
       --diff-command=<command>
                           An external command (probably other than standard unix
                           diff) to be used for computing the diff.
       --compare-to-serial
                           Do not compare data generated by parallel tests with
                           files in the respective reference folder, but with
                           files generated by the serial test. Use this to test
                           the accordance of certain files without adding them
                           to the repository.
       --clean             Clean up the test directories, then exit

Additionally there are some GKW specific options:

       --gkw-check-cleanup=<script>
                           Immediately after every run, the given script is 
                           invoked and the remaining files are compared to 
                           the testcase basefolder, to see if the cleanup 
                           script misses something.
       --gkw-generate-parallelisations=<absolutepath>
                           GKW-specific: Considering the parallelisation
                           scheme of the respective testcases, generate
                           testcases for certain sub-schemes (e.g. only in
                           species, only in velspace, etc., see source code
                           of $cmdname to see which ones)
                           and put them into the given folder. To run them,
                           a second invocation of this script is necessary.

Mandatory arguments to long options are also mandatory
for any corresponding short options.

Keywords:
For special purposes, keywords contained in the input file are evaluated.

gkw_run_tests:do_not_cleanup_before_run : If this is contained somewhere
  in the input file, the files in the run folders are not deleted before
  copying the new files and running the executable. This can be used in
  testcases to test restarting.
       
Examples:
  # run all tests in the default folder:
  $cmdname
  # run all tests in the given folders:
  $cmdname extra/
  # run tests using the executables and
  # test files found in /tmp/gkw:
  $cmdname -d /tmp/gkw

  # run all parallel nonlinear testcases of the
  # standard test set:
  $cmdname --parallel \$(grep -ril 'non_linear[^fF]*[tT]' \$GKW_HOME/tests/standard | sed -e 's\|[^/]*\$\|\|' | uniq)
  # run all serial testcases which produce the 'kxspec' file
  $cmdname --serial \$(find \$GKW_HOME/tests/standard -name 'kxspec*' | sed -e 's\|[0-9]*/[^/]*\$\|\|' | uniq)
    
  # list all serial testcases below the current directory:
  $cmdname --serial -l .

  # With GKW's perform.dat output in the reference folder, one can
  # check if execution performance of testcases has changed, using a
  # tool like numdiff. The ignore option must be specified,
  # because performance output files are catched by the default
  # ignore-pattern:
  $cmdname --serial --diff-command="numdiff -V -r 0.05 -a 1.0 --strict" --ignore='^\$' \$GKW_HOME/tests/more/performance

  # run all tests in the default folder. Compare the data of the serial runs to
  # the reference files, and the data of the parallel runs to the serial data.
  $cmdname --compare-to-serial

  # run serial tests and check if the cleanup script will remove everything
  # that is produced
  $cmdname --gkw-check-cleanup=gkw_clean_run --serial

  # generate GKW testcases with some other parallelisation sub-schemes,
  # derived from the standard testsuite, run them, and compare the
  # data to the serial runs
  $cmdname --serial
  $cmdname --gkw-generate-parallelisations=\$GKW_HOME/tests/more
  $cmdname --parallel --compare-to-serial \$GKW_HOME/tests/more

  # run the tests parallelised through valgrind. This requires valgrind to be
  # compiled with MPI. Evtl. replace \$SHAREDLIBS with the place to
  # your valgrind directory, adjust the platform of the library and
  # pick the correct gkw executable.
  GKW_MPIRUNCMD="LD_PRELOAD=\$SHAREDLIBS/valgrind-current/lib64/valgrind/libmpiwrap-amd64-linux.so mpirun" gkw_run_tests --executable="\$SHAREDLIBS/valgrind-current/bin/valgrind --leak-check=full --show-leak-kinds=all \$GKW_HOME/run/gkw_btppx-gnu-DP_0.3-b2-139-g9368ca7-dirty.x"

  # test Arne's RMHD solver
  $cmdname -e ~/linar/linar -i 'linar_input.dat' ~/linar/test

Testcases:
  A testcase is a folder which contains
     * the file $input_filename
     * all further necessary files, like restart files. These files should not
       have the suffixes .txt or .info .
     * a subfolder reference/ with all output files the test run is to compare
       against. When you construct a new testcase, you first generate the
       reference files by running with an old (trusted) excutable, then copy
       them to the reference folder.
       Note that binary files cause difficulties in comparison, you should
       consider testing them just via the md5 checksum.
     * From a trusted run, one can generate md5 checksums with the command
         md5sum onefiletocheck.dat otherfiletocheck > checksums.md5
       If this file containing the checksums is put to
       reference/checksums.md5, then the script will compare the md5 sums
       of the testruns with them.
     * one or several (empty) run folders with numbers as names. This causes
       test runs with the respective number of MPI processes.
       Example: the subfolders 1/ and 8/ exist. Thus, this testcase is once run
       serially and once with 8 MPI processes.
  Note that files with the suffixes *.txt or *.info like COMMENTS.TXT are not
  copied to the run folders.

EOF
    exit;
}

#------------------------------------------------------------------------------
sub perform_tests
{

    print "\n\n";

    if ($just_list) {
        print "Numbers in brackets denote the number of MPI processes.\n\n";
    }

    print "$barrier";
    printf "| %2d | Check if source code is ASCII and does not contain tabulators:\n", -1;
    print "+----+\n";
    if(check_code_is_ascii_notabs()) {
        print "\n OK - only ASCII letters, no tabs.\n";
    } else {
        print "\n BAD - there are non-ASCII letters or tabs!\n";
    }

    # loop over the given testcase directories.
    foreach my $file (@testcases) {
        if ( -d $file) {
            #print "I will now look for testcases in $file\n";
            
            # $file might be either a folder containing a testcase
            # or a folder which contains subfolder containing several testcases
            # (hence a testcase set)

            visit_folder($file);
            
        } elsif ( -f $file ) {
            print "$barrier";
            print "$file is a file, not a folder. Give me folders!\n\n";
        }
    }

    print_summary();

    # exit with nonzero exit code if there were failed runs or tests
    exit ($ntests_failed > 0);
}

sub visit_folder
{
    my $tstdir = shift;
    #print "** VISIT $tstdir\n";

    # check if it contains an input file
    if ( is_a_testcase($tstdir) ) {
        # This means the folder $tstdir contains one testcase.

        if ($just_list) {
            # We only list it.
            my @nprocs_list = get_procs($tstdir);
            my $nprocs = join(',', @nprocs_list);
            chomp($nprocs);
            print " $tstdir  \($nprocs\)\n";
        } elsif ($gkw_generate_parallelisations ne "none") {
            print "$barrier";
            printf "| %2d | In folder '%s'\n", $ntests_completed, $tstdir;
            print "+----+\n";
            
            chdir($tstdir);

            #generate testcases for particular schemes and put
            #them into the given folder
            my @par_scheme;
            # only one direction
            @par_scheme = ("n_procs_x");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_s");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_sp");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_mu");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_vpar");
            generate_other_parallelisation($tstdir, @par_scheme);

            # only position space
            @par_scheme = ("n_procs_x", "n_procs_s");
            generate_other_parallelisation($tstdir, @par_scheme);
            # only velspace
            @par_scheme = ("n_procs_mu", "n_procs_vpar");
            generate_other_parallelisation($tstdir, @par_scheme);

            # species and one other direction
            @par_scheme = ("n_procs_sp", "n_procs_x");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_sp", "n_procs_s");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_sp", "n_procs_mu");
            generate_other_parallelisation($tstdir, @par_scheme);

            # all directions but one
            @par_scheme = (              "n_procs_mu", "n_procs_vpar", "n_procs_s", "n_procs_x");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_sp",               "n_procs_vpar", "n_procs_s", "n_procs_x");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_sp", "n_procs_mu",                 "n_procs_s", "n_procs_x");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_sp", "n_procs_mu", "n_procs_vpar",              "n_procs_x");
            generate_other_parallelisation($tstdir, @par_scheme);
            @par_scheme = ("n_procs_sp", "n_procs_mu", "n_procs_vpar", "n_procs_s"             );
            generate_other_parallelisation($tstdir, @par_scheme);

            #(parallelisation in certain directions only)
            # perform_one_run($dir, $sfac)
            chdir("..");
        } else {
            # Execute the requested runs.
            perform_one_test($tstdir);
        }
        
    }
    # check if the subfolders contain testcases
    opendir(my $DIR, $tstdir)
        || die "get_test_directories: can't open $tstdir to find tests\n";
    # get a list of subfolders, excluding such which begin with a dot or a number
    my @subfolders = grep { /^[^\.0-9]/ && -d "$tstdir/$_" } readdir($DIR);
    closedir $DIR;

    # sort lexically
    @subfolders = sort @subfolders;

    # look into every subfolder
    foreach my $subfolder (@subfolders) {
        #print " I am visiting $tstdir/$subfolder\n";
        if ( -d "$tstdir/$subfolder" ) {
            # this might be either a folder containing a testcase
            # or a folder which contains subfolder containing several testcases
            # (hence a testcase set)
            visit_folder("$tstdir/$subfolder");
        }
    }
}

            
sub smallest_factor
{
    my $elems_tot = shift;
    my $greater_than = shift;

    # look for the largest number which divides $elems_tot but is
    # unequal $elems_tot itself
    for (my $elems_pp = $elems_tot-1; $elems_pp >= $greater_than; $elems_pp-- )
    {
        if ($elems_tot % $elems_pp == 0)
        {
            #if each process works on elems_pp points, then there is
            #this number of processes:
            return $elems_tot/$elems_pp;
        }
    }
    return $elems_tot/$greater_than;
}

sub are_all_parallel
{
    my $par_schemeRef = shift;
    my $pRef = shift;

    my $ret = 1;
    foreach my $dim (@{$par_schemeRef})
    {
        $ret = ($ret && ($$pRef{$dim} > 1));
    }

    return $ret;
}

sub generate_other_parallelisation
{
    my $tstdir = shift;
    my (@par_scheme) = @_;

    # minimum number of gridpoints per process in the respective
    # dimensions
    my %grid_min = ("n_s_grid", 2,
                    "nx", 2,
                    "n_mu_grid", 2,
                    "number_of_species", 1,
                    "n_vpar_grid", 2);

    my %corresp = ("n_procs_s", "n_s_grid",
                   "n_procs_x", "nx",
                   "n_procs_mu", "n_mu_grid",
                   "n_procs_sp", "number_of_species",
                   "n_procs_vpar", "n_vpar_grid");
    my %p;
    
    # read the grid and parallelisation scheme parameters from the
    # input file
    my @paramlist=("n_procs_s", "n_procs_x", "n_procs_vpar", "n_procs_mu", "n_procs_sp", "nx", "n_s_grid", "n_mu_grid", "n_vpar_grid", "number_of_species");
    foreach my $paramname (@paramlist)
    {
        my $val = get_param_from_gkw_inputfile($paramname, 1);
        $p{$paramname} = $val;
    }
    
    if(are_all_parallel(\@par_scheme,\%p))
    {
        print ("\n");
        print("-> generate testcase with parallelisation only in: @par_scheme..\n");

        #determine the new parallelisation
        my $nprocs = 1;
        my $sed_args = "";
        my %prepl;
        my $repl_string="";
        foreach my $dim (@par_scheme)
        {
            my $corresp_grid = $corresp{$dim};
            $prepl{$dim} = smallest_factor($p{$corresp_grid},$grid_min{$corresp_grid});
            $nprocs *= $prepl{$dim};
            $sed_args .= " -e 's/".$dim."[^p].*/".$dim." = ".$prepl{$dim}."/i'";
            $repl_string .= "_".$dim."_".$prepl{$dim};
        }

        if($nprocs > 1)
        {
            #first: generate a plain serial testcase:

            # copy necessary files of this testcase into the run folder.
            opendir(my $DIR, '.');
            # version 1: just copy all files
            #my @files_to_copy = grep { -f "$dir/$_" } readdir($DIR);
            # version 2: copy all files, except those whose name
            # contains the given strings. This allows to put a COMMENTS.txt next to
            # the file necessary for the testcase.
            my @files_to_copy = grep { !/COMMENTS.txt|info/ && -f "./$_" } readdir($DIR);
            closedir $DIR;
            my $complete_path = "$gkw_generate_parallelisations/".basename($tstdir)."/$repl_string/";
            mkpath($complete_path);
            mkpath("$complete_path/$nprocs");
            foreach my $file_to_copy (@files_to_copy)
            {
                copy ("./$file_to_copy","$complete_path/$file_to_copy");
            }
            # create a symlink to the references
            unlink getcwd()."/reference";
            unlink getcwd()."/1";
            my $symlink_exists;
            $symlink_exists = eval { symlink(getcwd()."/reference","$complete_path/reference"); 1 };
            $symlink_exists = eval { symlink(getcwd()."/1","$complete_path/1"); 1 };
            #make it a serial testcase:
            my $sed_out;
            $sed_out = `sed -i -e 's/n_procs_mu.*/n_procs_mu = 1/i' -e 's/n_procs_vpar.*/n_procs_vpar = 1/i' -e 's/n_procs_s[^p].*/n_procs_s = 1/i' -e 's/n_procs_sp.*/n_procs_sp = 1/i' -e 's/n_procs_x.*/n_procs_x = 1/i' $complete_path/$input_filename`;
            print($sed_out);
            

            #second: explicitely set certain n_procs_ parameters
            $sed_out = `sed -i $sed_args $complete_path/$input_filename`;
            print($sed_out);

            print("$complete_path\n");
        }

    }
}


sub is_a_testcase
{
    my $tstdir = shift;
    return ( -e "$tstdir/$input_filename" && -d "$tstdir/reference" );
}

#--------------------------------------------------------------------------------

sub perform_one_test
{
    my $dir = shift;
    
    chdir($dir) || die "cannot cd to $dir \n";

    print "$barrier";
    printf "| %2d | In folder %s:\n", $ntests_completed, $dir;
    print "+----+\n";
    my @procs           = get_procs($dir);

    # perform one simulation run for each requested parallelisation
    # scheme.
    RUN: foreach my $nprocs (@procs)
    {
        perform_one_run($dir, $nprocs)
    }

    # Successful or not, this test is completed
    $ntests_completed++;

}

sub get_param_from_gkw_inputfile
{
    my $param = shift;
    my $defaultval = shift;

    my $retval = $defaultval;
    
    open(INPUT,"<$input_filename") || die "can not open file $input_filename\n";
    while( my $line = <INPUT> )
    {
        chop($line);
        $line = lc $line;
        if ($line =~ m/^\s*$param\s*=\s*[^!\,\s]*/i )
        {
            $line =~ s/^\s*$param\s*=\s*([^!\,\s]*).*/$1/g;
            $retval = $line;
            #print("\$$param = $retval;\n");
            last;
        }
    }
    close(INPUT);
    return $retval;
}

sub perform_one_run
{
    my $dir = shift;
    my $nprocs = shift;
        
    # if the input file contains this snippet somewhere then do not
    # delete existing files before running. This can be used for
    # special testcases, e.g. to test restarting.
    my $grep_for_keyword = `grep 'gkw_run_tests:do_not_cleanup_before_run' $input_filename`;
    if ($opts{'clean'} || $? ne 0)
    {
        # delete all preexisting files
        clean_dir("$dir/$nprocs");
    }
    if ($opts{'clean'})
    {
        print "cleaned up $dir/$nprocs\n";
        return;
    }
        
    my $code_output = "$successful_run_line";

    # go into the run folder
    chdir("$dir/$nprocs");

    # copy necessary files of this testcase into the run folder.
    opendir(my $DIR, $dir);
    # version 1: just copy all files
    #my @files_to_copy = grep { -f "$dir/$_" } readdir($DIR);
    # version 2: copy all files, except those whose name
    # contains the given strings. This allows to put a COMMENTS.txt next to
    # the file necessary for the testcase.
    my @files_to_copy = grep { !/txt|info/ && -f "$dir/$_" } readdir($DIR);
    closedir $DIR;
    foreach my $file_to_copy (@files_to_copy)
    {
        copy ("$dir/$file_to_copy","$dir/$nprocs/$file_to_copy");
    }

    ##########################################################
    # Run the simulation
    print "running $dir with $nprocs processes\n";
    #die 0;
    $code_output = run_code($nprocs);
    ##########################################################

    chdir("$dir");
    
    if (check_run($code_output))
    {
        print "The run terminated successfully.\n";

        check_forbidden_files($nprocs);
        
        my $refdir = get_reference_folder($nprocs);
        my @reference_files = get_reference_files($refdir);

        print "The avail. reference files are:\n";
        print "@reference_files\n";

        if (perform_checks($nprocs,$refdir,@reference_files))
        {
            print "Test $dir with $nprocs processes failed.\n\n";
            $failed_tests_output[$nfailed_output] = "$dir $nprocs";
            $nfailed_output = $nfailed_output + 1;
            $ntests_failed  = $ntests_failed  + 1;
        }
        else
        {
            print "Test $dir with $nprocs processes successful.\n\n";
        }
        
        # Evtl run the cleanup script and check what's left
        if($gkw_check_cleanup ne "none")
        {
            my $cleanup_output = `(cd $nprocs; $gkw_check_cleanup)`;
            my $diff_output = `(diff --exclude=reference --exclude=resource_usage --exclude="[0-9]*" --exclude=out --exclude=.keepfolder --exclude=input.out --brief ./ $nprocs/ )`;
            if($? ne 0)
            {
                print "Cleanup script does not clean properly:\n";
                print "$diff_output";
                print "\n";
            }
            
        }

    }
    else
    {
        print "#  RUN $dir with $nprocs processes FAILED! (check output below)\n";
        print "\n\n$code_output\n\n";
        print "#  RUN $dir with $nprocs processes FAILED! (check the output above)\n";
        
        $failed_tests_running[$nfailed_running] = "$dir with $nprocs processes";
        $nfailed_running = $nfailed_running + 1;
        $ntests_failed   = $ntests_failed   + 1;
    }
}

#------------------------------------------------------------------------------
# clean up the local run directory
# (uses gkw_clean_run for now)
#

sub clean_dir
{
    my $dir = shift;
    my $cleanstat = 1;
    my $cleanout  = `rm -f $dir/*`;
    #or like that?
    #unlink glob "$dir/*";
    return
}

#------------------------------------------------------------------------------
# run the code and pipe both stdout and stderr to the file called 'out'.
#
sub run_code
{

    my $procs   = shift;
    my $output;
    my $runstat = 0;
    my $runcmd;

    # these depend much on the conditions when running, so leave them out:
    #elapsed real time: %E
    #time in kernel mode: %S
    #CPU: %P
    #filesystem inputs:%I
    #filesystem outputs: %O
    #major pagefaults: %F
    #minor pagefaults: %R
    #swaps: %W
    my $memory_and_timing;
    if($time_avail)
    {
        $memory_and_timing='/usr/bin/time -f "time in user mode: %U\nmemory(kB,avg.text): %X\nmemory(kB,avg.data): %D\nmemory(kB,peak max): %M\nmemory(kB,avg total): %K" --output=resource_usage ';
    } else {
        $memory_and_timing='';
    }
        

    # set and print run command
    if ($procs == 1 && $force_serialruns_w_mpirun == 0)
    {
        $runcmd = "$exe 2>&1";
    }
    elsif ($procs == 1 && $force_serialruns_w_mpirun == 1)
    {
        $runcmd = "$mpiruncmd -np 1 $exe 2>&1";
    }
    else
    {
        $runcmd = "$mpiruncmd -np $procs $exe 2>&1";
    }

    print "... with $runcmd\n";

    # run the test, putting the std output into the variable $output.
    
    # One would like to schedule the run with low priority to allow
    # backups processes and working. However, using 'nice' inhibits
    # misusing the GKW_MPIRUNCMD to set env variables (e.g. for valgrind)
    # Therefore the whole script renices itself in the very beginning.
    $output = `($memory_and_timing $runcmd | tee out)`;

    return ($output);

}

#------------------------------------------------------------------------------
# check the run exited successfully
#
sub check_run
{
    my $output  = shift;
    return ($output =~ /$successful_run_line/);
}

#------------------------------------------------------------------------------
# check if certain files exist, which should not
#
sub check_forbidden_files
{
    my $nprocs = shift;

    opendir(DIR, "$nprocs");
    my @files = readdir(DIR);
    closedir(DIR);

    my $failed = 0;
    
    foreach my $file (@files)
    {

        # check if there are files due to Fortran IO mistakes
        my @forbidden_patterns = map { qr/$_/i } qw( fort\.[0-9]* );
        
        if ("$file" ~~ @forbidden_patterns)
        {
            print "The file $nprocs/$file exists and matches a forbidden pattern.\n";
            $failed = 1;
        }

        # furthermore, check if there is any file with size 0, i.e. a file
        # that has been opened but never been written to.
        if ( "$file" ne ".keepfolder" && ! -s "$nprocs/$file")
        {
            print "The file $nprocs/$file exists and is empty.\n";
            $failed = 1;
        }
    }
    return $failed;

}

#------------------------------------------------------------------------------
# perform checks on the various reference files
#
sub perform_checks
{
    my $nprocs           = shift;
    my $refdir           = shift;
    my @reference_files  = @_;
    my $nreference_files = @reference_files;

    my $failed           = 0;

    print "performing checks...\n";
    if ($nreference_files eq 0)
    {
        print "(no files to check)\n";
    }
    foreach my $file (@reference_files)
    {
        if($file =~ /\.h5/ && $h5diff_avail)
        {
            print "comparison $refdir/$file <-> $nprocs/$file : ";

            # --relative comparison will complain about comparison of zero if very small nonzero values
            #my $h5difference = `h5diff --relative=1.0e-8 $nprocs/$file $refdir/$file`;
            # --delta comparison will complain if very small values are significant
            my $h5difference = `h5diff --delta=1.0e-12 $refdir/$file $nprocs/$file 2>&1`;
            if($? == 0)
            {
                print "OK - no large differences in the HDF5 datasets\n";
            }
            else
            {
                print "BAD!\n";
                print "$h5difference";
                # h5diff reports the actual differences with the
                # option --report .
                $failed = 1;
            }
        }
        elsif($file !~ qr/(\.md5|$input_filename|$filename_regexp_not_to_diff)/)
        {
            print "comparison $refdir/$file <-> $nprocs/$file : ";
            my $difference = `$diff_command $refdir/$file $nprocs/$file 2>&1`;
            #if(compare("$nprocs/$file", "$refdir/$file") == 0)
            if($? == 0)
            {
                print "OK - no difference\n";
            }
            else
            {
                print "BAD!\n";
                print "$difference";
                $failed = 1;
            }
        }
    }
    # check md5 checksums
    my $checksumfile = "$refdir/checksums.md5";
    if( -e $checksumfile)
    {
        chdir($nprocs);
        my $checkmsg = `md5sum --check ../$checksumfile 2>&1`;
        print "MD5 checksums:\n";
        print $checkmsg;
        if ( $? != 0 )
        {
            $failed = 1;
        }
    }
    
    return $failed;
}

#------------------------------------------------------------------------------
# obtain the reference files for this test
#
sub get_reference_folder
{
    my $nprocs = shift;
    my $refdir;

    if ($compare_to_serial == 1 && $nprocs ne "1")
    {
        $refdir = "1";
    }
    else
    {
        $refdir = "reference";
    }

    return $refdir
}



#------------------------------------------------------------------------------
# obtain the reference files for this test
#
sub get_reference_files
{
    my $refdir = shift;
    my @files;
    my @to_test;

    if (!opendir(DIR, "$refdir"))
     {
        print "No reference files for this run\n";
        return 1;
     }
     else
     {
        @to_test = readdir(DIR);
        closedir(DIR);
     }

    foreach my $file (@to_test)
    {
        #GKW specific
        my @gkw_files_not_to_compare = map { qr/$_/i } qw( out input\.out FDS.* FD[0-9].* DM[0-9].* nfs.* gkwdata\.meta \.keepfolder mrad_[Gl].* kx_connect.*  );
        #my @files_not_to_compare = ( "out", "input.out", /FDS.*/, /nfs.*/, "gkwdata.meta", ".keepfolder");
        
        # if it is a regular file and not called 'out'
        if (-f "$refdir/$file" && ! ("$file" ~~ @gkw_files_not_to_compare))
        {
            push(@files, "$file");
        }
    }

    return @files;
}


#------------------------------------------------------------------------------
# obtain the number of processors for this test
#
sub get_procs
{
    # HOW IT IS DONE AT THE MOMENT:
    # return a list of subfolders whose names are made up of numbers only.

    my $dir        = shift;
    my @nprocs;

    opendir(DIR, $dir)
      || die
      "get_procs: can't open $dir\n";
      #
      while (my $file = readdir(DIR))
      {
        if ($file =~ m/^[0-9]+$/ )
        {
            #print " found $file\n";
            push(@nprocs, "$file");
        }
    }
    closedir DIR;

    # # ALTERNATIVE:
    # # parse the line of the input file which contains
    # # the keyword 'nmp'
    # my $nprocs_label        = "nmp";           # no. of procs string in input file


    # my $procs_line = `grep nmp $dir/$input_filename || echo "1"`;

    # $procs_line =~ s/.*\!{2,}.*$//g;    # remove lines with 2 or more "!"
    # $procs_line =~ s/.*nmp\=//g;        # remove everything before "="
    # chomp($procs_line);

    # @nprocs = split(/,/, $procs_line);

    if ($force_serialparallel == $FORCE_SERIAL)
    {
        # remove all elements from the array @nprocs which are different than "1"
        @nprocs = grep { $_ == "1" } @nprocs;
    }

    if ($force_serialparallel == $FORCE_PARALLEL)
    {
        # remove all elements from the array @nprocs which are equal to "1"
        @nprocs = grep { $_ != "1" } @nprocs;
    }

    # make sure the serial run is always executed before the parallel runs
    @nprocs = sort @nprocs;
    
    # return an ARRAY of strings, each element being a number of processors
    return @nprocs;
}


#---------------------------------------------------------------------------
# print a summary of the tests
#
sub print_summary
{

    print "\n------------------ SUMMARY -----------------------------\n";
    print "\n";
    print "# failures: $ntests_failed\n";
    print "   (running): $nfailed_running\n";
    print "  (checking): $nfailed_output\n";
    print "\n";
    if ($ntests_failed > 0)
    {

        if ($nfailed_running > 0)
        {
            print "failed to run:\n";
            print "\n";
            foreach my $testfile (@failed_tests_running)
            {
                print "* $testfile\n";
            }
            print "\n";
        }
        if ($nfailed_output > 0)
        {
            print "failed tests:\n";
            print "\n";
            foreach my $testfile (@failed_tests_output)
            {
                print "- $testfile\n";
            }
        }
    }
    else
    {
        print "All tests successful\n";
    }

}

#---------------------------------------------------------------------------------
# find the GKW executable to test; $GKW_HOME/run
#
sub find_gkw_testing_executables
{

    my $dir = shift;
    my @testing_executables;

    #print "looking for testable executables...\n\n";

    opendir(DIR, $dir)
      || die
      "find_testing_executables: can't open $dir to find executables\n";

    # Look for gkw*.x in the directory and see if it is executable;
    # if so, add it to the list of executables.
    while (my $file = readdir(DIR))
    {
        if ($file =~ /gkw.*\.x/ && -x "$dir/$file")
        {

            #            print " found $file\n";
            push(@testing_executables, "$dir/$file");
        }
    }
    closedir DIR;

    # sort the executables from new to old in modification time
    @testing_executables = sort { -M $a <=> -M $b } @testing_executables;

    return @testing_executables;
}

#---------------------------------------------------------------------------------
# get the right mpirun command
#
sub get_mpiruncmd
{

    # set mpiruncmd via environment variable if present
    if (exists $ENV{GKW_MPIRUNCMD})
    {
        print "Setting mpiruncmd via \$GKW_MPIRUNCMD:\n";
        $mpiruncmd = "$ENV{GKW_MPIRUNCMD}";
        print "mpiruncmd = $mpiruncmd\n";
        return 0;
    } else {
        print "** you may need to set the environment variable\n";
        print "** GKW_MPIRUNCMD to the correct run command\n";
        print "** using default: $mpiruncmd\n";
        #die "done";
    }
}

#-----------------------------------------------------------------------------
#
# This is a 'special testcase' : The external command checks if the GKW source
# code contains only ASCII letters without tabulators, and complains otherwise.
#
sub check_code_is_ascii_notabs
{
    # print `find $code_home -regextype posix-basic -name '[^\#]*[Ff]90'`;
    # print '-----------------';
    # check everything below the code_home folder:
     
    my $extcmd_output = `LC_ALL=C find $code_home/ -regex '.*\/?[a-zA-Z0-9_]*\\.[Ff]90' -exec grep --with-filename --line-number '[^[:print:][:space:]]' '{}' \\;`;
    
    print $extcmd_output;
    return ($extcmd_output eq "");
}

#-----------------------------------------------------------------------------
sub get_executable_capabilities
{

    my ($clean_directory, $executable) = @_;
    chdir($clean_directory) || die "cannot cd to $clean_directory\n";
    if (-d "$input_filename") { unlink "$input_filename"; }
    my $basic_output = `$executable`;
    #print "$basic_output\n";

    #GKW_SPECIFIC:
    if ("$basic_output" =~ m/GKW/)
    {
        print "The executable\n" . " $executable\n" . "appears to be fine.\n";

        #        return 1;
    }
    else
    {
        print "ERROR: The executable "
          . "$executable"
          . " is either broken or is\n";
        print "not GKW. Please check and try again!\n";
        #print "output:\n";
        #print "$basic_output\n";
        die;
    }

    # rather than print the basic output, the information should be put into
    # variables.
    my $capvar;
    my $capval;
    ##while my $line ($basic_output) {
    ##  chomp($line);
    ##    if ( "$line" =~ m/\:/ ) {
    ##        ($capvar, $capval) = split(/:/, $line);
    ##        print "capvar: $capvar - capval: $capval\n";
    ##    }
    ##}
    open(EXOUT, "$executable|");
    my @lines = <EXOUT>;
    close(EXOUT);
    foreach my $line (@lines)
    {
        $line =~ s/\s//g;
        if ("$line" =~ m/\:/)
        {
            ($capvar, $capval) = split(/:/, $line);
            #print "##$line\n";
            chomp($capvar);
            chomp($capval);
            print "$capvar\=$capval\n";
            $cap{$capvar} = $capval;
        }
    }
    #GKW_SPECIFIC:
    my @neededparams = ("MPI", "OPENMP", "MPIRUNCMD", "FFTLIB", "REVISION",
        "REAL_PRECISION", "Eigenvaluesolver");
    foreach my $param (@neededparams)
    {
        die "Missing $param\n" unless defined $cap{$param};
        $cap{$param} = simplify_val($param,$cap{$param});
          print "$param \=$cap{$param}\n";
    }

    #%cap = simplify_hash(%cap);


}

#-----------------------------------------------------------------------------

sub simplify_val
{
    my $param = shift;
    my $val   = shift;

    $val =~ tr/a-z/A-Z/;

    if ( $param eq "MPI" || $param eq "OPENMP" )
    {
        my @acceptable_values = ( "MPI2", "TRUE", "FALSE", "0", "OFF", "ON",
                                   "YES",   "NO", "1", "T", "F", "Y", "N" );
#        die "Value for $param ($val) is bad, acceptable values are"
#            . "@acceptable_values" unless ( grep ($val,@acceptable_values) eq "" );
        $val =~ s/MPI2/1/;
        $val =~ s/TRUE/1/;
        $val =~ s/FALSE/0/;
        $val =~ s/OFF/0/;
        $val =~ s/YES/1/;
        $val =~ s/ON/1/;
        $val =~ s/NO/0/;
        $val =~ s/T/1/;
        $val =~ s/Y/1/;
        $val =~ s/F/0/;
        $val =~ s/N/0/;
    }

    if ( $param eq "FFTLIB" )
    {
        $val =~ s/OFF/0/;
        $val =~ s/NONE/0/;
    }

    if ( $param eq "REAL_PRECISION" )
    {
        $val =~ s/\,.*$//;
        $val =~ s/8/DOUBLE/;
        $val =~ s/4/SINGLE/;
    }

    if ( $param eq "Eigenvaluesolver" )
    {
        $val =~ s/T/1/;
        $val =~ s/N/0/;
    }

    return $val;
}

