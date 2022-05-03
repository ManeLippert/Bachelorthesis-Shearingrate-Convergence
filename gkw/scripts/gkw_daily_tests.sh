#!/bin/bash

##############################################################################
##############################################################################
#
##############################################################################
# Summary:
# This is a script for enabling periodic tests (e.g. daily) of gkw, when run
# in a cronjob. See below for an example how to create those.
#
##############################################################################
# Description
#
# This script is intended for automatic testing of gkw.
# For this purpose it is capable of automaticly cloning the git repository
# (read only)/pulling it. The develop branch is checked out. Then either the
# default config file is used or looped over all of them for the machine.
# The standard tests are run. At the end it is tried to make the
# documentation.
# The results of all the makes and tests are logged. The results of the tests
# are also mailed.
# Command line arguments can alter the behaviour in some respects. See the
# 'Options' section below.
#
##############################################################################
# Basic usage:
# Create a cronjob at the machine where the tests should run, for an example
# how to do this see below.
# A .bashrc in the home directory is expected to set the needed environment
# variables. Besides of default ones like PATH,... also the gkw specific
# GKW_HOME, has to be set.
# The latter is used to determine, where to put the folder for the tests
# ($GKW_HOME/../gkw-daily-tests
# Additionaly also EMAIL_ADDRESS should be set, to an adress, where to send
# the mails with the result of this script.
#
# Comands/programs that must be available, besides of those required by gkw:
# grep
# mail
# rm
# sed
#
##############################################################################
# Options
#
# \attention Long option names not supported yet.
#
# -a --all-branches   Check all the branches.
# -b --branch-blacklist file  Check every branch, except those listed in the
#                     file given as argument.
# -c --all-config     Check with every config file in the used config-path.
#                     Suppose the makefile will use as default
#                     config/name/default.mk, then every config file in
#                     config/name/ will be checked for compiling and testing.
# -h --help           Print help for usage and quit.
# -q --quiet          Quiet mode, no general logfile.
# -u --just-update    Just update the branches instead of forcing a clone.
#    --no-clone
# -w --branch-whitelist file  Check every branch listed in the file given as
#                     argument. The filename
#
# Note: only one of a,b,w should be given. If multiple are given the last one
#       will be in effect.
#
##############################################################################
# Example cronjob:
# To add this script as a cronjob use the command 'crontab -e' (e=editing, to
# not delete existing cronjobs).
# Then add a line for calling this script, e.g.
#
# '40   11     *   *    * source ~/.bashrc && nice -n 1 gkw_daily_tests.sh'
#
# This would run the script every day at 11:40 am.
#
# Remainder for the meaning of the colums
#  *     *     *   *    * command to be executed
#  -     -     -   -    -
#  |     |     |   |    |
#  |     |     |   |    +----- day of week (0 - 6) (Sunday=0)
#  |     |     |   +------- month (1 - 12)
#  |     |     +--------- day of        month (1 - 31)
#  |     +----------- hour (0 - 23)
#  +------------- min (0 - 59)
# a '*' means every legal value, list with entries separated by comma are
# possible.
# The 'nice -n 1' sets the priority of the job down, to be still able to use
# the machine during the tests.
# (Note that gkw_run_tests also increases the priority, thus we use only 1
# here.)
#
##############################################################################
# Todos
#
# \todo   9 Enable configuration of BPATH. At least make it user independent,
#           ${GKW_HOME}/../gkw-daily-tests might be appropriate.
#           Or make a check if it is defined, if not, then use
#           ~/gkw-daily-tests.
#           +++ Changed to the second, seems to work.
# \todo   9 get_config should be made more robust.
#
# \todo  10 If compiling fails, check previous revisions (up to last one
#           checked) until a working one is found.
#           Might also be usefull for testing, but as there are usally errors,
#           this won't work yet.
# \todo  10 Error handling in clone.
# \todo  10 Working iteration over config files, even when there are spaces in
#           the pathname.
#           Setting IFS to just linebreaks and sort them again into lines
#           would work. If IFS could be reset just for setting up the loop,
#           then this would be an acceptable solution.
#           Alternativly use an array?
# \todo  10 Enable testing of a list of branches (all branches).
#           It should also be possible to give a list of branches to check
#           (via file), a list of branches that should not be checked (via
#           file).
#           Must be ensured that branches will not appear twice, because of
#           already existing 'links'.
#           +++ Add function for getting list of branches.
#           +++ Add commandline arguments for selecting which branches to
#               check.
#           +++ Add functionality for black/white lists. At the moment, the
#               path has to be relative to $GKW_HOME.
#           +++ Not using a branch twice, should be avoided by skipping HEAD
#               and using only the names of remote branches.
#           But looping over the branches is so far not done. Especially the
#           messages have to be adapted.
#           --- Inserted corresponding code, but not tested yet.
#
# \todo  11 Implement option for verbosity.
#
# \todo  12 Implement option -q.
#               Need to be tested.
#               Current implementation broken, as always $BPATH/$LOGFILE is
#               used and the latter is set to /dev/null.
# \todo  12 Implement version of getting config files, that does not depend on
#           the name of the default config file.
#
##############################################################################
# Assumed as done:
# \todo -10 Check that the setting of the nice level is working.
#           --- Not as expected, nice level is 19, not 12.
#           --- A value of 2, results in a nice leve of 17.
#           +++ Setting of the nice level is working, reason for the higher
#               than expected results, is an increase of the nice level in
#               gkw_run_tests.py for running the tests.
# \todo -10 Finish the switch to git.
#           - Check if removing is complete
#           +++ Should be the case.
#           - Getting revision number.
#           +++ Changed from revision number to commit hash.
#           - Check if subdirectory is still needed (probably, so:better name).
#           +++ Name changed.
#           - Adapting messages (mail).
#           +++ Done.
# \todo   5 Test not only if the folder is there, but also if there is already
#           a checkout (and not only an empty folder).
#           --- Code included, but it is not tested yet.
#           +++ Seems to work, as the current version has used this successful.
# \todo   5 Try to get the mail-adress from somewhere else. (Otherwise this
#           will lead to trouble every time someone commits his changes.)
#           --- Implemented a check for a global variable. If these is not set,
#               the mail will be send to the mail adress 'mail' is using.
#               Name generation seems to work, the rest is not tested yet.
#           +++ Default case with name generation tested successfuly.
#           +++ Using of the environment variable successful.
# \todo   0 Config dependent flags for compiling and testing needed. So far
#           the global ones would be reset (if first fails and last one is
#           successfull, this would be reported as success).
#           This includes setting of the mail-text.
#           +++ Actually mail-text seems to be ok.
#           +++ The compiling text for successful compiling still overwrote
#               the text, has been fixed.
#           --- One failed compile still prevents other tests.
#           +++ Seems to be fixed.
# \todo   9 Send mail also if all is okay (containing at least the checked
#           revision number).
#           --- Code included, but it is not tested yet. It is expected that
#               the mail will be send, as this should not differ from the
#               case of failure, but if the text is correct has to be checked.
# \todo   9 Consider also the config file for the output/messages of the
#           functions.
#           +++ Done, seems to work.
# \todo  10 Do not send more than one mail?
#               Implemented, seems to work, but more testing needed.
#               (Espacially the cases that testing is not broken, and that
#               making doc is broken.)
# \todo  10 Enable testing using different configuration files.
#           Think of how the get the list of configuration files. Everything
#           that is in the corresponding config folder?
#           FILE=`make info | grep "CONFIGFILE =" | sed -e 's/CONFIGFILE = //'`
#           ${FILE%/*}
#           Do one run, the using `grep -A 1 "Attempting to build with" trunk_make.log | grep "/"`
#           to get the config file used, including location. List the files in
#           the base directory, ignore default.mk.
#           Alternatively use
#           make --dry-run 2> /dev/null | grep -A 1 "Attempting" | grep 'echo "   ' | sed 's/echo \"   //'
#           to get the used default config file. Get path from these and use
#           all files in the path.
#           +++ Implemented. Using dry run to get the path and ls to get the files.
# \todo  10 Switch to enforce checkout (and deleting of the old source before).
#           +++ At the moment checkout is done every time, as the old code is
#               deleted. Might even be preferable.

##############################################################################
# To have environment variables defined (e.g. PATH).
source ~/.bashrc

##############################################################################
# Variables

# Programs and theire options.

RM=rm
RM_RECURSIV="-r"

# If the various tests have been successful (=0). The 'TEMP' variants are for
# the test with different config files.
# Starting with false value (=1) would be safer, but more complicated for
# compile and tests, because of the possibility of multiple configuration
# files (failure could be overwritten).
FLAG_SCOMPILE=0
FLAG_SCOMPILE_SINGLE=0
FLAG_STESTS=0
FLAG_SDOC=1

# Logfile name, might be changed due to comandline arguments.
LOGFILE="out.log"

# Revisionfile is a parameter, the revisionnumber determined later.
REVISIONNUMBER=""
REVISIONFILE=".gkw-daily-tests-rev"

# Base path to use
#BPATH="${HOME}/Programme/gkw-daily-tests"
BPATH="${GKW_HOME}/../gkw-daily-tests"

# To whom should the mails go.
if [[ "x$EMAIL_ADDRESS" == "x" ]] ; then
  # Environment variable not set, use a default instead.
  TOMAILADDRESS="`whoami`@`uname -n`"
else
  # Environment variable is set, use it.
  TOMAILADDRESS="${EMAIL_ADDRESS}"
fi

# For creating the subject of the mail.
SUBJECTBASE="[gkw-daily-tests]"
SUBJECT=""

# For storing the text for the different parts.
MAILTEXT_COMPILE=""
MAILTEXT_TESTS=""
MAILTEXT_DOC=""
# These seperators should enhance readibility, and are also used for the
# output of the script itself (which is usally written to a file).
MAILTEXT_SEP_PART="============================================="
MAILTEXT_SEP_BRANCH="#############################################"
MAILTEXT_SEP_CONF="------------------------------"

CHECKTRUNK="true"
CHECK_ALL_CONFIG="false"
DEFAULT_CONFIG="default.mk"
CONFIG_FILES=""
TEST_BRANCH="develop"
# Default branch for checkout and testing the code.
DEFAULT_BRANCH="develop"
BRANCH_LIST=""
BRANCH_FILE=""
# Name of the subdirectory to create. Using this makes removing easier, also
# so there are only the main-log files in the main directory.
TEST_DIRECTORY="work_dir"

# Variables for commandline options.
JUST_UPDATE="false"
CHECK_PREVIOUS="false"
TEST_ON_FAILED_COMPILE="false"
CHECK_ALL_BRANCHES="false"

# Error return values.
ERROR_BRANCH_DOES_NOT_EXIST=5

##############################################################################
# functions
#
# Functions are used as they structure the code and thus provide and advantage
# even if they are used only once.

# Set options.
function set_options {

  # Process the options. Set variables.
  while getopts ":ab:chquw:" opt; do
    case $opt in
      # Check compiling and testing for all branches.
      a)
        CHECK_ALL_BRANCHES="true"
        ;;
      # Check compiling and testing for all branches, with exceptions listed
      # in a file.
      b)
        CHECK_ALL_BRANCHES="black"
        BRANCH_FILE="$OPTARG"
        ;;
      # Check compiling and testing for all the configuration files in the
      # default folder.
      c)
        CHECK_ALL_CONFIG="true"
        ;;
      # Print help.
      h)
        echo "-h was triggered!"
        help # function call.
        exit 0
        ;;
      # Activate quiet mode.
      q)
        # quiet mode. Instead of logging to a file, discard the output.
        #echo "-q was triggered!"
        LOGFILE="/dev/null"
        ;;
      # Just update the repository, do not remove it, to force a new 'clone'.
      u)
        #echo "-u was triggered!"
        JUST_UPDATE="true"
        ;;
      # Check compiling and testing for all branches listed in a file.
      w)
        CHECK_ALL_BRANCHES="white"
        BRANCH_FILE="$OPTARG"
        ;;
      \?)
        #echo "Invalid option: -$OPTARG"
        ;;
      :)
        #echo "Option -$OPTARG requires an argument."
        exit 1
        ;;
    esac

  done

}

# Print help: Simply print the head of this file.
function help {
  head -n 91 "$0"
}

# Check if a certain path does exist and change to it.
# If the path does not exist, create it, before changing to it.
function check_for_path {
  CHECK_PATH="$1"
  if [ ! -d "$CHECK_PATH" ] ; then
    mkdir -p "${CHECK_PATH}"
    if [ ! -d "$CHECK_PATH" ] ; then
      echo "Error: Could not create Path ${CHECK_PATH}"
      exit -1
    fi
    echo "$CHECK_PATH had to be created."
  fi
  cd ${CHECK_PATH}
  pwd
}

# Cloning repository or updating.
function clone_rep {
  BRANCH="$1"

  # `git log ...` is called, (1) to have the hash in the logfile, as it
  # might be helpfull for debugging, (2) to check if there is already a
  # checked out version of the code.
  git log --pretty=format:'HASH %H' -n 1 2> /dev/null
  SUCESSFULLINFO="$?"
  echo "" # No automatic endline is printed after the hash.
  if [[ "x$SUCESSFULLINFO" == "x0" ]] ; then
    echo "Existing repository found, updating it."
    git pull 2>&1
  else
    echo "No existing repository found, checking out new one."

    # Check out the trunk if it is explicitly requested, or if BRANCH is empty.
    # Better use the https version, the ssh version requires a key.
    # Passing -q to surpress progress report.
    git clone -q https://bitbucket.org/gkw/gkw.git ./
    # Make sure that the hash ('revisionnumber') is in the log.
    git log --pretty=format:'HASH %H' -n 1 2> /dev/null
    echo "" # No automatic endline is printed after the hash.
  fi

  checkout_branch "$BRANCH" # function call
  #if [[ "x$BRANCH" == "x" ]] ; then
    ## Per default develop is checked, not master.
    ##git checkout -b "$DEFAULT_BRANCH" "origin/$DEFAULT_BRANCH" 2> /dev/null
    #git checkout "$DEFAULT_BRANCH" 2> /dev/null
  #else
    #git branch -a | grep "$BRANCH" 1> /dev/null
    #BRANCH_EXISTS="$?"
    #if [[ "x$BRANCH_EXISTS" == "x0" ]] ; then
      ##git checkout -b $BRANCH origin/$BRANCH
      #git checkout $BRANCH
    #else
      #echo "ERROR checkout: Branch '$BRANCH' not found."
      #exit $ERROR_BRANCH_DOES_NOT_EXIST
    #fi
  #fi

  # Actually used by gkw to define the 'revision', but more difficulty to
  # extract. Thus it is set here directly.
  git describe --always --tags
  REVISIONNUMBER=`git describe --always --tags 2> /dev/null`
}

# Checkout a branch. If none is given a default is used.
function checkout_branch {
  BRANCH="$1"

  if [[ "x$BRANCH" == "x" ]] ; then
    # Per default develop is checked, not master.
    #git checkout -b "$DEFAULT_BRANCH" "origin/$DEFAULT_BRANCH" 2> /dev/null
    git checkout "$DEFAULT_BRANCH" 2> /dev/null
  else
    git branch -a | grep "$BRANCH" 1> /dev/null
    BRANCH_EXISTS="$?"
    if [[ "x$BRANCH_EXISTS" == "x0" ]] ; then
      #git checkout -b $BRANCH origin/$BRANCH
      git checkout $BRANCH 2> /dev/null
    else
      echo "ERROR checkout: Branch '$BRANCH' not found."
      exit $ERROR_BRANCH_DOES_NOT_EXIST
    fi
  fi
}

# This function will get the available configuration files and store the
# result in the global variable CONFIG_FILES, if CHECK_ALL_CONFIG is set (-c).
# Otherwise it will be set only to the default.
# \attention At the moment this relies heavyly on specifics of the make
#   process and the output.
function get_configs {
  # Assumes that default.mk is actually the default for the config file name,
  # and that a config file of this name exists.
  if [[ "x$CHECK_ALL_CONFIG" == "xtrue" ]] ; then
    CONFIG_PATH_TEMP=`make --dry-run 2> /dev/null | grep -A 1 "Attempting" | grep 'echo "   ' | sed 's/echo \"   //' | sed 's/default.mk\"//'`
    CONFIG_FILES_TEMP=`ls --ignore=compiler "$CONFIG_PATH_TEMP/"`
  else
    CONFIG_FILES_TEMP=""
  fi

  #CONFIG_PATH_TEMP="/home/heinz/programm gkw/"
  LPATH=`pwd`
  CONFIG_PATH_TEMP=`echo "$CONFIG_PATH_TEMP" | sed "s%$LPATH/%%"`

  IFS_OLD="$IFS"
  IFS="
"
  for CONFIG_TEMP in $CONFIG_FILES_TEMP ; do
    CONFIG_FILES="$CONFIG_FILES$CONFIG_PATH_TEMP/$CONFIG_TEMP
"
  done
  IFS="$IFS_OLD"

  echo "List of config files to use:"
  echo "$MAILTEXT_SEP_CONF"
  echo "$CONFIG_FILES"
  echo "$MAILTEXT_SEP_CONF"
}

# Compile the code. If this fails and it is requested, go to previous versions
# until compiling is succesfull.
function make_compile {
  GKW_HOME_TEMP="$1"
  BRANCH="$2"
  CONFIGFILE="$3"

  TEMP_REV="$REVISIONNUMBER"

  TEMP_BASE=`basename "$CONFIGFILE"`

  make clean
  make GKW_HOME="$GKW_HOME_TEMP" \
       CONFIG="$CONFIGFILE" 1> "${GKW_HOME_TEMP}/make_$TEMP_BASE.log" 2>&1
  CSUCESS="$?"

  if [[ "x$CSUCESS" == "x0" ]]; then

    FLAG_SCOMPILE_SINGLE=0
    # Set the variables for the text of the mail.
    MAILTEXT_COMPILE="$MAILTEXT_COMPILE
Compiling with CONFIG=$CONFIGFILE successful."

  else
    FLAG_SCOMPILE=1
    FLAG_SCOMPILE_SINGLE=1
    MAILTEXT_COMPILE_TEMP=`tail "${GKW_HOME_TEMP}/make_${TEMP_BASE}.log"`
    MAILTEXT_COMPILE="$MAILTEXT_COMPILE
Compiling with CONFIG=$CONFIGFILE broken and no compiling revision found.
CHECK_PREVIOUS=$CHECK_PREVIOUS
${MAILTEXT_COMPILE_TEMP}"

    MAILTEXT_TESTS="
Tests can only be made if compiling was successful.
"

    echo "" >> "${GKW_HOME_TEMP}/make_$TEMP_BASE.log"
    echo "---- Compiling broken. ----" >> "${GKW_HOME_TEMP}/make_$TEMP_BASE.log"
  fi

  cat "make_$TEMP_BASE.log" # Output of this function.
}

# Checking a single revision of the code.
# Assuming path and revision (of the code) have been set.
function make_tests {
  GKW_HOME_TEMP="$1"
  BRANCH="$2"
  CONFIGFILE="$3"

  TEMP_BASE=`basename "$CONFIGFILE"`

  # The make testcommand must set GKW_HOME properly, otherwise the lwrong
  # executional is taken.
  # Also we write the output to a special file and to the generic logfile of
  # this script.
  make tests GKW_HOME=${GKW_HOME_TEMP} 1> "${GKW_HOME_TEMP}/tests_$TEMP_BASE.log" 2>&1
  TSUCESS=`grep "# failures: [1-9]" "${GKW_HOME_TEMP}/tests_$TEMP_BASE.log"`
  cat "${GKW_HOME_TEMP}/tests_$TEMP_BASE.log" # Output of this function.

  # Print a message about the result of the tests and set the variables for
  # the mail.
  if [[ "x$TSUCESS" == "x" ]]; then
    echo "++++ Tests with CONFIG=$CONFIGFILE successfull. ++++"
    MAILTEXT_TESTS="$MAILTEXT_TESTS
Tests with CONFIG=$CONFIGFILE successful."
  else
    FLAG_STESTS=1
    echo ""
    echo "---- Testing with CONFIG=$CONFIGFILE: ----"
    echo "---- One/Some test(s) is/are broken. ----"

    # Append to the mailtext the summary part of the tests output. The sed
    # command will extract this.
    MAILTEXT_TESTS_TEMP=`cat ${GKW_HOME_TEMP}/tests_$TEMP_BASE.log | sed -n -e '/SUMMARY/,$p'`
    MAILTEXT_TESTS="$MAILTEXT_TESTS

   Tests with CONFIG=$CONFIGFILE broken:
${MAILTEXT_TESTS_TEMP}"
    #alternative: awk 'BEGIN{ found=0} /SUMMARY/{found=1}  {if (found) print }'

  fi
}

# Check if making of the documentation works.
function make_docs {
  GKW_HOME_TEMP="$1"
  BRANCH="$2"

  # To be sure, the repository is on the latest version.
  echo "=== Updating for creating docs ===" "" >> ${BPATH}/${LOGFILE}
  clone_rep "$BRANCH" >> ${BPATH}/${LOGFILE} 2>&1 # function call

  make docs 1> "${GKW_HOME_TEMP}/${BRANCH}_docs.log" 2>&1
  DSUCESS="$?"
  #cat ${GKW_HOME_TEMP}/docs.log # Output of this function.

  # Set the variables for the mail.
  if [[ "x$DSUCESS" == "x0" ]]; then
    echo "docs successfully made"
    FLAG_SDOC=0
    MAILTEXT_DOC="make doc successful."
  else

    MAILTEXT_DOC=`tail ${GKW_HOME_TEMP}/${BRANCH}_docs.log`
    MAILTEXT_DOC="
make doc broken.
${MAILTEXT_DOC}"
  fi
}

# Get a list with the available branches.
# Checks the settings that have been made and sets the variable BRANCH_LIST
# accordingly.
# This function will print a message which option is used.
function get_branches {
  ## This version would not see any branches, if there is not already a local
  ## version. Using -a would lead to having branches twice in the list, if
  ## there is a local copy.
  #BRANCH_LIST=`git branch | sed 's/*//' | sed 's/ *//'`
  if [[ "x$CHECK_ALL_BRANCHES" == "xtrue" ]] ; then
    echo "Branches: Using all."
    # As all branches, is defined here, all the items of branch -a, that are
    # at the remote side, with the special entry HEAD removed (as this refers
    # to master and thus doubles another entry).
    BRANCH_LIST=`git branch -a | grep "remotes" | sed 's/remotes\/origin\///' | grep -v "HEAD"`
  elif [[ "x$CHECK_ALL_BRANCHES" == "xblack" ]] ; then
    echo "Branches: Using blacklist from '${GKW_HOME}/$BRANCH_FILE'."
    # First get all the branches, then for each element in the file, remove
    # it, if it is in the list. This way nonexistent or wrong spelled
    # branches are ignored.
    BRANCH_LIST=`git branch -a | grep "remotes" | sed 's/remotes\/origin\///' | grep -v "HEAD"`
    for REM_BRANCH in `cat ${GKW_HOME}/$BRANCH_FILE` ; do
      BRANCH_LIST=`echo "$BRANCH_LIST" | grep -v "$REM_BRANCH"`
    done
  elif [[ "x$CHECK_ALL_BRANCHES" == "xwhite" ]] ; then
    echo "Branches: Using whitelist from '${GKW_HOME}/$BRANCH_FILE'."
    # The branch list is simply the content of the file. This means, that
    # one should take care, that the list contains, no non-existent entries,
    # i.e. because of a spelling mistake.
    BRANCH_LIST=`cat "${GKW_HOME}/$BRANCH_FILE"`
  else
    echo "Branches: Using default."
    # No option regarding the branch to check has been given, so use develop
    # as default.
    BRANCH_LIST="$DEFAULT_BRANCH"
  fi
  ## Testoutput
  #for LOC_BRANCH in $BRANCH_LIST ; do
    #echo "branch: $LOC_BRANCH"
  #done
}

##############################################################################
# main part

set_options "$@" > ${BPATH}/${LOGFILE} # function call, gets commandline options passed to be able to parse them in the function.

STARTDATE=`date`
echo "$STARTDATE" >> ${BPATH}/${LOGFILE}

# Test if BPATH already exists.
# - yes, change to it.
# - no,  create it, then change into.
check_for_path "$BPATH" >> ${BPATH}/${LOGFILE} # function call

# Remove the branch, to force a clean cloning of the repository.
if [[ "x$JUST_UPDATE" == "xfalse" ]] ; then
  # Remove only if it exists.
  if [ -d "$TEST_DIRECTORY" ] ; then
    #echo "removing branch"
    $RM $RM_RECURSIV "$TEST_DIRECTORY/" >> ${BPATH}/${LOGFILE}
  fi
fi

# Check if the folder for the trunk already exists and change to it.
check_for_path "$TEST_DIRECTORY" >> ${BPATH}/${LOGFILE} # function call.

clone_rep >> ${BPATH}/${LOGFILE} # function call

get_branches >> ${BPATH}/${LOGFILE} # function call

# FOR every BRANCH
# 1. Checkout the branch
# 2. Compile the code all new.
# 3. Run tests
for TEST_BRANCH in $BRANCH_LIST ; do
  checkout_branch "$TEST_BRANCH" >> ${BPATH}/${LOGFILE} # function call
  echo "" >> ${BPATH}/${LOGFILE}
  echo "$TEST_BRANCH" >> ${BPATH}/${LOGFILE}
  echo "$MAILTEXT_SEP_BRANCH" >> ${BPATH}/${LOGFILE}
  MAILTEXT_COMPILE="${MAILTEXT_COMPILE}

$TEST_BRANCH
$MAILTEXT_SEP_BRANCH
"
  MAILTEXT_TESTS="${MAILTEXT_TESTS}

$TEST_BRANCH
$MAILTEXT_SEP_BRANCH
"
  MAILTEXT_DOC="${MAILTEXT_DOC}

$TEST_BRANCH
$MAILTEXT_SEP_BRANCH
"

  ## As the revision number may be needed for naming files, etc., it is stored
  ## in a variable.
  #REVISIONNUMBER=`cat ${BPATH}/${LOGFILE} | grep "HASH " | grep -o "[0-9a-f]*"`
  ## At least temporary moved into clone_rep funktion.

  get_configs >> ${BPATH}/${LOGFILE} # function call

  ## Testoutput
  #echo "#####" >> ${BPATH}/${LOGFILE}
  #echo "$CONFIG_FILES" >> ${BPATH}/${LOGFILE}
  #echo "#####" >> ${BPATH}/${LOGFILE}

  if [[ "x$CONFIG_FILES" == "x" ]] ; then
    echo "==================================================" >> ${BPATH}/${LOGFILE}
    echo "Using now '$CONFIG' for compiling and testing."     >> ${BPATH}/${LOGFILE}
    echo "==================================================" >> ${BPATH}/${LOGFILE}

    make_compile "${BPATH}/$TEST_DIRECTORY" "$TEST_BRANCH" "" >> ${BPATH}/${LOGFILE} # function call

    if [ ${FLAG_SCOMPILE_SINGLE} -eq 0 ] ; then
      make_tests "${BPATH}/$TEST_DIRECTORY" "$TEST_BRANCH" "" >> ${BPATH}/${LOGFILE} # function call
    fi
  else
    # Looping over all the configuration files.
    IFS_OLD="$IFS"
    IFS="
"
    for CONFIG in $CONFIG_FILES ; do
      echo "==================================================" >> ${BPATH}/${LOGFILE}
      echo "Using now '$CONFIG' for compiling and testing."     >> ${BPATH}/${LOGFILE}
      echo "==================================================" >> ${BPATH}/${LOGFILE}

      make_compile "${BPATH}/$TEST_DIRECTORY" "$TEST_BRANCH" "$CONFIG" >> ${BPATH}/${LOGFILE} # function call

      if [ ${FLAG_SCOMPILE_SINGLE} -eq 0 ] ; then
        make_tests "${BPATH}/$TEST_DIRECTORY" "$TEST_BRANCH" "$CONFIG" >> ${BPATH}/${LOGFILE} # function call
      fi

    done

    IFS="$IFS_OLD"
  fi

  #make_docs "${BPATH}/$TEST_DIRECTORY" "$TEST_BRANCH" >> ${BPATH}/${LOGFILE} # function call
  make_docs "${BPATH}/$TEST_DIRECTORY" >> ${BPATH}/${LOGFILE} # function call

done # for every branch

# Tests are done, now the mail is build and send.

# Set the subject. Obviously the whole test was only successful only if each
# of the subtests was successful.
if [ ${FLAG_SCOMPILE} -eq 0 ] && [ ${FLAG_STESTS} -eq 0 ] && [ ${FLAG_SDOC} -eq 0 ] ; then
  SUBJECT="${SUBJECTBASE} r${REVISIONNUMBER} Tests successful"
else
  SUBJECT="${SUBJECTBASE} r${REVISIONNUMBER} Tests failed"
fi

mail -s "${SUBJECT}" ${TOMAILADDRESS} << EOF
${MAILTEXT_COMPILE}
${MAILTEXT_SEP_PART}
${MAILTEXT_TESTS}
${MAILTEXT_SEP_PART}
${MAILTEXT_DOC}
EOF

#date >> ${BPATH}/${LOGFILE}
ENDDATE=`date`
echo "$ENDDATE" >> ${BPATH}/${LOGFILE}

cd ${BPATH}
if [[ ! "x${LOGFILE}" == "x/dev/null" ]] ; then
  mv "${LOGFILE}" "r${REVISIONNUMBER}_${LOGFILE}"
fi

##############################################################################
##!/bin/bash
#set -x
#function pass_back_a_string() {
#    eval "$1='foo bar rab oof'"
#}
#
#return_var=''
#pass_back_a_string return_var
#echo $return_var

# ``Flow chart'', items marked with a '+' have a function/code (might still be untested).
# + check the settings
# + create the main directory if neccessary
# + remove the code, if not told otherwise
# + check out the code/update it.
# - go through every branch (loop)
#   - for every config file (for this machine) do
#     + compile the code.
#     + run the tests. (If compiling was not successful the first time: only when asked so.)
#   + make the documentation. (As this does not depend on compiling/config, do it in every case.)
# - send a mail with the summary.
