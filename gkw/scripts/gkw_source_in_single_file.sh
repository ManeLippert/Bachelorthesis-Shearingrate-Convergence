#!/bin/sh
one_file='gkw_in_one_file.F90'
rm -f $one_file

# concatenate lines ending with backslash | grep certain line | remove everything up to the equal sign
for objfile in $(perl -p -e 's/\\\n//' $GKW_HOME/src/objfiles.mk  | grep 'F90OBJ =' | sed 's/^.*=//')
do
    cat ${objfile%o}[fF]90 >> $one_file
done

echo "All GKW source files were concatenated to $one_file"
