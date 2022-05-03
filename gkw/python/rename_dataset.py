#!/usr/bin/env python3

import h5py
import sys
import re
import functools

def print_usage():
    print("""  Usage:
        %s [--dry] <h5filename> <orig_dataset_regex> <new_dataset_name>

    This renames all datasets matching the regular expression.
    The new name may contain groups.

    Example:
      # rename all datasets ending on 'weigh', to end on 'weight'
      ~/gkw/python/rename_dataset.py gkwdata.h5 '(.*)weigh$' '\\1weight'
    
"""
          % (sys.argv[0]))


class rename_visitor:

    objects_to_rename = {}
    
    def __call__(self, name, object, orig_dsetname_prog, new_dsetname_pattern, dry):
        if(isinstance(object, h5py.Dataset) and orig_dsetname_prog.match(name)):
            new_dsetname = orig_dsetname_prog.sub(new_dsetname_pattern, name)
            print("%s ... %s" % (name,new_dsetname))
            # it does not seem to be possible to rename the dataset
            # while visiting. So just make a list and rename later
            self.objects_to_rename[name] = new_dsetname

    def execute_renaming(self, h5file, dry):
        for name in self.objects_to_rename:
            print("%s => %s" % (name, self.objects_to_rename[name]))
            if(not dry):
                h5file[self.objects_to_rename[name]] = h5file[name]
                del h5file[name]

if __name__ == "__main__":
    argv = sys.argv[1:]

    if(argv[0] == '--help'):
        print_usage()
        exit(0)

    
    if(argv[0] == '--dry'):
        print("Dry run, no changes are applied.")
        dry = True
        argv = argv[1:]
    else:
        dry = False
    
    if(len(argv) == 3):
        filename = argv[0]
        orig_dsetname = argv[1]
        new_dsetname_pattern = argv[2]
        try:
            h5file = h5py.File(filename, "r+")
            orig_dsetname_prog = re.compile(orig_dsetname)
            visitor = rename_visitor()
            h5file.visititems(functools.partial(visitor,
                                                orig_dsetname_prog=orig_dsetname_prog,
                                                new_dsetname_pattern=new_dsetname_pattern,
                                                dry=dry))
            visitor.execute_renaming(h5file, dry)
            h5file.close()
        except Exception as err:
            print(err)
            exit(1)
    else:
        print_usage()
        exit(1)
