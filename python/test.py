def find_string_in_file(filename, string):
    
    with open(filename) as f:
        if string in f.read():
            return True
        else:
            return False

def get_error_type(filename):
    
    slurm_errors = {"executable":["error on file ./gkw.x (No such file or directory)", "No executable found"],
                    "walltime":["process killed (SIGTERM)", "Exceeded wall time"],
                    "timeout":["DUE TO TIME LIMIT", "Exceeded time limit"],
                    "config":["couldn't open config directory", "Config not loading"],
                    "hdf5":["HDF5-DIAG", "Writing h5 file failed"]}
    
    for key in slurm_errors:
        if find_string_in_file(filename, slurm_errors[key][0]):
            return slurm_errors[key][1]
    
    return "Unknown error occurred"

print(get_error_type("slurm_nofile.out"))