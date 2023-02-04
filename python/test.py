def get_value_of_variable_from_file(file, file_index, relative_index, string):
    content = [i.strip().split() for i in open(file).readlines()]
    index = [idx for idx, s in enumerate(content) if string in s][file_index]
    
    value = content[index][relative_index]
    
    return value

def find_string_in_file(file, string):
    
    with open(file) as f:
        if string in f.read():
            return True
        else:
            return False

outputCriteria = ["0", "Run successfully completed"]

outputContent = open("slurm_success.out").readlines()[-5].replace("\n","")

test = get_value_of_variable_from_file("output.dat", 0, 0, "successfully")

        
runSuccess = find_string_in_file("output.dat", outputCriteria[1])
        
if outputCriteria[0] in outputContent and runSuccess:
    print("Success")