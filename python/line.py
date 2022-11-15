filename = 'status.txt'

def delete_line_in_file(file, start=None, end=None):
    f = open(file)

    lines = f.readlines()[start:end]
    
    try:
        #lines[-1] = lines[-1].split('\n')[0]

        content = ''.join(lines)

        f = open(file, 'w+')
        f.write(content)
        f.close()
    except IndexError:
        pass
    

delete_line_in_file(filename, end=-4)