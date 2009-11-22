




def read_file(file_string):
    file_ptr = file(file_string)
    lines = file_ptr.readlines()
    data = []
    [data.append(float(line.strip())) for line in lines]
    return data
    





if __name__ == '__main__':
    for i in read_file(open(sys.argv[1])):
        print i

