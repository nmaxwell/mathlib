
#http://github.com/nathanial/pylas/blob/master/las/parser.py

def line_parse( string, delims = [',', '\t' ], sep = [' '] ):
    
    for c in delims:
        string = string.replace( c ,' ' )
    
    try:
        return float(string)
    except:
        R = []
        for x in string.split():
            try:
                R.append(float(x))
            except:
                pass
        
        return R


def read_file(file_string, delims = [',', '\t' ], sep = [' '] ):
    file_ptr = file(file_string)
    lines = file_ptr.readlines()
    
    R = [ line_parse(line, delims, sep ) for line in lines ]
    
    while [] in R:
        R.remove([])
    
    return R

def read_text( text, delims = [',', '\t' ], line_delim = [';' ] ):
    
    if isinstance(text,str):
        for c in line_delim:
            text = text.replace( c ,'\n' )
        R = [ line_parse(line, delims, sep ) for line in text.split('\n') ]
        
        while [] in R:
            R.remove([])
        
        return R




if __name__ == '__main__':
    for i in read_file(open(sys.argv[1])):
        print i

