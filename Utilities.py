import datetime

def log(caller,message,indent=0):

    now = datetime.datetime.now().strftime('%d/%m/%y, %H:%M:%S')
    ind = ''
    for _ in range(indent):
        ind += "  "
    print('[%s]: %s%s%s >> %s' % (now, "",ind,caller,message))
    
def ORlC(filePath,popNo=0):
    with open(filePath,'r') as f:
        lines = f.readlines()
        [ lines.pop(0) for _ in range(popNo) ]
    return lines
    

