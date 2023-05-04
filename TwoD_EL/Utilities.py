import datetime
import sys

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
    

def convertStrToType(strValue, typeName):
    try:
        if typeName == 'str':
            return strValue.strip()
        elif typeName == 'int':
            return int(strValue)
        elif typeName == 'float':
            return float(strValue)
        elif typeName == 'list':
            return [ float(coord) for coord in strValue.split(',') ] 
        else:
            sys.exit(__name__ + ": convertStrToType: ERROR: undefined type: " + typeName)
    except:
        sys.exit(__name__ + ": convertStrToType: ERROR: cannot convert: " + strValue + " to " + typeName)
