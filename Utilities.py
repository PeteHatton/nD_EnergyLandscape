import datetime

def log(caller,message,indent=0):

    now = datetime.datetime.now().strftime('%d/%m/%y, %H:%M:%S')
    ind = ''
    for _ in range(indent):
        ind += "  "
    print('[%s]: %s%s%s >> %s' % (now, "",ind,caller,message))
    

