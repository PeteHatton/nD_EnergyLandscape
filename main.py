#!/opt/homebrew/Caskroom/miniconda/base/bin/python3

import TwoD_EL.Minimize as mn
import TwoD_EL.Dimer as dm
import TwoD_EL.Utilities as ut

import sys
import importlib
import textwrap


availableCommands = {
    "Minimize": "Minimize",
    "Dimer": "Dimer",
}

commandInfo = {
    "Minimize": "Minimize a system",
    "Dimer": "Run the Dimer method",
}

requiredLibs = [
    ("numpy", None, None),
    ("matplotlib", None, None),
]

################################################################################

def usage():
    availCmdString = ""
    for cmd in sorted(availableCommands.keys()):
        if cmd == "Client":
            pass
        else:
            if cmd in commandInfo.keys():
                availCmdString += "%s - %s\n                      " % (cmd, commandInfo[cmd])
            else:
                availCmdString += "%s\n                      " % cmd
    
    usageString = """
                  Usage: main.py NAME_OF_MODULE [OPTIONS] [ARGUMENTS]
                       : main.py NAME_OF_MODULE -h will print help"
                  
                  NAME_OF_MODULE must be one of:
                      
                      %s""" % availCmdString
    
    return textwrap.dedent(usageString)

################################################################################

def checkRequirements():
    """
    Check required libs are installed
    
    """
    errors = []
    for lib, libver, minver in requiredLibs:
        try:
            _module = importlib.import_module(lib)
        
        except ImportError:
            errors.append("Could not find required lib: %s" % lib)
        
        else:
            if libver is not None and minver is not None:
                getver = getattr(_module, libver, None)
                if getver is not None:
                    if callable(getver):
                        installedver = getver()
                    else:
                        installedver = getver
                    
                    ok = versionTuple(installedver) >= versionTuple(minver)
                    
#                     print "LIB %s: installed ver %s; min ver %s; ok %s" % (lib, installedver, minver, str(ok))
                    
                    if not ok:
                        errors.append("Min version not satisfied for lib: %s (%s < %s)" % (lib, installedver, minver))
    
    if len(errors):
        sys.stderr.write("ERROR: requirements not satisfied (details below):\n\n")
        sys.exit("\n".join(errors) + "\n")

################################################################################

def main():
    checkRequirements()
    
    if len(sys.argv) == 2 and (sys.argv[1] == "-v" or sys.argv[1] == "--version"):
        print("LKMC %s" % LKMC.__version__)
        sys.exit(0)
    
    if len(sys.argv) < 2 or sys.argv[1] == "-h" or sys.argv[1] not in availableCommands.keys():
        sys.exit(usage())
    
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
        print(usage())
        sys.exit(0)
    
    if sys.argv[1] != "Client":
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('~~~~~~~~~~~~~~~~~~~~~~ 2D - Energy Landscape ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        
    
    _module = importlib.import_module(".%s" % (availableCommands[sys.argv[1]],), package="TwoD_EL")
    
    modname = sys.argv.pop(1)
    sys.argv[0] += " %s" % modname
    
    if modname != "Client":
        ut.log("2D_EL", "Running: %s" % modname, 0)
    
    if hasattr(_module, "main"):
        status = _module.main()
    
    else:
        print("ERROR: cannot run: %s" % modname)
        status = 254
    
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    sys.exit(status)

################################################################################

if __name__ == '__main__':
    main()
