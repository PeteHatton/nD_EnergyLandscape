#!/opt/homebrew/Caskroom/miniconda/base/bin/python3

import TwoD_EL.Minimize as mn
import TwoD_EL.Dimer as dm
import sys

"""
To-Do:
    - use argparse or something to allow for more command line options
"""
if sys.argv[1] == 'Minimize':
    mn.main()
elif sys.argv[1] == 'Dimer':
    dm.main()

