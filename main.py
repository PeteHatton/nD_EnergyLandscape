import TwoD_EL.Minimize as mn
import TwoD_EL.Dimer as dm
import sys

if sys.argv[1] == 'Minimize':
    mn.main()
elif sys.argv[1] == 'Dimer':
    dm.main()

