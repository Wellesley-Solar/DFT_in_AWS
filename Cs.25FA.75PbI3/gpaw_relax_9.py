import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
print(currentdir)
parentdir = os.path.dirname(currentdir)
print(parentdir)
sys.path.append(parentdir)

from utilsDFT import *

myGpaw = sqs2QE(False)
myGpaw.load_sqs("./")#Cs.25FA.75PbI3")
myGpaw.writeGPAW()

# print(myGpaw.gcell)
# print(myGpaw.pos)

from gpaw import restart
from ase.parallel import paropen as open
from ase.optimize import QuasiNewton
from gpaw import GPAW, FermiDirac, PW
from ase import Atoms
from ase.visualize import view


atoms = Atoms(symbols=myGpaw.symstr,
             pbc=[ True,  True,  True],
             cell=myGpaw.gcell,
             positions=myGpaw.pos
                    )
# view(atoms)
kpts={'size': (9, 9, 9), 'gamma': True}
calc = GPAW(#h=0.16,
            mode=PW(300),
            kpts=kpts,
            xc='PBE',
            txt=myGpaw.symstr+'_9.out',
            occupations=FermiDirac(width=0.05),
          #   nbands=-2,
            convergence={'energy': 0.0005,  # eV / electron
               'density': 1.0e-4,
               'eigenstates': 4.0e-8,  # eV^2 / electron
               'bands': 'CBM+2.0',
               'forces': float('inf')} # eV / Ang Max
            )

atoms.calc = calc
e2 = atoms.get_potential_energy()
d0 = atoms.get_distance(0, 1)

fd = open('optimization9.txt', 'w')
print('experimental bond length:', file=fd)
print('solid energy: %5.2f eV' % e2, file=fd)
print('bondlength              : %5.2f Ang' % d0, file=fd)

# # Find the theoretical bond length:
# relax = QuasiNewton(atoms, logfile='qn.log')
# relax.run(fmax=0.05)

# e2 = atoms.get_potential_energy()
# d0 = atoms.get_distance(0, 1)

# print(file=fd)
# print('PBE energy minimum:', file=fd)
# print('hydrogen molecule energy: %5.2f eV' % e2, file=fd)
# print('bondlength              : %5.2f Ang' % d0, file=fd)
# fd.close()
