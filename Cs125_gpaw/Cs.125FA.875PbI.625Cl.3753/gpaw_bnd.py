from utilsDFT import *
from gpaw import restart
from ase.parallel import paropen as open
from ase.optimize import QuasiNewton
from gpaw import GPAW, FermiDirac, PW
from ase import Atoms
from ase.visualize import view
import pandas as pd
from ase.io import read

myGpaw = sqs2QE(False, sprcel_dim=[2.0,2.0,2.0])
myGpaw.load_sqs("./", isqs="bestsqs0.out")#Cs.25FA.75PbI3")
myGpaw.writeGPAW()

# label=myGpaw.symstr
# atoms = read(label + '_rerun.gpw')

atoms = Atoms(symbols=myGpaw.symstr,
             pbc=[ True,  True,  True],
             cell=myGpaw.gcell,
             positions=myGpaw.pos
                    )

kpts={'size': (5, 5, 5), 'gamma': True}
calc = GPAW(#h=0.16,
            mode=PW(600),
            kpts=kpts,
            xc='GLLBSC',
            txt=myGpaw.symstr+'band.out',
            occupations=FermiDirac(width=0.05),
          #   nbands=-2,
            convergence={'energy': 0.0005,  # eV / electron
               'density': 1.0e-4,
               'eigenstates': 4.0e-8,  # eV^2 / electron
               'bands': 'CBM+2.0',
               'forces': float('inf')}, # eV / Ang Max
            )

atoms.calc = calc
e2 = atoms.get_potential_energy()
calc.write(myGpaw.symstr+'_band.gpw')

# Band structure calculation with fixed density
bs_calc = calc.fixed_density(#nbands=10,
                             kpts={'path': 'GXMGRX,MR', 'npoints': 60},
                             symmetry='off',
                             convergence={'bands': 'CBM+2.0'},
                             txt=myGpaw.symstr+'_bs.out')

# Plot the band structure
bs = bs_calc.band_structure().subtract_reference()
bs.plot(filename=myGpaw.symstr+'bs.png', emin=-6, emax=6)

# Get the accurate HOMO and LUMO from the band structure calculator
homo, lumo = bs_calc.get_homo_lumo()

# Calculate the discontinuity potential using the ground state calculator and
# the accurate HOMO and LUMO
response = calc.hamiltonian.xc.response
dxc_pot = response.calculate_discontinuity_potential(homo, lumo)

# Calculate the discontinuity using the band structure calculator
bs_response = bs_calc.hamiltonian.xc.response
KS_gap, dxc = bs_response.calculate_discontinuity(dxc_pot)

# Fundamental band gap = Kohn-Sham band gap + derivative discontinuity
QP_gap = KS_gap + dxc

fd = open('band.txt', 'w')
print(file=fd)
print(f'Kohn-Sham band gap:         {KS_gap:.5f} eV', file=fd)
print(f'Discontinuity from GLLB-sc: {dxc:.5f} eV', file=fd)
print(f'Fundamental band gap:       {QP_gap:.5f} eV', file=fd)
fd.close()