from ase import Atoms
from ase.build import bulk
from flapw_python.flapw_calculator import FlapwCalculator

a = 2.95 
c= 4.68
atoms = bulk('Ti', crystalstructure='hcp', a = a, c=c)
# atoms = Atoms()
atoms.calc = FlapwCalculator(a=a, write_lapwin = True)

energy = atoms.get_potential_energy()
print(f"{energy=} eV")