import os
import numpy as np
from ase.calculators.calculator import FileIOCalculator

class FlapwCalculator(FileIOCalculator):
    implemented_properties = ['energy']
    def __init__(self, 
                 restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 atoms=None,
                 lapwin_format=None, 
                 lapwout=None,
                 a=None,
                 command='./run_flapw.sh', 
                 label='flapw', **kwargs):
        self.atoms = None
        self.lapwin_format=lapwin_format
        self.lapwout=lapwout
        self.a = a
        
        # self.command = command 

        FileIOCalculator.__init__(self, 
                                  restart=restart, 
                                  ignore_bad_restart_file=ignore_bad_restart_file, 
                                  label=label, 
                                  atoms=atoms, 
                                  command=command,
                                  **kwargs)
       
    def write_input(self, atoms, properties=None, system_changes=None):
        atom = atoms
        FileIOCalculator.write_input(self, atoms, properties,system_changes)
        with open(self.lapwin_format,'r') as file:
            lines=file.readlines()

        with open('lapwin','w') as file:
            if self.a == None:
                a=np.linalg.norm(atom.cell[0])
            else:
                a=self.a
            w_pos=False
            for i in range(len(lines)):
                if i==0:
                    line=lines[i].split(':')
                    file.write(f"{line[0]}: {atom.get_chemical_formula()}\n")
                    continue
                if i==3:
                    file.write(f"{a}/0.529177\n")
                    continue
                if i==4:
                    line=''.join(f' {x}   ' for x in atom.cell[0]/a)
                    file.write(line+'\n')
                    continue
                if i==5:
                    line=''.join(f' {x}   ' for x in atom.cell[1]/a)
                    file.write(line+'\n')
                    continue
                if i==6:
                    line=''.join(f' {x}   ' for x in atom.cell[2]/a)
                    file.write(line+'\n')
                    continue
                if i==10:
                    for label,row in zip(atom.get_chemical_symbols(),atom.get_scaled_positions()):
                        line=label+'  '+'          '.join(str(x) for x in row)
                        file.write(line+'\n')
                else:
                    file.write(lines[i])

    def read_results(self):
        print("Read Results")
        hartree_to_ev = 27.211386245988
        with open(self.lapwout, 'r') as file:
            for line in file:
                if 'total energy for' in line:
                    parts = line.strip().split()
                    try:
                        energy_htr = float(parts[-2])
                        self.results['energy'] = energy_htr * hartree_to_ev
                        break
                    except (IndexError, ValueError):
                        raise RuntimeError('Failed to parse total energy from lapwout')


# from ase.build import bulk
# a=2.866
# atoms = bulk('Fe',crystalstructure='bcc',a=a)
# atoms.calc = FlapwCalculator(a=a)

# energy = atoms.get_potential_energy()
# print(f"{energy=} eV")