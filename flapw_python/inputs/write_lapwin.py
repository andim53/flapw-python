from ase.build import bulk
import numpy as np

atom=bulk('Ti','hcp')

with open('lapwin_format','r') as file:
    lines=file.readlines()

with open('lapwin','w') as file:
    a=np.linalg.norm(atom.cell[0])
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