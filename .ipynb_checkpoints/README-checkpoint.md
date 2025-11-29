# Introduction

We developed **FLAPW-PYTHON Wrapper v0.1** to provide a Python interface for running FLAPW commands.

In this demo, we use script to optimize the lattice parameters and calculate the Density of States (DOS) for bulk materials with fcc/bcc, and hcp crystal structures.

First, clone the repository from GitHub:

```
git clone https://github.com/andim53/flapw-python.git
```

The repository includes all necessary components, such as input files (`infiles`), binaries (`bin` and `src`, containing FLAPW executables for each operating system), and the `flapw_python` module, which provides Python scripts for running the FLAPW calculator, plotting DOS, and optimizing lattice parameters.

Inside the `flapw_python` directory, bash scripts are included to execute the underlying FLAPW binaries. Two scripts are provided: one for running the self-consistent field (SCF) calculation, and another for the DOS calculation.

The SCF bash script (`run_flapw.sh`) includes the following commands:

```bash
export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
./src_macos/flapw
./bin/FLcopy SCF
./bin/FLclean
cp lapwinSCF lapwin
```

The DOS bash script (`run_flapw_dos.sh`) contains:

```bash
export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib:$DYLD_LIBRARY_PATH
./bin/FLrst SCF
sed -i '' 's/ 12  12  12 / 33  33  33 /' lapwin
./src_macos/flapw
./bin/FLcopy SCF
./bin/FLclean
cp lapwinSCF lapwin
sed -i '' 's/Density of states:F/Density of states:T/' lapwin
cp ./infiles/dosin .
./bin/FLrst SCF
./src_macos/flapw
./src_macos/xdos
mkdir DOS/
mv $(find . -maxdepth 1 -type f -name '*in') DOS/
mv $(find . -maxdepth 1 -type f -name '*out') DOS/
mv $(find . -maxdepth 1 -type f -name '*xy') DOS/
./bin/FLclean
```

The SCF script runs the SCF calculation, cleans up temporary files, and keeps only the \*SCF file (FLAPW output).

The DOS script on the other hand uses the existing \*SCF file, updates the k-point to larger grid, enables DOS in the `lapwin` file, and performs DOS calculation and organize the DOS output.

# Python Scripts

To run the Python script, we provide an example Jupyter Notebook, `main.ipynb`. We provide a set of example materials along with their experimental lattice parameters. 

Below is the lattice data:

```python
lattice_data = {
    'Ti': {'structure': 'hcp', 'a': 2.95, 'c': 4.68},  # Titanium, hexagonal close-packed
    'Cr': {'structure': 'bcc', 'a': 2.88},             # Chromium, body-centered cubic
    'Ni': {'structure': 'fcc', 'a': 3.52},             # Nickel, face-centered cubic
    'Zr': {'structure': 'hcp', 'a': 3.23, 'c': 5.15},  # Zirconium, hexagonal close-packed
    'Mo': {'structure': 'bcc', 'a': 3.15},             # Molybdenum, body-centered cubic
    'Pd': {'structure': 'fcc', 'a': 3.89},             # Palladium, face-centered cubic
    'Hf': {'structure': 'hcp', 'a': 3.19, 'c': 5.05},  # Hafnium, hexagonal close-packed
    'W':  {'structure': 'bcc', 'a': 3.16},             # Tungsten, body-centered cubic
    'Pt': {'structure': 'fcc', 'a': 3.92}              # Platinum, face-centered cubic
}
```

Our script runs using the ASE (Atomic Simulation Environment) package, so make sure it is installed:

```bash
pip install ase
```

Before running the notebook, ensure that the following directories and files in the same directory as main.ipynb:

- `flapw_python/` (Python module for FLAPW)
- `infiles/` (contains input templates like `lapwinNOSPIN`)
- `bin/` (compiled FLAPW utilities)
- `src_*/` (FLAPW binaries)

The first step is to import the necessary modules from ASE and FLAPW-PYTHON.

```python
# Import the bulk structure builder from ASE
from ase.build import bulk

# Import the FLAPW calculator wrapper
from flapw_python.flapw_calculator import FlapwCalculator

# Define the LAPW input template to be used (non-spin-polarized)
lapwin_format = './infiles/lapwinNOSPIN'

# Specify the name of the SCF output file
lapwout = './lapwoutSCF'

# Define the path to the shell script that will execute FLAPW commands
command = './flapw_python/run_flapw.sh'
```

In this setup:

- `lapwin_format` points to the LAPW input template (e.g. without spin polarization),
- `lapwout` is where the SCF output will be saved,
- `command` is the shell script that manages the FLAPW execution (SCF or DOS).

Note: Before running the scripts, make sure to grant executable permission to the bash files. Also, keep in mind that executing bash commands on Windows will differ, so adjustments may be needed to ensure compatibility.  

```bash
chmod +x flapw_python/run*
```

# Running the Python Scripts

For the fcc/bcc demo, we use palladium (Pd) and by optimizing the lattice constant `a`.

```python
# Select the element and retrieve its structure and initial lattice constant
symbol = 'Pd'
data = lattice_data[symbol]
crystal = data['structure']
a = data['a']  # Initial lattice parameter (experimental)
```

The following script performs the lattice optimization:

```python
from flapw_python import optimize_lattice as opt

# Optimize lattice constant 'a' for fcc/bcc structure
a_opt, e_opt, a_list, energies = opt.optimize_lattice(
    symbol=symbol, 
    structure=crystal, 
    a0=a,           # Initial guess for lattice constant
    da=0.1,         # Search range: ±0.1 Å from a0
    n=4,            # Number of points to sample
    lapwin_format=lapwin_format,
    lapwout=lapwout,
    command=command
)

# Plot the energy vs. lattice constant curve and save the figure
opt.plot_optimization_result(
    structure=crystal, 
    a_or_c_list=a_list, 
    energies=energies, 
    a_opt=a_opt, 
    e_opt=e_opt, 
    figsize=(4, 3),
    output_file=f'./attachments/fig_latt_opt_{symbol}.png'
)
```

In this scripts:

- We define `symbol` as Pd and set `a0 = 3.89 Å` (from experimental data).
- `da = 0.1` specifies the variation range to explore around `a0` (±0.1 Å).
- `n = 4` sets the number of evaluation points.
- `a_opt` and `e_opt` are the optimized lattice constant and corresponding minimum energy, obtained via quadratic fitting.
- `a_list` and `energies` contain the scanned lattice constants and their respective total energies. 

This figure shows optimized lattice:

![[./attachments/fig_latt_opt_Pd.png]]

Note: If the minimum is not well-converged (e.g., if the energy curve appears too flat or the minimum lies near the edge of the scanned range), we can rerun the optimization using the latest lattice (`a_opt`) as the new starting point (`a0`). Additionally, to perform a finer optimization, we reduce the step size (`da`) to 0.05 Å while keeping `a_opt` as the `a0`.

Once we have obtained the optimized lattice constant `a_opt`, we can proceed with the DOS calculation.

```python 
# Perform DOS calculation for fcc/bcc system
# Use the bash script specifically for DOS generation
command_dos = './flapw_python/run_flapw_dos.sh'

# Build the bulk structure using the optimized lattice constant
atoms = bulk(symbol, crystalstructure=crystal, a=a_opt)

# Assign the FLAPW calculator with DOS-specific command
atoms.calc = FlapwCalculator(
    lapwin_format=lapwin_format,
    lapwout=lapwout,
    a=a_opt,
    command=command_dos
)

# Run the calculation and retrieve the total energy
energy = atoms.get_potential_energy()
print(energy)
```

In this script:

- We switch to `run_flapw_dos.sh` to enable the full DOS workflow.
- The structure is built using `a_opt`, the optimized lattice constant.
- The assigned calculator uses the updated command and runs the necessary steps for DOS.
- Finally, we compute the total energy, which also ensures that the DOS-related files are generated.

Once calculation is finished, we run the DOS plotting script:


```python
from flapw_python import dos_plotter

output_file = './attachments/'

# Generate and save the DOS plot for fcc/bcc system
dos_plotter.plot_dos(
    file_path='DOS/dos.xy',              # Path to the DOS data
    number_of_atoms=len(atoms),          # Normalize DOS by number of atoms
    figsize=(3, 3),                       # Set figure size
    output_file=output_file + f'fig_dos_{symbol}.png'
)
```


Below is the resulting Density of States (DOS) plot for Pd:
![[./attachments/fig_dos_Pd.png]]

# Calculating the HCP structure

Unlike bcc or fcc structures, hcp requires optimization of both lattice constant `a` and `c`. In this demo, we perform the calculation for titanium (Ti).

```python
# Select hcp element and extract lattice parameters
symbol = 'Ti'
data = lattice_data[symbol]
crystal = data['structure']
a = data['a']
c = data['c']
```

We start by optimizing the lattice constant `a`, while keeping `c` fixed.

```python
# Optimize lattice constant 'a' for hcp (with fixed c)
a_opt, c0, e_opt, a_list, energies = opt.optimize_lattice(
    symbol=symbol, 
    structure=crystal, 
    a0=a,
    c0=c,
    da=0.1,        # ±0.1 Å range around initial a
    n=4,           # Number of points to sample in each direction
    hcp_mode='a_only',  # Enable a-only optimization mode for hcp
    lapwin_format=lapwin_format,
    lapwout=lapwout,
    command=command
)

# Plot and save the result
opt.plot_optimization_result(
    structure=crystal, 
    a_or_c_list=a_list, 
    energies=energies, 
    a_opt=a_opt, 
    e_opt=e_opt, 
    figsize=(4, 3),
    output_file=f'./attachments/fig_latt_opt_{symbol}_aOnly.png'
)
```

In this script, we enable `hcp_mode='a_only'` to vary `a` while keeping `c` constant. This is the first step in a two-stage optimization process for hcp systems.

Below is the result of the `a-only` optimization:

![[./attachments/fig_latt_opt_Ti_aOnly.png]]

After obtaining the optimized lattice constant `a_opt`, we now use it to optimize the `c` parameter:

```python
# Optimize lattice constant 'c' for hcp (with fixed a)
a0, c_opt, e_opt, c_list, energies = opt.optimize_lattice(
    symbol=symbol, 
    structure=crystal, 
    a0=a_opt,
    c0=c,
    da=0.1,            # ±0.1 Å range around initial c
    n=4,               # Number of points to sample in each direction
    hcp_mode='c_only', # Enable c-only optimization mode for hcp
    lapwin_format=lapwin_format,
    lapwout=lapwout,
    command=command
)

# Plot and save the result
opt.plot_optimization_result(
    structure=crystal, 
    a_or_c_list=c_list, 
    energies=energies, 
    a_opt=c_opt, 
    e_opt=e_opt, 
    figsize=(4, 3),
    output_file=f'./attachments/fig_latt_opt_{symbol}_cOnly.png'
)
```

The figure below shows the result of the `c-only` optimization:

![[./attachments/fig_latt_opt_Ti_cOnly.png]]

With both `a_opt` and `c_opt` determined, we can now proceed with the DOS calculation and plot the DOS:

```python
# Calculate DOS for hcp system
command_dos = './flapw_python/run_flapw_dos.sh'

# Build the bulk hcp structure with optimized a and c
atoms = bulk(symbol, crystalstructure=crystal, a=a_opt, c=c_opt)

# Assign the FLAPW calculator with DOS-specific command
atoms.calc = FlapwCalculator(
    lapwin_format=lapwin_format,
    lapwout=lapwout,
    a=a_opt,
    command=command_dos
)

# Run the calculation and get the total energy
energy = atoms.get_potential_energy()
print(energy)
```

```python
from flapw_python import dos_plotter

output_file = './attachments/'

# Plot and save the DOS figure
dos_plotter.plot_dos(
    file_path='DOS/dos.xy', 
    number_of_atoms=len(atoms), 
    figsize=(3, 3),
    output_file=output_file + f'fig_dos_{symbol}.png'
)
```

Below is the resulting DOS for Ti:

![[./attachments/fig_dos_Ti.png]]

