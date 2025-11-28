import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from flapw_python.flapw_calculator import FlapwCalculator
from scipy.optimize import minimize
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

def build_structure(symbol, structure, a, c=None):
    """
    Build atomic structure for a given crystal type and lattice parameters.

    Parameters:
    - symbol (str): Atomic symbol (e.g., 'Fe')
    - structure (str): Crystal structure type ('fcc', 'bcc', 'hcp', etc.)
    - a (float): Lattice constant a
    - c (float, optional): Lattice constant c (used for 'hcp')

    Returns:
    - ASE Atoms object
    """
    if structure == 'hcp':
        return bulk(symbol, crystalstructure='hcp', a=a, c=c)
    return bulk(symbol, crystalstructure=structure, a=a)

def get_energy(symbol, structure, a, c=None, lapwin_format=None, lapwout=None, command='./run_flapw.sh'):
    """
    Calculate total potential energy for a structure using FLAPW.

    Parameters:
    - symbol (str): Atomic symbol
    - structure (str): Crystal structure
    - a, c (float): Lattice constants
    - lapwin_format, lapwout (str): Paths for FLAPW input/output files

    Returns:
    - Total energy (float)
    """
    atoms = build_structure(symbol, structure, a, c)
    atoms.calc = FlapwCalculator(lapwin_format=lapwin_format, lapwout=lapwout, a=a, command=command)
    return atoms.get_potential_energy()

def optimize_lattice(symbol, structure, a0, c0=None, da=0.05, dc=0.1, n=5,
                     hcp_mode='a_only',
                     lapwin_format='./NECESSITY/lapwin_format', lapwout='./lapwoutSCF', command='./run_flapw.sh'):
    """
    Optimize lattice constant(s) by fitting a quadratic to energy evaluations.

    Parameters:
    - symbol (str): Atomic symbol
    - structure (str): Crystal structure ('fcc', 'bcc', 'hcp')
    - a0, c0 (float): Initial guess for a and c
    - da, dc (float): Perturbation step size for a and c
    - n (int): Number of sampling points
    - hcp_mode (str): 'a_only', 'c_only', or 'both' for hcp mode
    - lapwin_format, lapwout (str): FLAPW input/output paths

    Returns:
    - Tuple with optimized constants and energies sampled
    """
    if structure in ['fcc', 'bcc']:
        a_list = np.linspace(a0 - da, a0 + da, n)
        energies = []

        for a in a_list:
            try:
                e = get_energy(symbol, structure, a, lapwin_format=lapwin_format, lapwout=lapwout,command=command)
                energies.append(e)
            except Exception as err:
                print(f"Skipped a={a:.3f}: {err}")

        # Fit quadratic curve to energy vs. a
        coeffs = np.polyfit(a_list, energies, 2)
        a_opt = -coeffs[1] / (2 * coeffs[0])
        e_opt = np.polyval(coeffs, a_opt)
        print(f'{structure.upper()} optimized: a = {a_opt:.4f} Å, energy = {e_opt:.6f} eV')

        return a_opt, e_opt, a_list, energies

    elif structure == 'hcp':
        if hcp_mode == 'a_only':
            a_list = np.linspace(a0 - da, a0 + da, n)
            energies = []

            for a in a_list:
                try:
                    e = get_energy(symbol, 'hcp', a, c0, lapwin_format=lapwin_format, lapwout=lapwout,command=command)
                    energies.append(e)
                except Exception as err:
                    print(f"Skipped a={a:.3f}: {err}")
                    energies.append(np.nan)

            a_list = np.array(a_list)
            energies = np.array(energies)
            valid = ~np.isnan(energies)

            coeffs = np.polyfit(a_list[valid], energies[valid], 2)
            a_opt = -coeffs[1] / (2 * coeffs[0])
            e_opt = np.polyval(coeffs, a_opt)

            print(f'HCP optimized a (c={c0} fixed): a = {a_opt:.4f} Å, energy = {e_opt:.6f} eV')
            return a_opt, c0, e_opt, a_list, energies

        elif hcp_mode == 'c_only':
            c_list = np.linspace(c0 - dc, c0 + dc, n)
            energies = []

            for c in c_list:
                try:
                    e = get_energy(symbol, 'hcp', a0, c, lapwin_format=lapwin_format, lapwout=lapwout,command=command)
                    energies.append(e)
                except Exception as err:
                    print(f"Skipped c={c:.3f}: {err}")
                    energies.append(np.nan)

            c_list = np.array(c_list)
            energies = np.array(energies)
            valid = ~np.isnan(energies)

            coeffs = np.polyfit(c_list[valid], energies[valid], 2)
            c_opt = -coeffs[1] / (2 * coeffs[0])
            e_opt = np.polyval(coeffs, c_opt)

            print(f'HCP optimized c (a={a0} fixed): c = {c_opt:.4f} Å, energy = {e_opt:.6f} eV')
            return a0, c_opt, e_opt, c_list, energies

        else:
            raise ValueError(f"Unknown hcp_mode '{hcp_mode}': choose 'a_only', 'c_only', or 'both'")

def plot_optimization_result(structure, a_or_c_list=None, energies=None, a_opt=None, c_opt=None, e_opt=None,
                             output_file='fig_lattice_opt.png', figsize=(2, 3)):
    """
    Plot energy vs lattice constant (a or c) with consistent styling and save as PNG.

    Parameters:
    - structure (str): Crystal structure name for title
    - a_or_c_list (array): List of a or c values
    - energies (array): Corresponding total energies
    - a_opt (float): Optimal lattice constant
    - e_opt (float): Corresponding minimal energy
    - output_file (str): Filename to save the figure
    - figsize (tuple): Figure size in inches (width, height)
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Convert to numpy arrays
    a_or_c_array = np.array(a_or_c_list)
    energies_array = np.array(energies)

    # Filter out NaNs
    mask = ~np.isnan(energies_array)
    a_or_c_array = a_or_c_array[mask]
    energies_array = energies_array[mask]

    # Set plot style
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # Define axis labels
    xlabel_plot = 'Lattice constant (Å)'
    ylabel_plot = 'Total energy (eV)'

    # Define plotting limits based on data range
    xmin_plot = min(a_or_c_array)
    xmax_plot = max(a_or_c_array)
    ymin_plot = min(energies_array)
    ymax_plot = max(energies_array)

    # Add some margin to the limits
    x_margin = (xmax_plot - xmin_plot) * 0.1
    y_margin = (ymax_plot - ymin_plot) * 0.1

    # Create plot
    fig = plt.figure(dpi=150, figsize=figsize)
    plt.plot(a_or_c_array, energies_array, 'o-', color='black', linewidth=0.8, alpha=1.0, label='Total energy')

    # Mark optimal point
    if a_opt is not None:
        plt.plot(a_opt, e_opt, 'o', color='blue', markersize=3, label=f'Minimum at {a_opt:.3f} Å')
    
    if c_opt is not None:
        plt.plot(c_opt, e_opt, 'o', color='blue', markersize=3, label=f'Minimum at {c_opt:.3f} Å')

    # Set limits and ticks
    plt.xlim(xmin_plot - x_margin, xmax_plot + x_margin)
    plt.ylim(ymin_plot - y_margin, ymax_plot + y_margin)
    plt.xticks(np.linspace(xmin_plot, xmax_plot, 5), fontsize=10)
    plt.yticks(np.linspace(ymin_plot, ymax_plot, 5), fontsize=10)

    # Axis labels and title
    plt.xlabel(xlabel_plot, fontsize=10)
    plt.ylabel(ylabel_plot, fontsize=10)
    plt.title(f'{structure.upper()} Lattice Optimization', fontsize=10)

    # Save figure
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_file, format='png', bbox_inches='tight', transparent=True)
    plt.show()
    plt.clf()
    plt.close()



# from flapw_python import optimize_lattice as opt

# symbol = 'Zr'
# crystal = 'hcp'
# a = 3.23 
# c = 5.15 

# lapwin_format = './NECESSITY/lapwin_format'
# lapwout = './lapwoutSCF'

# # Only optimize a, keep c fixed
# a_opt, c0, e_opt, a_list, energies = opt.optimize_lattice(symbol, crystal, a0=a, c0=c, n=3,hcp_mode='c_only')
# opt.plot_optimization_result(crystal, a_or_c_list=a_list, energies=energies, a_opt=a_opt, e_opt=e_opt, figsize=(6,3))
