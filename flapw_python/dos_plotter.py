import os
import numpy as np
import matplotlib.pyplot as plt

def read_dos(file_path, x_column, y_column):
    """Reads DOS data from file and returns x and y arrays."""
    x_value, y_value = [], []
    with open(file_path, 'r') as file:
        for line in file:
            data = line.strip().split()
            if data and data[0] != '#':
                if len(data) >= max(x_column, y_column):
                    x_value.append(float(data[x_column - 1]))
                    y_value.append(float(data[y_column - 1]))
    return x_value, y_value

def plot_dos(file_path='dos.xy', number_of_atoms=1, output_file='fig_dos.png', figsize=(1.5, 3)):
    """Plots the Density of States from a .xy file and saves it as a PNG."""
    # Set plot style
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['mathtext.fontset'] = 'cm'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    # Check if spin-polarized by reading header
    with open(file_path, 'r') as file:
        header = file.readline().strip()
        columns = header.split()
        lgc_spin = len(columns) == 8

    # Set plot limits and labels
    xlabel_plot = "DOS (/eV-atom)"
    ylabel_plot = "Energy (eV)"
    xmin_plot, xmax_plot, xstp_plot = 0.0, 5.0, 1.0
    ymin_plot, ymax_plot, ystp_plot = -6, 5.0, 2.0

    # Set columns and styles
    xcol_data, ycol_data = [], []
    color_data, style_data = [], []

    if lgc_spin:
        # Spin-up
        xcol_data.append(2)
        ycol_data.append(1)
        color_data.append('blue')
        style_data.append('-')
        # Spin-down
        xcol_data.append(3)
        ycol_data.append(1)
        color_data.append('red')
        style_data.append('-')
    else:
        # Non-spin-polarized
        xcol_data.append(2)
        ycol_data.append(1)
        color_data.append('black')
        style_data.append('-')

    # Read and normalize data
    x_plot, y_plot = [], []
    for i in range(len(xcol_data)):
        x_data, y_data = read_dos(file_path, xcol_data[i], ycol_data[i])
        x_plot.append(np.array(x_data) / number_of_atoms)
        y_plot.append(y_data)

    print('xmin =', np.min(x_plot), 'xmax =', np.max(x_plot))
    print('ymin =', np.min(y_plot), 'ymax =', np.max(y_plot))

    # Create figure and plot
    fig = plt.figure(dpi=150, figsize=figsize)
    for i in range(len(x_plot)):
        x = y_plot[i]
        y = x_plot[i]
        plt.plot(x, y, color=color_data[i], linestyle=style_data[i], linewidth=0.5, alpha=1.0)

    # Set limits and labels
    plt.xlim(ymin_plot, ymax_plot)
    plt.ylim(xmin_plot, xmax_plot)
    plt.yticks(np.arange(xmin_plot, xmax_plot + xstp_plot, xstp_plot), fontsize=10)
    plt.xticks(np.arange(ymin_plot, ymax_plot + ystp_plot, ystp_plot), fontsize=10)
    plt.ylabel(xlabel_plot, fontsize=10)
    plt.xlabel(ylabel_plot, fontsize=10)

    # Add zero lines
    plt.axvline(0, color='black', linestyle='--', linewidth=0.5)
    plt.axhline(0, color='black', linestyle='--', linewidth=0.5)

    # Save and close
    plt.savefig(output_file, format='png', bbox_inches='tight', transparent=True)
    plt.show()
    plt.clf()
    plt.close()
