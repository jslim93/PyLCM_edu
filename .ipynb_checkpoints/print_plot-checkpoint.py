import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors
from IPython.display import clear_output
from parameters import *
from micro import *
from chem import *
from aero_init import *
from parcel import *
from condensation import *
from collision import *
from analysis import *

def spec_plot(ax, spectra_arr, nt, log_edges):
    # Clear the previous plot
    clear_output(wait=True)

    # Set the y-scale to logarithmic
    ax.set_yscale("log")

    # Create the contour plot with a logarithmic color scale
    contour = ax.contourf(np.arange(nt+1), log_edges[1:], spectra_arr.T, norm=matplotlib.colors.LogNorm())

    # Create a logarithmic colorbar
    cbar = plt.colorbar(contour, ax=ax, format="%.0e")

    # Add labels and title
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Radius ($\mu m$)")
    ax.set_title("DSD time evolution")
    
    plt.show()
    
    
def print_output(t,dt, z_parcel, T_parcel, q_parcel, qa, qc, qr, na, nc, nr):
    # Clear previous output
    clear_output(wait=True)

    # Print the initial variable names
    print("value: {:<8}  {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}".format(
        "Time (s)", "z (m)", "T (K)", "qv (g/kg)", "QA", "QC (g/kg)", "QR (g/kg)", "NA", "NC", "NR"))

    # Print the updated output
    print("after: {:<8.1f}  {:<8.2f} {:<8.2f} {:<9.2f} {:<8.5f} {:<8.5f} {:<8.5f} {:<8.2f} {:<8.2f} {:<8.2f}".format(
        t * dt + 1, z_parcel, T_parcel, 1e3 * q_parcel, 1e3 * qa, 1e3 * qc, 1e3 * qr, na / 1e6, nc / 1e6, nr / 1e6))