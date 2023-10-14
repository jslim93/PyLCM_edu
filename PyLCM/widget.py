# widget.py
# Functions for input and output (e.g. widgets settings, file output settings)
import numpy as np
import ipywidgets as widgets
from IPython.display import display, clear_output

# initialization functions: display widgets and read in parameters given by the user
def model_steering_input():
    # Adjust the style that descriptions are displayed in full length
    style = {'description_width': 'initial'}
    
    # Create the widgets for the variables
    dt_widget      = widgets.BoundedFloatText(description='dt (s):', min = 0.0001, max = 5.0, value = 0.5, style=style)
    nt_widget      = widgets.IntText(description='nt:', value = 3600, style=style)

    Condensation_widget = widgets.Checkbox(description='Condensation:', value=True, style=style)
    Collision_widget = widgets.Checkbox(description='Collision:', value=False, style=style)

    n_particles_widget = widgets.BoundedIntText(description='n_particles:', min=100, max=1000, value=500, style=style)
   
    max_z_widget = widgets.BoundedFloatText(description='z_max (m):', min = 0.0, max = 2000.0, step = 0.1, value=1500.0, style=style)
    z_widget = widgets.BoundedFloatText(description='z_0 (m):', min = 0.0, max = 800.0, step = 0.1, value=0.0, style=style)

    # Display the widgets
    display('Model steering parameters',dt_widget, nt_widget, Condensation_widget,Collision_widget, n_particles_widget, max_z_widget) 
    
    return dt_widget, nt_widget, Condensation_widget, Collision_widget, n_particles_widget, max_z_widget

def parcel_info_input():
    style = {'description_width': 'initial'}
    # Section for widgets for parcel info
    T_widget = widgets.BoundedFloatText(description='T (K):', min = 200.0, max = 320.0, step = 0.1, value=293.2, style=style)
    P_widget = widgets.BoundedFloatText(description='P (Pa):', min = 950.0E2, max = 1050.0E2, step = 1, value=1013.0E2, style=style)
    RH_widget= widgets.BoundedFloatText(description='RH (-):', min = 0.01, max = 0.99, step = 0.01, value=0.88, style=style)
    w_widget = widgets.BoundedFloatText(description='w (m/s):', min = 0.0, max = 10, step = 0.1, value=0.5, style=style)
    z_widget = widgets.BoundedFloatText(description='z (m):', min = 0.0, max = 20000.0,  step = 100.0, value=0.0, style=style)


    # Display widgets
    display('Parcel initial parameters: ', T_widget, P_widget, RH_widget, w_widget, z_widget)
    
    return T_widget, P_widget, RH_widget, w_widget, z_widget

def ascending_mode_input():
    # Users can choice the ascending mode with this widget.
    ascending_mode_widget = widgets.ToggleButtons(options=['linear', 'sine', 'in_cloud_oscillation'], value='linear', description='Mode', layout={'width': 'max-content'}, disabled=False)
    display(ascending_mode_widget)
    
    return ascending_mode_widget

def aero_mode_input():
    # Widgets for aerosol initialisation
    mode_aero_init_widget = widgets.ToggleButtons(options=['weighting_factor', 'random'], value='weighting_factor', layout={'width': 'max-content'}, disabled=False)
    print("\n")
    # Display widgets
    display('Aerosol initialisation mode: ', mode_aero_init_widget)
    
    return mode_aero_init_widget

# Initialization function for the grid widget for aerosol modes
def grid_modes_input():
    # N, mu, sigma for each mode
    print('N_aero: number of aerosols per cubic centimeter, mu: mean droplet radius, sigma: std of mu')

    # Arrange the widgets in a 5x4 layout
    gridwidget = widgets.GridspecLayout(5, 4)
    
    # Adjust the style that descriptions are displayed in full length
    style = {'description_width': 'initial'}

    # Column 0 for mode 1
    gridwidget[0, 0]= widgets.Button(description='Mode #1', disabled=False, button_style='info', tooltip='Heading', style=style)
    gridwidget[1, 0] = widgets.BoundedFloatText(description='N_aero (cm⁻3)', min = 0.0, max = 5000.0, step = 0.1, value=1000.0, style=style)
    gridwidget[2, 0] = widgets.BoundedFloatText(description='mu (µm)', min = 0.0, max = 5.0, step = 0.001, value=0.008, style=style)
    gridwidget[3, 0] = widgets.BoundedFloatText(description='sigma (-)', min = 0.0, max = 3.0, step = 0.1, value=1.6, style=style)
    gridwidget[4, 0] = widgets.BoundedFloatText(description='Hygroscopicity parameter', min = 0.0, max = 100.0, step = 0.1, value=1.6, style=style)
    
    # Column 1 for mode 2
    gridwidget[0, 1]= widgets.Button(description='Mode #2', disabled=False, button_style='info', tooltip='Heading')
    gridwidget[1, 1] = widgets.BoundedFloatText(description='N_aero (cm⁻3)', min = 0.0, max = 5000.0, step = 0.1, value=800.0, style=style)
    gridwidget[2, 1] = widgets.BoundedFloatText(description='mu (µm)', min = 0.0, max = 5.0, step = 0.001, value=0.034, style=style)
    gridwidget[3, 1] = widgets.BoundedFloatText(description='sigma (-)', min = 0.0, max = 3.0, step = 0.1, value=2.1, style=style)
    gridwidget[4, 1] = widgets.BoundedFloatText(description='Hygroscopicity parameter', min = 0.0, max = 100.0, step = 0.1, value=1.6, style=style)
    # Column 2 for mode 3
    gridwidget[0, 2]= widgets.Button(description='Mode #3', disabled=False, button_style='info', tooltip='Heading')
    gridwidget[1, 2] = widgets.BoundedFloatText(description='N_aero (cm⁻3)', min = 0.0, max = 5000.0, step = 0.01, value=0.72, style=style)
    gridwidget[2, 2] = widgets.BoundedFloatText(description='mu (µm)', min = 0.0, max = 5.0, step = 0.1, value=0.46, style=style)
    gridwidget[3, 2] = widgets.BoundedFloatText(description='sigma (-)', min = 0.0, max = 3.0, step = 0.1, value=2.2, style=style)
    gridwidget[4, 2] = widgets.BoundedFloatText(description='Hygroscopicity parameter', min = 0.0, max = 100.0, step = 0.1, value=1.6, style=style)
    # Column 3 for mode 4
    gridwidget[0, 3]= widgets.Button(description='Mode #4', disabled=False, button_style='info', tooltip='Heading')
    gridwidget[1, 3] = widgets.BoundedFloatText(description='N_aero (cm⁻3)', min = 0.0, max = 5000.0, step = 0.1, value=0, style=style)
    gridwidget[2, 3] = widgets.BoundedFloatText(description='mu (µm)', min = 0.0, max = 5.0, step = 0.1, value=0, style=style)
    gridwidget[3, 3] = widgets.BoundedFloatText(description='sigma (-)', min = 0.0, max = 3.0, step = 0.1, value=0, style=style)
    gridwidget[4, 3] = widgets.BoundedFloatText(description='Hygroscopicity parameter', min = 0.0, max = 100.0, step = 0.1, value=1.6, style=style)
    # An additional button will be added later below of these widgets for the plot of the cumulative distribution.
    # plotbutton_widget = widgets.Button(description='Plot cumulative distr.', disabled=False, button_style='info', tooltip='Plot cumulative distribution of all modes')

    display(gridwidget)
    # display(plotbutton_widget)
    
    return gridwidget

# Function for kohler activation radius widget
def kohler_settings():
    # Adjust the style that descriptions are displayed in full length
    style = {'description_width': 'initial'}
    # Set if activation radius uses kohler critical radius
    kohler_widget = widgets.Checkbox(description='Koehler critical radius', value=False, style=style, layout={'width': 'max-content'})
    
    display(r"Set activation radius using Koehler critical radius, otherwise activation radius is 1 µm",kohler_widget)
     
    return kohler_widget

# Function for kohler activation radius widget
def kappa_settings():
    # Adjust the style that descriptions are displayed in full length
    style = {'description_width': 'initial'}
    # Set if activation radius uses kohler critical radius
    kappa_widget = widgets.Checkbox(description='Hygroscopicity', value=False, style=style, layout={'width': 'max-content'})

    display(r"Use hygroscopicity parameter from kappa-kohler theory",kappa_widget)

    return kappa_widget

# Function for time step widgets
def timestep_display_mode_settings():
    # Widget to set the mode of the continuously updated output during model run time: either output of variables in text form (fast, recommended) = 'text'
    # or: continuously updated plot using plotly (slower, memory consuming) = 'graphics'
    mode_displaytype_widget = widgets.ToggleButtons(options=['text_fast', 'graphics'], value='graphics', layout={'width': 'max-content'}, disabled=False)
    display('Display mode: ', mode_displaytype_widget)
    
    return mode_displaytype_widget

# Function for plot widgets
def plot_widgets_settings(nt):
    # Users can choice if they want the plots time-dependent or height-dependent (except of DSD and particle densities always being time dependent)
    # Adjust the style that descriptions are displayed in full length
    style = {'description_width': 'initial'}
    mode_plots_widget = widgets.ToggleButtons(options=['time-series', 'vertical profile'], value='time-series', description='Plots are:', layout={'width': 'max-content'}, disabled=False, style=style)
    # Widget to set in which increment lines for the spectra are drawn
    increment_widget = widgets.BoundedIntText(description='For droplet spectra: increment: ', value=20, min = 1, max = nt, step=1, style=style)
    
    # Display the widgets
    display(mode_plots_widget, 'For droplet spectra: set increment (e.g. 20: every 20th timestep a line is drawn): ', increment_widget)
    
    return mode_plots_widget, increment_widget