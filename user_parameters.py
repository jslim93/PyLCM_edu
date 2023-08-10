# user_parameters.py
# functions for input and output
import numpy as np
import ipywidgets as widgets
from IPython.display import display, clear_output

# initialization functions: displays widgets and reads in parameters given by the user
def model_steering_input():
    # Create the widgets for the variables
    dt_widget      = widgets.BoundedFloatText(description='dt [s]:', min = 0.0001, max = 5.0, value = 0.5)
    nt_widget      = widgets.IntText(description='nt:', value = 3600)

    Condensation_widget = widgets.Checkbox(description='Condensation:', value=True)
    Collision_widget = widgets.Checkbox(description='Collision:', value=False)

    n_particles_widget = widgets.BoundedIntText(description='n_particles:', min=500, max=1000, value=500)
    n_particles_slider = widgets.IntSlider(description=' ', min=500, max=1000, value=500)
    # link slider and textbox
    mylink = widgets.jslink((n_particles_widget, 'value'), (n_particles_slider, 'value'))

    collision_start_t_widget = widgets.BoundedIntText(description='collision start time [s]:', min = 0,max=10000,value = 0)
    
    max_z_widget = widgets.BoundedFloatText(description='z_max [m]:', min = 0.0, max = 1400.0, step = 0.1, value=1500.0)
    # max value of collision_start_t_widget not updated when nt changed in the textbox! It remains at the default value of nt_widget

    # Display the widgets
    display('Model steering parameters',dt_widget, nt_widget, Condensation_widget,Collision_widget, n_particles_widget, n_particles_slider, \
            collision_start_t_widget, max_z_widget) 
    
    return dt_widget, nt_widget, Condensation_widget, Collision_widget, n_particles_widget, n_particles_slider, collision_start_t_widget, max_z_widget

def parcel_info_input():
    # section for widgets for parcel info
    T_widget = widgets.BoundedFloatText(description='T [K]:', min = 200.0, max = 320.0, step = 0.1, value=293.2)
    P_widget = widgets.BoundedFloatText(description='P [Pa]:', min = 950.0E2, max = 1050.0E2, step = 1, value=1013.0E2)
    RH_widget= widgets.BoundedFloatText(description='RH [-]:', min = 0.01, max = 0.99, step = 0.01, value=0.88)
    w_widget = widgets.BoundedFloatText(description='w [m/s]:', min = 0.0, max = 10, step = 0.1, value=0.5)


    # Display widgets
    display('Parcel initial parameters: ', T_widget, P_widget, RH_widget, w_widget)
    
    return T_widget, P_widget, RH_widget, w_widget