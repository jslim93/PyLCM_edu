# module for animated plots using plotly
# version 1 for 4 plots (testing)
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import time
from IPython.display import clear_output

def animation_init(dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array):
    # sets all initial parameters for the animation
    time_array = np.arange(nt+1)*dt
    
    # Create subplot layout
    fig = make_subplots(rows=1, cols=4)

    # Add initial traces (lines) to subplots
    fig.add_trace(go.Scatter(x=time_array, y=RH_parcel_array, mode='lines'), row=1, col=1)
    fig.add_trace(go.Scatter(x=time_array, y=q_parcel_array, mode='lines'), row=1, col=2)
    fig.add_trace(go.Scatter(x=time_array, y=z_parcel_array, mode='lines'), row=1, col=3)
    fig.add_trace(go.Scatter(x=time_array, y=T_parcel_array, mode='lines'), row=1, col=4)
    
    
    fig.update_layout(height=500, width=1000, title='Continuosly updating RH(t), q_v(t), z(t), T(t) subplot', xaxis_title='Time [s]', yaxis_title='Relative Humidity RH [-]', yaxis_range=[RH_parcel_array[0],1.04], yaxis2_title='q_v', yaxis3_title='Height z [m]', yaxis4_title='Temperature T [K]')
    fig.update_traces(line_color='green')
    
    
    return fig
    
def animation_call(fig, time_array, dt, nt,rm_spec, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, T_parcel_array, RH_parcel_array, q_parcel_array, z_parcel_array):
    # to be included inside the loop
    
    # Clear previous figures
    clear_output(wait=True)
    
    # Update traces (lines) in all subplots
    fig.data[0].x = time_array
    fig.data[0].y = RH_parcel_array
    fig.data[1].x = time_array
    fig.data[1].y = q_parcel_array
    fig.data[2].x = time_array
    fig.data[2].y = z_parcel_array
    fig.data[3].x = time_array
    fig.data[3].y = T_parcel_array
    
    
    # Display updated figure
    fig.show()

    # Pause for a while
    # time.sleep(0.0001)