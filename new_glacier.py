"""
script that attempts to calculate the flow of ice in a 2-D grid space
using LandLab. This code has been inspired by the mapping Hedge example 
written by Gregory M. Tucker. This code has also been inspired by the 
component for the thin ice approximation in LandLab.

for more informaiton please go to LandLab.io

@author: ccp127

written by: Cole C. Pazar March 15th, 2016
"""
from __future__ import division
from landlab import RasterModelGrid
import numpy as np

# other libraries that may be needed...

#from landlab.io import read_esri_ascii
#from landlab.plot.imshow import imshow_node_grid
#from matplotlib import pyplot as plt

# define parameters

t_min = 0.0          # inital time [=] years
t_max = 100.0        # final time [=] years
dt = 0.001           # time step [=] years
x_min = 0.0          # minimum distance [=] m
x_max = 10000.0      # maximum distance [=] m
dx = 1000.0          # node spacing [=] m
z_max = 4100.0       # maximum elevation of valley floor [=] m
side_S = 0.2         # slope of the valley sides
S = 0.02             # slope of glacier elevation (linear function)
ELA_mean = 3600.0    # mean equilibrium-line altitude [=] m
gamma = 0.01         # mass balance coefficient [=] m/yr/m
b_cap = 1.5          # cap on the mass balance [=] m/yr
Amp = 500            # Amplitude of the oscillations [=] m
P = 2000             # period of the ELA osillations [=] years
A = (2.1*10**(-16))  # Flow-Law parameter [=] Pa^-3 yr^-1
N_flow = 3           # 3rd order rheology
rho_ice = 917        # density of ice [=] kg m^-3
g = 9.81             # gravity [=] m s^-2
Ho = 10.0

# create a grid manually:

num_rows = 5
num_cols = 7
mg = RasterModelGrid(num_rows, num_cols, dx)
# imshow_node_grid(mg, name='elevation')

"""
# alternativley, you can based it on a DEM with something like this... 
(mg, z) = read_esri_ascii('upper_arkansas_10m.asc', name='elevation')
"""

# create data fields

z = mg.add_empty('node', 'elevation')
H = mg.add_zeros('node', 'ice_thickness')
z_ice = mg.add_empty('node', 'ice_elevation')
dZicedx = mg.add_zeros('link', 'ice_surface_slope')

# initialize 2-D elevations (simple, open book)

z[:] = z_max - mg.node_x * S
z +=  side_S * np.abs(mg.node_y - mg.dx * ((num_rows - 1) / 2.0))

mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
core_nodes = mg.core_nodes

# setting up ice thickness
z_ice[:] = z + H 

class Glacier:

    # ELA
    def ELA_func(self, Amp, P, ELA_mean):
       self.ELA = ELA_mean + (Amp*np.sin(2*np.pi*(self.timestep[t]/P)))
       return self.ELA

    # mass balance
    def b_func (self, gamma, Amp, P, ELA, b_cap):
        self.b = gamma*(self.z-new_glacier.ELA_func(Amp, P, ELA_mean))
        return self.b 

    # loop

    for i in range(2):

        # ice-surface slope
        dZicedx[mg.active_links] = mg.calculate_gradients_at_active_links(z_ice)

        # thickness at links
        H_edge = mg.map_mean_of_link_nodes_to_link('ice_thickness')

        # discharge
        def dQdx_func(self, dx):
            self.dQdx=-np.sign(np.diff(new_glacier.Q(A, rho_ice, g, N_flow, dx))/dx)
            return self.dQdx
        Q = A * ((rho_ice * g * S)**N_flow) 

        # thickness of the ice
        def dHdt_func(self, dt):
            self.dHdt=-(np.diff(new_glacier.H(A, rho_ice, g, N_flow, dt))/dt)
            return self.dHdt

        # testing what was written above...
        print(H)
        print(H_edge)

        # update ice-surface elevation...
        z_ice[:] = z + H

        # more testing...
        print(z_ice)

# run 
def run (self, gamma, Amp, P, ELA_mean, dt, dx):
    global t  # this makes it global to be bale to use in in the ELA_func index
    for t in range(len(self.timestep)):
        self.dQdx=new_glacier.dQdx_func(dx)
        self.dHdt=new_glacier.b_func(gamma, Amp, P, ELA_mean)-self.dqdx
        self.H+=self.dHdt*dt
        self.z=self.H+self.zb
        self.z=self.z
        for i in range(0, len(self.z)):#bottom limit of z can't go below zb
            if self.z[i]<self.zb[i]:
                self.z[i]=self.zb[i]
            if self.timestep[t] % 1 ==0:

# finalize
                def finalize(self):
                    self.fb.savefig('2D_glacier_attempt.jpg')

new_glacier=Glacier()

#initialize
new_glacier.start_conds(t_min, t_max, dt, x_min, x_max, dx, z_max, S, Ho)

#run
new_glacier.run(gamma, Amp, P, ELA, dt, dx)

#finalize
new_glacier.finalize()

"""
end of code
"""
