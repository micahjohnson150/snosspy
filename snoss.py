 t#Work flow.


#1. Run smrf over the area of interest to get forcing data
#   * Precip
#   * Air temp
#   * Snow density

#2. SETUP MODELING
# * Create empty data sets to work in
# * Calculate slope from DEM
# * Establish a max depth domain and discretize


# 3. Enter SNOSS Loop
    # a. calculate inital depth from precip and density add this depth to depth domain
    # for d in Depth_domain
    #   calculate new density in time at each point
    #

import numpy as np
from netCDF4 import Dataset
import os
from matplotlib import pyplot as plt


class SNOSS(object):

    def __init__(self,angle,air_temp,precip, new_density, dt = 3600):
        """
        angle - Slope angle
        air_temp - air_temp time series numpy array
        precip - precip in mm time series
        snow_density - initial snow density as it fell (time series)
        dt - hours change in time for modeling
        """
        self.time_range = [dt*i for i in range(len(precip))]
        self.dt = dt

        self.precip = precip
        self.air_temp = air_temp
        self.new_density = new_density

        self.rho_ice = 917.0 # Density of ice
        self.g = 9.81 # Gravity
        self.b1 = 6.5*10**-7 # Viscosity constant
        self.b2 = 19.3 # Viscosity constant
        self.E = 67.3 # Activation energy
        self.R = 0.0083 # Gas constant

        self.angle = angle
        self.surface_layer_count = 0

        #Identify layers to be tracked
        layers = 0
        for r in (precip > 0.0):
            if r:
                layers+=1
        print('Tracking {0} layers'.format(layers))


        self.rho = np.zeros((len(self.precip),layers))
        self.deposition = np.zeros(layers)
        self.basal_temp = np.zeros(layers)

        # Time savers
        self.g_component = self.g*(np.cos(self.angle))**2.0

    def solve(self):
        # time domain
        for tstep,t in enumerate(self.time_range):
            # Layer tracking
            if self.precip[tstep] > 0.0:

                # Amount on the ground for overburden
                self.deposition[self.surface_layer_count] = self.precip[tstep]/1000.0

                # Assign air temp for bonding strength
                self.basal_temp[self.surface_layer_count] = self.air_temp[tstep]+273

                # Assign initial density
                self.rho[tstep,self.surface_layer_count] = self.new_density[tstep]

                self.surface_layer_count += 1

            # Iterate through the snow domain by layer
            #print("Timestep = {0}".format(tstep))
            depth = 0
            for layer in range(0,self.surface_layer_count):
                rho = self.rho[tstep,layer]
                T_z = self.basal_temp[layer]

                # Overburden
                ob = self.g * (np.cos(self.angle)**2) * np.sum(self.deposition[layer:])

                #while
                # Viscosity
                term1 = self.b1*np.exp(self.b2*(rho/self.rho_ice))
                term2 = np.exp(self.E/(self.R*T_z))
                visc = term1*term2

                # C = 1.0
                # visc = C * np.exp(self.E/(self.R * T_z)) * rho ** 4.0

                # Explicit density equation
                rho_1 = rho + self.dt * (rho * (1.0 / visc) * (75.0 + ob))

                #print("\tLayer = {0}".format(layer))
                depth += rho*self.deposition[layer]
                try:
                    self.rho[tstep+1,layer] = rho_1
                except Exception as e:
                    print(e)
                    print('Weird hour thing')
            print(depth)
if __name__ == '__main__':
    dem_ds = Dataset('./forcings/topo.nc')
    at_ds = Dataset('./forcings/air_temp.nc')
    pp_ds = Dataset('./forcings/precip.nc')
    rho_ds = Dataset('./forcings/snow_density.nc')

    #Calculate Slop angle
    sx,sy = np.gradient(dem_ds.variables['dem'])
    slope_im = np.arctan2(sy,sx)

    x = 5
    y = 5

    angle = np.abs(slope_im[y,x])
    angle_deg = angle*180/3.14

    if angle_deg >90:
        angle_deg = 180-angle_deg

    print("Slope Angle  = {0}".format(angle_deg))

    time = len(at_ds.variables['air_temp'][:])
    air_temp = (at_ds.variables['air_temp'])[0:time,y,x]
    precip = (pp_ds.variables['precip'])[0:time,y,x]
    rho_new = (rho_ds.variables['snow_density'])[0:time,y,x]

    sp = SNOSS(angle, air_temp, precip,rho_new)
    sp.solve()
    print(sp.rho)
    for i in range(len(sp.rho[0,:])):
        plt.plot(sp.rho[:,i])
    plt.show()
    # s.setup_domain()
    # s.setup_outputs('./forcings/stability.nc')
    # s.run_model()
    # s.close()
