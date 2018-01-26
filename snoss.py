#Work flow.


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

class SpatialSNOSS(object):
    def __init__(self,dem,air_temp,precip, snow_density):
        self.topo = Dataset(dem)
        self.temperature_ds = Dataset(air_temp)
        self.precip_ds = Dataset(precip)
        self.snow_rho_ds = Dataset(snow_density)

        self.rho_ice = 917.0 # density of ice
        self.g = 9.81 #Gravity
        self.b1 = 6.5*10**-7 #viscosity constant
        self.b2 = 19.3 #viscosity constant
        self.E = 67.3 #Activation energy
        self.R = 0.0083 #Gas constant


    def setup_outputs(self,out_f):
        """
        sets up the output netcdf file
        """
        self.out = Dataset(out_f, "w", format="NETCDF4",clobber=True)

        self.out.createDimension('time',self.nt)
        self.out.createDimension('x',self.nx)
        self.out.createDimension('y',self.ny)
        self.out.createDimension('z',self.nz)

        self.out.createVariable('time','f',dimensions=('time'))
        self.out.createVariable('x','f',dimensions=('x'))
        self.out.createVariable('y','f',dimensions=('y'))
        self.out.createVariable('z','f',dimensions=('z'))

        self.out.createVariable('depth','f',dimensions=('time','y','x'))
        self.out.createVariable('min_stability','f',dimensions=('time','y','x'))

        self.out.createVariable('stability','f',dimensions=('time','y','x','z'))
        self.out.createVariable('snow_density','f',dimensions=('time','y','x','z'))
        self.out.createVariable('basal_temp','f',dimensions=('time','y','x','z'))

        #Assign domain value to output
        self.out.variables['x'] = self.topo.variables['x'][:]
        self.out.variables['y'] = self.topo.variables['y'][:]
        self.out.variables['z'] = np.array(np.arange(0,self.max_depth+self.resolution,self.resolution))


    def setup_domain(self, max_depth=5.0, resolution = 0.1):
        """
        Takes in the dem filename and the max depth to be modeled
        resolution is  number of layers to track over the potential depth domain
        """
        self.resolution = resolution # Z cell size
        self.max_depth = max_depth

        #discretize
        self.nx = self.topo.dimensions['x'].size
        self.ny = self.topo.dimensions['y'].size
        self.nz = int(self.max_depth/self.resolution)

        self.nt = self.precip_ds.variables['time'].size

        #Assumes we start at 0 hours
        self.dt = (self.precip_ds.variables['time'][-1])/self.precip_ds.variables['time'].size


    def print_intro(self):
        msg = "SPATIAL SNOSS"

        print("="*len(msg))
        print(msg)
        print("="*len(msg))
        print("\nModels the stability index and outputs a single netcdf"
              "\ncontaining the modeled:"
              "\n\n\t * depth (2D)"
              "\n\t * stability index (3D)"
              "\n\t * minimum stability index (2D)")


    def print_overview(self):
        msg = "MODEL OVERVIEW:"
        print("\n{0}".format(msg))
        print("-"*len(msg))

        print("Model Time: {0} Hours".format(self.nt))
        print("Domain in cells: {0} X {1} X {2}".format(self.nx,self.ny,self.nz))
        print("Total Cells: {0}".format(self.nx*self.ny,self.nz))

    def calculate_viscosity(self,t,x,y,z):
        rho_z = self.rho[t,y,x,z]
        T_z = self.basal_temp
        term1 = self.b1*np.exp(self.b2*(rho_z/self.rho_ice))
        term2 = np.exp(self.E/(self.R*T_z))
        T_s = self.basal_temp[y,x,z]

        return term1*term2

    def calculate_overburden(self,t,x,y,z):
        if z <= self.depth[t,y,x]:

            over_head_mass = np.sum(self.precip[z:t,y,x])
            result = self.g*self.cos(self.slope[y,x])
        else:
            result = 0.0

        return result

    def initialize_model(self):
        #Print out info
        self.print_intro()
        self.print_overview()

        #Make variables easy to use
        self.air_temp = self.temperature_ds.variables['air_temp']
        self.precip = self.precip_ds.variables['precip']
        self.rho = self.snow_rho_ds.variables['snow_density']

        self.depth = self.out.variables['depth']
        self.stability = self.out.variables['stability']
        self.basal_temp = self.out.variables['basal_temp']
        self.z = self.out.variables['z']

        #Calculate first time step:
        #Calculate the new snow depth
        snow_id = self.rho[0] > 0
        data = np.zeros((self.ny,self.nx))
        data[snow_id] =  0.001*(self.precip[t][snow_id])/self.rho[t][snow_id]
        self.depth[t] = data

        self.depth_index =

        for y in range(self.ny):
            for x in range(self.nx):
                for z in range(self.nz):
                    n_zz = self.calculate_viscosity(t,y,x,z)
                    sigma_zz = self.calculate_overburden(t,y,x,z)
                    sigma_m = 75.0
                    self.rho[t,y,x,z] =

    def run_model(self):


        #Iterate through all time,space
        for t in range(1,self.nt):

            #Calculate the new snow depth
            snow_id = self.rho[t] > 0
            data = np.zeros((self.ny,self.nx))
            data[snow_id] =  0.001*(self.precip[t][snow_id])/self.rho[t][snow_id]
            if t != 0:
                self.depth[t] = data+self.depth[t-1]
            else:
                self.depth[t] = data

            self.depth_index =

            for y in range(self.ny):
                for x in range(self.nx):
                    for z in range(self.nz):
                        n_zz = self.calculate_viscosity(t,y,x,z)
                        sigma_zz = self.calculate_overburden(t,y,x,z)
                        sigma_m = 75.0
                        self.rho[t,y,x,z] =
    def close(self):
        self.temperature_ds.close()
        self.precip_ds.close()
        self.snow_rho_ds.close()
        self.out.close()

if __name__ == '__main__':
    s = SpatialSNOSS('./forcings/topo.nc',
                     './forcings/air_temp.nc',
                     './forcings/precip.nc',
                     './forcings/snow_density.nc')

    s.setup_domain()
    s.setup_outputs('./forcings/stability.nc')
    s.run_model()
    s.close()
