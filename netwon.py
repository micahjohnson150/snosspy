import numpy as np


def snoss_netwon():
    """
    Solves the next timestep in the snoss equation using the
    netwon solver

    rho_1 = rho_0 + f(rho_0)/f(rho_0)'
    """




class SpatialSNOSS(object):
    def __init__(self,dem,air_temp,precip, snow_density):
        self.topo = Dataset(dem)
        self.temperature_ds = Dataset(air_temp)
        self.precip_ds = Dataset(precip)
        self.new_rho_ds = Dataset(snow_density)

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
        print("Setting up outputs file at :\n {0}".format(out_f))
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
        self.out.createVariable('basal_temp','f',dimensions=('y','x','z'))

        #Assign domain value to output
        self.out.variables['x'] = self.topo.variables['x'][:]
        self.out.variables['y'] = self.topo.variables['y'][:]
        self.out.variables['z'] = np.array(np.arange(0,self.max_depth+self.resolution,self.resolution))


    def calculate_slope(self):
        print('Calculating slope...')
        sx,sy = np.gradient(self.topo.variables['dem'])
        slope = np.arctan2(sy,sx)
        return slope


    def setup_domain(self, max_depth=5.0, resolution = 0.1):
        """
        resolution is z domain spacing to track over the potential depth domain
        """
        print("Setting up model domain...")
        self.resolution = resolution # Z cell size
        self.max_depth = max_depth

        #discretize
        self.nx = self.topo.dimensions['x'].size
        self.ny = self.topo.dimensions['y'].size
        self.nz = int(self.max_depth/self.resolution)

        self.nt = self.precip_ds.variables['time'].size

        #Assumes we start at 0 hours
        self.dt = (self.precip_ds.variables['time'][-1])/self.precip_ds.variables['time'].size
        self.calculate_slope()
        self.slope = self.topo.variables['dem']


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


    def calculate_viscosity(self,t,y,x,z):
        #print("Calculating viscosity...")
        rho_z = self.rho[t,y,x,z]
        T_z = self.basal_temp[y,x,z]
        term1 = self.b1*np.exp(self.b2*(rho_z/self.rho_ice))
        term2 = np.exp(self.E/(self.R*T_z))
        return term1*term2


    def calculate_overburden(self,t,y,x,z):
        #print("Calculating sigma_zz...")

        if self.z[z] <= self.depth[t,y,x]:
            over_head_mass = np.sum(self.precip[z:t,y,x])
            result = over_head_mass * self.g*(np.cos(self.slope[y,x]))**2

        else:
            result = 0.0


        return result


    def calculate_density(self,t,y,x,z):
        #print("Calculating density...")

        n_zz = self.calculate_viscosity(t-1,y,x,z)
        sigma_zz = self.calculate_overburden(t-1,y,x,z)
        sigma_m = 75.0
        #RHS
        rhs = self.rho[t-1,y,x,z]*((sigma_zz+sigma_m)/n_zz)

        return self.dt*60.0*60.0*rhs+self.rho[t-1,y,x,z]


    def initialize_model(self):
        #Print out info
        print("Initializing model run...")
        self.print_intro()
        self.print_overview()

        #Make variables easy to use
        self.air_temp = self.temperature_ds.variables['air_temp']
        self.precip = self.precip_ds.variables['precip']
        self.new_rho = self.new_rho_ds.variables['snow_density']

        self.depth = self.out.variables['depth']
        self.stability = self.out.variables['stability']
        self.basal_temp = self.out.variables['basal_temp']
        self.rho = self.out.variables['snow_density']

        self.z = self.out.variables['z']

        #Calculate first time step:
        #Calculate the new snow depth
        snow_id = self.new_rho[0] > 0
        data = np.zeros((self.ny,self.nx))
        data[snow_id] =  0.001*(self.precip[0][snow_id])/np.array(self.new_rho[0][snow_id])
        self.depth[0] = data

        self.rho[0,:,:,0] = self.new_rho[0]
        self.basal_temp[:,:,0] = self.air_temp[0]
        self.slope = self.calculate_slope()


    def run_model(self):
        print("Running simulation...")
        self.initialize_model()

        #Iterate through all time,space
        for t in range(1,self.nt):
            print("Solving for hour: %s" % t)

            #Calculate the new snow depth
            data = np.zeros((self.ny,self.nx))

            data = self.precip[t]/self.new_rho[t]

            self.depth[t] = data + self.depth[t-1]

            for z in range(self.nz):
                self.basal_temp[:,:,z] = self.air_temp[t]
                # for y in range(self.ny):
                #     for x in range(self.nx):
                #         if self.z[z] <= self.depth[t,y,x]:
                x = 7
                y = 7
                #print(x,y,z)
                self.rho[t,y,x,z] = self.calculate_density(t,y,x,z)

        plt.plot([np.average(self.depth[t]) for t in range(self.nt)])
        plt.show()

    def close(self):
        self.temperature_ds.close()
        self.precip_ds.close()
        self.new_rho_ds.close()
        self.out.close()
