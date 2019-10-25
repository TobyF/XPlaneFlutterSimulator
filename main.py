import time
from time import sleep
import xpc
from scipy import signal
from scipy.interpolate import interp1d
from scipy import integrate
import h5py
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import collections
from geopy import distance
import socket

#TODO sort out bank and heading index values and check vertical motion is as expected - Done

def ex():
    print("Setting up simulation")

    with xpc.XPlaneConnect(timeout=1000) as client:
        # Verify connection
        try:
            # If X-Plane does not respond to the request, a timeout error
            # will be raised.
            client.getDREF("sim/test/test_float")
        except:
            print("Error establishing connection to X-Plane.")
            print ("Exiting...")
            return

        print("Entering Full Simulation")

        pe = PhysicsEngine(client)
        count = 0
        start_t = time.time()
        state = 'xp' #or 'sharpy'
        while True:  # main loop
            count += 1
            #print("Count:"+str(count))
            #time.sleep(0.001)
            try:
                if pe.controls[5] == 0:
                    if state =='xp':
                        pe.client.pauseSim(True)
                        pe.get_position()
                        pe.reset_position()
                        pe.update_position()
                        state = 'sharpy'
                    pe.get_position()
                    pe.get_controls()
                    pe.calculate_position()
                    pe.update_position()
                else:
                    if state =='sharpy':
                        state = 'xp'
                        print("Not Simulating, raise flaps, press 1")
                    pe.client.pauseSim(False)
                    pe.get_controls()
                    pe.get_position()
                    time.sleep(0.5)

            except socket.timeout:
                print("Missed an update")
                #sleep(10)
                #break  # if user pressed a key other than the given key the loop will break
            if count > 10000:
                break
        finish_t = time.time()
        print("Elapsed time")
        print(finish_t - start_t)
        print("Frequency")
        print(1000.0/(finish_t-start_t))

        input("Press any key to exit...")

def newLatLon(lat1,lon1,bearing,distance_m):
    E_radius = 6371000
    lat1R = np.deg2rad(lat1)
    lon1R = np.deg2rad(lon1)
    brgR  = np.deg2rad(bearing)
    d_ratio = distance_m/E_radius

    lat2R = np.arcsin(np.sin(lat1R)*np.cos(d_ratio)+np.cos(lat1R)*np.sin(d_ratio)*np.cos(brgR))

class PhysicsEngine:
    def __init__(self, client,mats = 'ss_horten.h5'):
        self.client = client
        self.position = []
        self.velocity = []
        self.controls = []
        self.last_update_time = time.time()
        self.times = []
        self.count = 0
        self.control_buffer = collections.deque([],maxlen = 5) #Only store 5 most recent inputs
        self.get_controls()
        time.sleep(1)
        self.get_controls()

        #Changes between angle and x-planes -1 to 1 scale
        #-1 is stick forward, positive elevator (down), nose down result.
        self.elev_angle_f = interp1d([1,-1],[np.deg2rad(-30),np.deg2rad(30)])

        #Initialises the matricies
        h_mats = h5py.File(mats,'r')
        #print(list(h_mats.keys()))
        self.a  = np.loadtxt('A.dat')#np.array(h_mats['A'])
        self.b  = np.swapaxes([np.loadtxt('B.dat')],0,1)#np.swapaxes(np.array([h_mats['B']]),0,1)
        self.c  = np.loadtxt('C.dat')#np.array(h_mats['C'])
        self.d  = np.swapaxes([np.loadtxt('D.dat')],0,1)#np.swapaxes(np.array([h_mats['D']]),0,1)
        self.dt = 0.008928 #Need a way of importing this from hdf5 file
        print(self.b.shape)
        self.x = np.zeros(len(self.b[1])) #initialise the equilibrium state


    def get_position(self):
        self.position = list(self.client.getPOSI())

        #print(self.position)


    def get_controls(self):
        #Gets control positions from x-Plane and saves them in a buffer (buffer not implemented)
        # 0 - Elevator
        # 1 - Ailerons...
        self.controls = self.client.getCTRL()
        #print(self.controls)
        self.control_buffer.append((time.time(),self.controls))
        print("Elevator State:"+str(self.controls[0]))

    def reset_position(self):
        self.position[3] = 0 # Sets pitch to 0
        self.position[4] = 0 # Sets bank to 0
        #self.position[5] = 0 # Sets hdg to North
        self.velocity = [28, 0, 0] #Sets speed to 28 m/s fwd, level flight
        print("Pitch Reset to 0, Speed set to 28m/s")

    def calculate_position(self):
        #Create a time array,t and input array u
        t = np.arange(self.control_buffer[-2][0],self.control_buffer[-1][0],self.dt)
        u = self.elev_angle_f(np.linspace(self.control_buffer[-2][1][0],self.control_buffer[-1][1][0],len(t)))
        #print("Simulating with "+str(len(u))+" timesteps")
        #print(f"Inputs: {u}")
        #print(f"Initial State: {self.x}")
        #Run the state space sim
        tout ,yout, xout = signal.dlsim((self.a,self.b,self.c,self.d,self.dt),u,x0 =self.x)
        self.x = xout[-1] #Stores the current state for the next time.
        #print(f"Yout: {yout}")
        #integrate the yout
        pitch = integrate.cumtrapz(np.array(yout)[:,4],dx = self.dt,initial=0) + np.deg2rad(self.position[3]) #Integrate the pitch rate for pitch change

        U_e = np.multiply(np.cos(pitch),28-np.array(yout)[:,0])-np.multiply(np.sin(pitch),np.array(yout)[:,2])
        V_e = np.multiply(np.sin(pitch),28-np.array(yout)[:,0])+np.multiply(np.cos(pitch),np.array(yout)[:,2])

        #U_e_vel = integrate.cumtrapz(dU_e,dx = self.dt,initial=0) + self.velocity[0]
        #V_e_vel = integrate.cumtrapz(dV_e,dx = self.dt,initial=0) + self.velocity[2]

        self.velocity[0] = U_e[-1]
        self.velocity[2] = V_e[-1]

        dX = integrate.simps(U_e,dx = self.dt)
        dZ = integrate.simps(V_e,dx = self.dt)
        print(f'Moving at speed {self.velocity[0]},Moved fwd {dX} m and up {dZ} m')

        #Update the pitch
        self.position[3] = np.rad2deg(pitch[-1])
        new_cords = distance.distance(meters=dX).destination((self.position[0],self.position[1]),self.position[5]) # Calculate new lat/lon
        #print(f'new_cords {new_cords}')
        self.position[0] = new_cords[0]
        self.position[1] = new_cords[1]
        self.position[2] += dZ
        #self.velocity = [U_vel[-1]]
        #print(f"New Position: {self.position}")
        #input()

        #Display timing info
        delay = time.time() - self.control_buffer[-1][0]
        #print("It has been "+str(delay)+ "seconds!")

        #Return new position tuple
        return self.position

    def update_position(self):
        t1 = timer()
        self.client.sendPOSI(self.position)
        t2 = timer()
        #print("Send Posi, took:"+str(t2-t1))
        self.times.append(t2-t1)
        self.count += 1
        #print("avg:" +str(sum(self.times)/self.count))


        print("New Pitch:"+str(self.position[3] ))

if __name__ == "__main__":
    ex()
