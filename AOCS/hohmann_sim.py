from astropy import units as u
from poliastro.bodies import Body
from astropy.constants import Constant
import math
import numpy as np
from poliastro.twobody import Orbit
import matplotlib.pyplot as plt
from poliastro.plotting import *
import plotly.io as pio
from astropy.time import Time, TimeDelta
from poliastro.twobody.propagation import propagate, cowell
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from poliastro.core.elements import rv2coe
import serial

density = 1700 # kg/m
density_est = 400 # +/- estimation
dens_max = density + density_est
dens_min = density - density_est
diameter_didymos = 780 # m 

didymos_mass = (4/3)*math.pi*((diameter_didymos/2)**3)*density

didymos_volume = (4/3)*math.pi*((diameter_didymos/2)**3)

print(didymos_mass)


a = ((4/3)*math.pi*((diameter_didymos/2)**3))**(1./3.0)
G = 6.67408*10**(-20) #km3/kg s2
m = didymos_mass
k = G*m*u.km**3/u.s**2
name = "Dydimos"
#Dydimos 65803
didy = Body(None,k, name)

## Determining the position vector
# 4 km is the distance to the main body
phi = 0.0

v = (G*m/4.0)**(1.0/2.0)
vx = v*math.cos(phi+(math.pi/2.0))
vx = -v*math.sin(phi)
vy = v*math.cos(phi+(math.pi/2.0))
vy = v*math.cos(phi)

r = [4.0,0.0,0.0] * u.km
v = [0.0,7.270456978812678e-05,4.19761216313734e-05]* u.km /u.s
orbit = Orbit.from_vectors(didy, r, v)


plt.ylabel("h(t)")
plt.xlabel("t, days")

for i in range(1,10000):
    prestate = i*100-100
    state = i*100
    tofs = TimeDelta(np.linspace(prestate * u.s,state * u.s, num=2))
    with serial.Serial('/dev/ttyACM0', 115200, timeout=0.01) as ser:
        line = ser.readline()
        # Starting maneuver if signal from microcontroller is received
        if line == b'orbit transfer\r\n':
            print(line)
            from poliastro.maneuver import Maneuver
            dv = [-0.001, 0, 0] * u.m / u.s
            man = Maneuver.impulse(dv)
            orbit = Orbit.from_vectors(Earth, [rx,ry,rz]*u.km,[vx,vy,vz]*u.km/u.s)
            orbit = orbit.apply_maneuver(man)
        else:
            None
    rr = propagate(orbit,tofs,method=cowell)
    plt.plot(tofs.value, rr.norm() - (0.780/2 * u.km), color='blue')
    altitude = rr.norm().value - 0.780/2
    line = altitude[1]
    with serial.Serial('/dev/ttyACM0', 115200, timeout=0.001) as ser:
        ser.write(line)
        ser.write(state)
        # sending to the microncontroller the time to perform the second transfer
    plt.pause(0.00005)

plt.show()  
