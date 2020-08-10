
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris

from poliastro.ephem import build_ephem_interpolant
from poliastro.core.elements import rv2coe

from poliastro.core.util import norm
from poliastro.core.perturbations import atmospheric_drag, third_body, J2_perturbation
from poliastro.bodies import Earth, Moon
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter3D
from plotly.graph_objs import Figure, FigureWidget
from poliastro.twobody.propagation import propagate, cowell
from poliastro.core.perturbations import atmospheric_drag, third_body, J2_perturbation #cube_perturbation
from poliastro.bodies import Earth, Moon
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from poliastro.core.elements import rv2coe
import math
from poliastro.bodies import Body


import plotly.io as pio
pio.renderers.default = "notebook_connected"

## Body properties

density = 1700 # kg/m
density_est = 400 # +/- estimation
dens_max = density + density_est
dens_min = density - density_est
diameter_didymos = 780 # m 
didymos_mass = (4/3)*math.pi*((diameter_didymos/2)**3)*density
didymos_volume = (4/3)*math.pi*((diameter_didymos/2)**3)
a = ((4/3)*math.pi*((diameter_didymos/2)**3))**(1./3.0)
G = 6.67408*10**(-20) #km3/kg s2
m = didymos_mass
k = G*m*u.km**3/u.s**2
name = "Dydimos"
didy = Body(None,k, name)

r = [4.0,0.0,0.0] * u.km
v = [0.0,7.270456978812678e-05,4.19761216313734e-05]* u.km /u.s
ss = Orbit.from_vectors(didy, r, v)

## Maneuver values

from poliastro.maneuver import Maneuver
dv = [25, 0, 0] * u.m / u.s
man = Maneuver.impulse(dv)
man = Maneuver((0 * u.s, dv))

# Plot and propagate

plt.ylabel("h(t)")
plt.xlabel("t, days")

import serial

def cube_perturbation(xp,yp,zp):

    density = 1700 # g/cm3

    diameter_didymos = 780 # m 

    didymos_mass = (4/3)*math.pi*(diameter_didymos^3)*density

    a = (4/3)*math.pi*(diameter_didymos^3)**(1/3)
    G = 6.67408*10**(-11)
    m = didymos_mass
    
    a_x = 8.319137045e-12
    a_y = 8.319137045e-12
    a_z = -8.319137045e-12
    #return np.array([a_x[0], a_y[1], a_z[2]])
    return np.array([a_x,a_y,a_z])

tofs = TimeDelta(np.linspace(0,10000000 * u.s, num=1000))
rr = propagate(ss,tofs,method=cowell,ad=cube_perturbation)
fig1 = plt.figure(figsize=(8, 7))
plt.ylabel("h(t)")
plt.xlabel("t,days")
plt.plot(tofs.value, rr.norm()-(0.780/2)*u.km)
plt.show()

for i in range(1,1000):
    prestate = i*100-100
    state = i*100
    tofs = TimeDelta(np.linspace(prestate * u.s,state * u.s, num=2))
    print(prestate)
    print(state)
    #print(tofs)
    ser = serial.Serial('/dev/ttyACM0')
    #print(ser.name)
    with serial.Serial('/dev/ttyACM0', 115200, timeout=0.1) as ser:
        line = ser.readline()
        if line == b'button pressed\r\n':
            print(line)
            ss = ss.apply_maneuver(man)
        else:
            None
    rr = propagate(ss,tofs,method=cowell,ad=cube_perturbation)
    plt.plot(tofs.value, rr.norm()-(0.780/2)*u.km, color='blue')
    plt.pause(0.00005)

plt.show()
