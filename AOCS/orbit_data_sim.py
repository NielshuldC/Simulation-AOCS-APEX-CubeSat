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
from poliastro.core.perturbations import atmospheric_drag, third_body, J2_perturbation #cube_perturbation
from poliastro.bodies import Earth, Moon
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


print("Cube Edge :",a)
print("Gravitational constant :",G)
print("Mass of the body :",m)


k = G*m*u.km**3/u.s**2
print(k)
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
#v = [0.0,0.00030501089648515324,0.00010501089648515324]* u.km /u.s
v = [0.0,7.270456978812678e-05,4.19761216313734e-05]* u.km /u.s
#print(v)
#v = [-3.457, 6.618, 2.533] * u.km /u.s
ss = Orbit.from_vectors(didy, r, v)

xp = 4.0
yp = 0.0
zp = 0.0

def cube_perturbation(xp,yp,zp):

    density = 1700 # g/cm3

    diameter_didymos = 780 # m 

    didymos_mass = (4/3)*math.pi*(diameter_didymos^3)*density

    a = (4/3)*math.pi*(diameter_didymos^3)**(1/3)
    G = 6.67408*10**(-11)
    m = didymos_mass
    a_x = (7.0*a**(4.0)*G*m*xp*(xp**(4.0)-5.0*xp**(2.0)*(yp**2.0+zp**2.0)+3*(yp**(4.0)-yp**(2.0)*zp**(2.0)+zp**(4.0))))/(6.0*(xp**(2.0)+yp**(2.0)+zp**(2.0))**(11.0/2.0))
    a_y = (7.0*a**(4.0)*G*m*yp*(3.0*xp**(4.0)+yp**(4.0)-5.0*yp**(2.0)*zp**(2.0)+3.0*zp**(4.0)-xp**(2.0)*(5.0*yp**(2.0)+3.0*zp**2.0)))/   (6.0*(xp**2.0+yp**2.0+zp**2.0)**(11.0/2.0))
    a_z = (7.0*a**(4.0)*G*m*zp*(3.0*xp**4.0 + 3.0*yp**4.0 -5.0*yp**(2.0)*zp**2.0+zp**4.0 -xp**(2.0)*(3.0*yp**2.0+5.0*zp**2.0)))/(6.0*(xp**2.0+yp**2.0+zp**2.0)**(11.0/2.0))
   
    return np.array([a_x,a_y,a_z])

plt.ylabel("h(t)")
plt.xlabel("t, days")

for i in range(1,10000):
    prestate = i*100-100
    state = i*100
    tofs = TimeDelta(np.linspace(prestate * u.s,state * u.s, num=2))
    print(prestate)
    print(state)
    print(tofs)
    rr = propagate(ss,tofs,method=cowell,ad=cube_perturbation)
    plt.plot(tofs.value, rr.norm() - (0.780/2 * u.km), color='blue')
    xp=rr.x[1]
    yp=rr.y[1]
    zp=rr.z[1]
    liner = rr.norm().value - 0.780/2
    line = liner[1]
    print(line)
    with serial.Serial('/dev/ttyACM0', 115200, timeout=0.001) as ser:
        ser.write(line)
        
    plt.pause(0.00005)

plt.show()  





