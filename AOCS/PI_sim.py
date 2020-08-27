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

#DIDIMOON


diameter_didymoon = 160 # m 

didymoon_mass = (4/3)*math.pi*((diameter_didymoon/2)**3)*density

didymoon_volume = (4/3)*math.pi*((diameter_didymoon/2)**3)

print(didymos_mass)

#a = ((4/3)*math.pi*((diameter_didymos/2)**3))**(1./3.0)
G = 6.67408*10**(-20) #km3/kg s2
mm = didymoon_mass
k_moon = G*mm*u.km**3/u.s**2
name_moon = "didymoon"
didymoon = Body(None,k_moon, name_moon)

moon_r = [1.18,0.0,0.0]*u.km
moon_v = [0.0,0.00015456828072977154,0.0]* u.km /u.s

rm = Orbit.from_vectors(didy, moon_r, moon_v)
fig = FigureWidget()
op = OrbitPlotter2D(fig)
op.plot(rm, label="didymoon")
op.plot(ss, label="APEX")

fig.show()
rm_rr = rm.propagate(11440* u.min)
print(rm_rr)


t = 5184000

def perturbation_body(t):
    r, v = rm.propagate(t * u.s).rv()
    return interpolate.interp1d(r.to_value(u.km),v.to_value(u.km / u.s))

moon_per = perturbation_body(t)


plt.ylabel("h(t)")
plt.xlabel("t, days")

for i in range(1,10000):
    prestate = i*100-100
    state = i*100
    tofs = TimeDelta(np.linspace(prestate * u.s,state * u.s, num=2))
    with serial.Serial('/dev/ttyACM0', 115200, timeout=0.01) as ser:
        line = ser.readline()
        # Starting maneuver if signal from microcontroller is received
        line == msg.decode('utf-8')
        from poliastro.maneuver import Maneuver
        dv = [line[0], line[1], 0] * u.m / u.s
        man = Maneuver.impulse(dv)
        orbit = Orbit.from_vectors(Earth, [rx,ry,rz]*u.km,[vx,vy,vz]*u.km/u.s)
        orbit = orbit.apply_maneuver(man)
    rr = propagate(
    orbit,
    tofs,
    method=cowell,
    ad=third_body,
    k_third=k_moon,
    third_body=moon_per,
    )
    plt.plot(tofs.value, rr.norm() - (0.780/2 * u.km), color='blue')
    altitude = rr.norm().value - 0.780/2
    line = altitude[1]
    line1 = rr.x.value[0]
    line2 = rr.y.value[0]
    line3 = (rr.norm().value[1] - 0.780/2)*(state-prestate)
    with serial.Serial('/dev/ttyACM0', 115200, timeout=0.001) as ser:
        ser.write(line)
        ser.write(line1)
        ser.write(line2)
        ser.write(line3)
    plt.pause(0.00005)

plt.show()


fig2 = plt.figure(figsize=(8, 7))
ax = fig2.add_subplot(111, projection='3d')

uu, vw = np.mgrid[0:4*np.pi:80j, 0:np.pi:80j]
xx = np.cos(uu)*np.sin(vw)
yy = np.sin(uu)*np.sin(vw)
zz = np.cos(vw)
ax.plot_wireframe(xx, yy, zz, color="grey")
ax.plot(rr.x, rr.y, rr.z)
plt.show()

fig3 = plt.figure(figsize=(8, 7))
plt.ylabel("h(t)")
plt.xlabel("t,days")
plt.plot(tofs.value, rr.norm()-(0.780/2)*u.km)
plt.show()

coords = rr.xyz.T.to(u.km).value
vv = rr.differentials["s"].d_xyz.T.to(u.km / u.s).value
raans = [rv2coe(k,r,v)[3] for r,v in zip(coords,vv)]
fig3 = plt.figure(figsize=(8, 7))
plt.ylabel("RAAN")
plt.xlabel("t,days")
plt.plot(tofs.value, raans)
plt.show()
