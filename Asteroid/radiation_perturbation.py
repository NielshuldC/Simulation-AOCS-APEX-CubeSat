from astropy import units as u
from poliastro.bodies import Body
from astropy.constants import Constant
from poliastro.bodies import Sun
import math
import numpy as np
from poliastro.twobody import Orbit
import matplotlib.pyplot as plt
from poliastro.plotting import *
import plotly.io as pio
from astropy.time import Time, TimeDelta
from poliastro.twobody.propagation import propagate, cowell
from poliastro.core.perturbations import atmospheric_drag, third_body, J2_perturbation, shadow_function #cube_perturbation
from poliastro.bodies import Earth, Moon
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from poliastro.core.elements import rv2coe
from numpy.linalg import norm


density = 1700 # kg/m
density_est = 400 # +/- estimation
dens_max = density + density_est
dens_min = density - density_est
diameter_didymos = 780 # m 

didymos_mass = (4/3)*math.pi*((diameter_didymos/2)**3)*density

didymos_volume = (4/3)*math.pi*((diameter_didymos/2)**3)

print(didymos_mass)

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
print(vx)
print(vy)


r = [4.0,0.0,0.0] * u.km
v = [0.0,7.270456978812678e-05,4.19761216313734e-05]* u.km /u.s
ss = Orbit.from_vectors(didy, r, v)
ss.plot()
plt.show()

def shadow_functions(r_sat, r_sun, R):
    
    r_sat_norm = np.sqrt(np.sum(r_sat ** 2))
    r_sun_norm = np.sqrt(np.sum(r_sun ** 2))

    theta = np.arccos(np.dot(r_sat, r_sun) / r_sat_norm / r_sun_norm)
    theta_1 = np.arccos(R / r_sat_norm)
    theta_2 = np.arccos(R / r_sun_norm)

    return theta < theta_1 + theta_2

def radiation_pressures(t0, state, k, R, C_R, A, m, Wdivc_s):

    r_star = 69769641#,-2.02846764e+08,-96346002.511]
    r_sat = state[:3]
    P_s = Wdivc_s / (69769641 ** 2)

    nu = float(np.linspace(shadow_functions(r_sat, r_star, R), endpoint=True))
    return [P_s * (C_R * A / m),P_s * (C_R * A / m),P_s * (C_R * A / m)]

rr = propagate(ss,tofs,method=cowell,ad=radiation_pressures,R=0.78,C_R=2,A=A_APEX,m=11.86,Wdivc_s=Wdivc)

fig2 = plt.figure(figsize=(8, 7))
ax = fig2.add_subplot(111, projection='3d')
rc = [-0.68, 0.68]
for s, e in combinations(np.array(list(product(rc, rc, rc))), 2):
    if np.sum(np.abs(s-e)) == rc[1]-rc[0]:
        ax.plot3D(*zip(s, e), color="grey")
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

