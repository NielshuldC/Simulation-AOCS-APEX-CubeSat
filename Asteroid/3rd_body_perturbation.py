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
from poliastro.plotting import OrbitPlotter2D
from plotly.graph_objs import Figure, FigureWidget
from poliastro.ephem import build_ephem_interpolant
from scipy import interpolate 

epoch = Time(
    2454283.0, format="jd", scale="tdb"
)  # setting the exact event date is important

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

#U = G*m/r - ((7*a^(4)*G*m)/(30*r^(9)))*(x^(4)+y^(4)-3*(x^(2)*y^(2)+x^(2)*z^(2)+y^(2)*z^(2)))

#x = 4000
#y = 0.0
#z = 0.0

#Fx = (7.0*a*4.0*G*m*(x**4.0-5.0*x**(2.0)*(y**2.0+z**2.0)+(y**4.0-y**(2.0)*z**(2.0)+z**(4.0))))/(6.0*(x**(2.0)+y**(2.0)+z**(2.0))**(11.0/2.0))
#print("FX",Fx)
#Fy = (7.0*a*4.0*G*m*y*(3.0*x**4.0+y**4.0-5.0*y**(2.0)*z**2.0+3.0*z**4.0-x**(2.0)*(5.0*y**2.0+3.0*z**2.0)))/(6.0*(x**2.0+y**2.0+z**2.0)**(11.0/2.0))
#Fz = (7.0*a*4.0*G*m*z*(3.0*x**4.0 + 3.0*y**4.0 -5.0*y**(2.0)*z**2.0+z**4.0 -x**(2.0)*(3.0*y**2.0+5.0*z**2.0)))/(6.0*(x**2.0+y**2.0+z**2.0)**(11.0/2.0))
#print("Fy",Fy)
#print("Fz",Fz)
#Ft = math.sqrt(Fx**2.0+Fy**2.0+Fz**2.0)

k = G*m*u.km**3/u.s**2
print(k)
name = "Dydimos"
#Dydimos 65803
didy = Body(None,k, name)

## Determining the position vector
# 4 km is the distance to the main body
phi = 0.0

#rx = 4*math.cos(phi)
#ry = 4*math.sin(phi)
#print(rx)
#print(ry)
v = (G*m/4.0)**(1.0/2.0)
vx = v*math.cos(phi+(math.pi/2.0))
vx = -v*math.sin(phi)
vy = v*math.cos(phi+(math.pi/2.0))
vy = v*math.cos(phi)
print(vx)
print(vy)

#vvv = 8.395206522348602e-05
#psi = 0.5236
#vvx = v*math.sin(psi)
#vvy = v*math.cos(psi)
#print(vvx)
#print(vvy)
#4.19761216313734e-05
#7.270456978812678e-05

r = [4.0,0.0,0.0] * u.km
#v = [0.0,0.00030501089648515324,0.00010501089648515324]* u.km /u.s
v = [0.0,7.270456978812678e-05,4.19761216313734e-05]* u.km /u.s
#print(v)
#v = [-3.457, 6.618, 2.533] * u.km /u.s
ss = Orbit.from_vectors(didy, r, v)
#plot(ss, label= 'sample orbit')
ss.plot()
plt.show()

#xp = 0.0
#yp = 0.0
#zp = 0.0

#def cube_perturbation(xp,yp,zp):

#    density = 1700 # g/cm3

#    diameter_didymos = 780 # m 

#    didymos_mass = (4/3)*math.pi*(diameter_didymos^3)*density

#    a = (4/3)*math.pi*(diameter_didymos^3)**(1/3)
#    G = 6.67408*10**(-11)
#    m = didymos_mass
#    a_x = (7.0*a**(4.0)*G*m*xp*(xp**(4.0)-5.0*xp**(2.0)*(yp**2.0+zp**2.0)+3*(yp**(4.0)-yp**(2.0)*zp**(2.0)+zp**(4.0))))/(6.0*(xp**(2.0)+yp**(2.0)+zp**(2.0))**(11.0/2.0))
#    print("AX",a_x)
#    a_y = (7.0*a**(4.0)*G*m*yp*(3.0*xp**(4.0)+yp**(4.0)-5.0*yp**(2.0)*zp**(2.0)+3.0*zp**(4.0)-xp**(2.0)*(5.0*yp**(2.0)+3.0*zp**2.0)))/   (6.0*(xp**2.0+yp**2.0+zp**2.0)**(11.0/2.0))
#    print(a_y)
#    a_z = (7.0*a**(4.0)*G*m*zp*(3.0*xp**4.0 + 3.0*yp**4.0 -5.0*yp**(2.0)*zp**2.0+zp**4.0 -xp**(2.0)*(3.0*yp**2.0+5.0*zp**2.0)))/(6.0*(xp**2.0+yp**2.0+zp**2.0)**(11.0/2.0))
#    print(a_y)
#    return np.array([a_x[0], a_y[1], a_z[2]])


a_1=780
b_1=480
f=(a_1-b_1)/a_1 # body oblatness
w_b= 0 # central body rotation rate
J = (2*f/3) # - (0.780**3**w_b**2)/(3*G*m) nodal precession


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
#moon_r = [1.18,0.0,0.0,0.0,0.00015456828072977154,0.0]
moon_v = [0.0,0.00015456828072977154,0.0]* u.km /u.s

rm = Orbit.from_vectors(didy, moon_r, moon_v)
fig = FigureWidget()
op = OrbitPlotter2D(fig)
op.plot(rm, label="didymoon")
op.plot(ss, label="APEX")

fig.show()
rm_rr = rm.propagate(11440* u.min)
print(rm_rr)

#body_r = build_ephem_interpolant(
#   Moon, 28 * u.day, (epoch.value * u.day, epoch.value * u.day + 60 * u.day)

t = 5184000

def perturbation_body(t):
    r, v = rm.propagate(t * u.s).rv()
    return interpolate.interp1d(r.to_value(u.km),v.to_value(u.km / u.s))

moon_per = perturbation_body(t)

tofs = TimeDelta(np.linspace(0,60 * u.day, num=1000))
#tofs = TimeDelta(np.linspace(prestate * u.s,state * u.s, num=2))
#rr = propagate(ss,tofs,method=cowell,ad=J2_perturbation,J2=J,R=0.780*u.km)
rr = propagate(
    ss,
    tofs,
    method=cowell,
    ad=third_body,
    k_third=k_moon,
    third_body=moon_per,
)
#rr2 = propagate(s0, tofs, method=cowell, rtol=1e-6)

#for i in range(1,100):
    #prestate = i*100-100
    #state = i*100
    #tofs = TimeDelta(np.linspace(prestate * u.s,state * u.s, num=2))
    #print(prestate)
    #print(state)
    #print(tofs)
    #x, z, y, G, m, a)
    #print(xp)
    #rr = propagate(ss,tofs,method=cowell)
    #rr = propagate(ss,tofs,method=cowell,ad=cube_perturbation)
    #rr = propagate(ss,tofs,method=cowell,ad=cube_perturbation,xp=xpp,yp=ypp,zp=zpp,G=G,m=m,a=a)
    #rr = propagate(ss,tofs,method=cowell,ad=atmospheric_drag,R=R,C_D=C_D,A=A,m=m,H0=H0,rho0=rho0)
    #rr = propagate(orbit,tofs,method=cowell)
    #plt.plot(tofs.value, rr.norm() - Earth.R, color='blue')
    #x = rr.x
    #y = rr.y
    #z = rr.z
    #print(rr.x)
    

#plt.pause(0.00005)
fig2 = plt.figure(figsize=(8, 7))
ax = fig2.add_subplot(111, projection='3d')
#rc = [-0.68, 0.68]
#for s, e in combinations(np.array(list(product(rc, rc, rc))), 2):
#    if np.sum(np.abs(s-e)) == rc[1]-rc[0]:
#        ax.plot3D(*zip(s, e), color="r")
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

