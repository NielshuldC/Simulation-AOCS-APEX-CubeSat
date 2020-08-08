from astropy.utils.data import conf
conf.dataurl
conf.remote_timeout
conf.remote_timeout = 10000
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris
solar_system_ephemeris.set("jpl")

from poliastro.bodies import Sun, Earth, Moon, Venus, Mercury
#from poliastro.ephem import Ephem
from poliastro.frames import Planes
from poliastro.twobody import Orbit
from poliastro.plotting import StaticOrbitPlotter
from poliastro.plotting.misc import plot_solar_system
from poliastro.plotting import OrbitPlotter2D
from poliastro.util import time_range
import matplotlib.pyplot as op
from plotly.graph_objs import Figure, FigureWidget
fig = FigureWidget()
op = OrbitPlotter2D(fig)

EPOCH = Time("2026-09-01 12:00:00", scale="tdb")
EPOCH_1 = Time("2026-12-01 12:00:00", scale="tdb")
C_FLORENCE = "#000"
C_MOON = "#999" 

#Earth.plot(EPOCH);

florence_osc = Orbit.from_sbdb("65803 Didymos")
#florence_os = Orbit.from_sbdb("65803 Didymos")
earth = Orbit.from_body_ephem(Earth, EPOCH_1)
venus = Orbit.from_body_ephem(Venus, EPOCH_1)
mercury = Orbit.from_body_ephem(Mercury, EPOCH_1)
florence_osc
florence_osc.epoch.iso
op.plot(florence_osc, label="65803 Didymos",color="purple")
#op.plot(florence_os, label="65803 Didymos",color="purple")
op.plot(earth, label="Earth",color="blue")
op.plot(venus, label="Venus", color="orange")
op.plot(mercury, label="Mercury", color="grey")
#op.plot(Earth,label="Initial orbit")

florence = florence_osc.propagate(EPOCH_1)

from astropy.coordinates import (
    ICRS, GCRS,
    CartesianRepresentation, CartesianDifferential
)
from poliastro.frames import HeliocentricEclipticJ2000
florence_heclip = HeliocentricEclipticJ2000(
    x=florence.r[0], y=florence.r[1], z=florence.r[2],
    v_x=florence.v[0], v_y=florence.v[1], v_z=florence.v[2],
    representation=CartesianRepresentation,
    differential_type=CartesianDifferential,
    obstime=EPOCH_1
)
print(florence_heclip)


florence_icrs_trans = florence_heclip.transform_to(ICRS)
florence_icrs_trans.representation = CartesianRepresentation
print(florence_icrs_trans)

florence_icrs = Orbit.from_vectors(
    Sun,
    r=[florence_icrs_trans.x, florence_icrs_trans.y, florence_icrs_trans.z] * u.km,
    v=[florence_icrs_trans.v_x, florence_icrs_trans.v_y, florence_icrs_trans.v_z] * (u.km / u.s),
    epoch=florence.epoch
)


op.plot(florence_icrs, label="65803 Didymos",color="purple")
fig.show()

## September
#<HeliocentricEclipticIAU76 Coordinate (obstime=2026-09-01 12:00:00.000): (x, y, z) in km
#    (67969641.02827488, -2.23652971e+08, -7722639.23892914)
# (v_x, v_y, v_z) in km / s
#    (18.8488624, 15.47623837, -0.80828651)>
#<ICRS Coordinate: (x, y, z) in km
#    (67756463.29392487, -2.02846764e+08, -96346002.51123594)
# (v_x, v_y, v_z) in km / s
#    (18.8488624, 14.52068901, 5.41450571)>
## December
#<HeliocentricEclipticIAU76 Coordinate (obstime=2026-12-01 12:00): (x, y, z) in km
#    (1.60778844e+08, -27515205.2200428, -9639133.51806481)
# (v_x, v_y, v_z) in km / s
#    (-0.91998542, 32.91995336, 0.61881203)>
#<ICRS Coordinate: (x, y, z) in km
#    (1.60647334e+08, -22097757.09316747, -20072491.72020712)
# (v_x, v_y, v_z) in km / s
#    (-0.91998542, 29.9573174, 13.66255436)>
