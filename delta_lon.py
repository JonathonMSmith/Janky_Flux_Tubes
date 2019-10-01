'''Field line tracing for delta lons for petitsat acs
'''
import datetime
import csv
import numpy as np
import pysatMagVect as pmv
#import matplotlib.pyplot as plt
import mayavi.mlab as mlab
from mayavi.sources.builtin_surface import BuiltinSurface

OBSERVATORIES = {'Haystack':(42.64, -71.45), 'McDonald':(30.67, -104.02),
                 'Arecibo':(18.3, -66.8), 'Zaquencipa':(5.6, -73.52),
                 'Jicamarca':(-11.95, -76.87), 'El Leoncito':(-31.8, -69.3),
                 'Mercedes':(-34.51, -59.4), 'SAAMER':(-53.79, -67.75),
                 'Mount John':(-43.99, 170.46), 'SAO':(-33.93, 18.48),
                 'Asiago':(45.87, 11.53)}


def draw_mag_equator():
    inc = []
    lat = []
    lon = []
    cnt = 0
    with open('igrfgridData.csv') as csvfile:
        igrf = csv.reader(csvfile, delimiter=',')
        for row in igrf:
            cnt += 1
            if cnt < 14:
                continue
    
            inc.append(np.abs(float(row[4])))
            lat.append(float(row[1]))
            lon.append(float(row[2]))
    inc = np.array(inc)
    lat = np.array(lat)
    lon = np.array(lon)
    
    idx, = np.where(inc < .25)
    lat = lat[idx]
    lon = lon[idx]

    mag_eq = np.array([lon, lat, inc])
    i = np.argsort(mag_eq[0])
    mag_eq = np.array([mag_eq[0][i], mag_eq[1][i]])

    x, y, z = pmv.geocentric_to_ecef(mag_eq[1], mag_eq[0], 0)
    sf = 6371 #scale_factor
    mlab.plot3d(x/sf, y/sf, z/sf, color=(0, 0, 0))

R_cont = 1 
#make an earth
sphere = mlab.points3d(0, 0, 0, scale_factor=2, color=(.67, .77, .93),
                       resolution=50, opacity=1, name="Earth")
sphere = mlab.points3d(0, 0, 0, scale_factor=2*(6771/6371),
                       color=(.1, .1, .1), resolution=50,
                       opacity=.35, name="Orbit")
Continents_src = BuiltinSurface(source="earth", name="Continents")
continents = mlab.pipeline.surface(Continents_src, color=(0, 0, 0))
                                   #extent=[-R_cont, R_cont, -R_cont,
                                           #R_cont, -R_cont, R_cont])
draw_mag_equator()

date = datetime.datetime(2015, 1, 1)
for observatory in OBSERVATORIES:
    print(observatory)
    lon = OBSERVATORIES[observatory][1]
    lat = OBSERVATORIES[observatory][0]
    x, y, z = pmv.geocentric_to_ecef(lat, lon, 100)

    if lat < 0:
        direction = 1
    else:
        direction = -1

    field_line = pmv.field_line_trace([x, y, z], date=date, height=10,
                                      direction=direction)

#    ax.plot3D(field_line[:, 2], field_line[:, 0], field_line[:, 1], 'k')

    sf = 6371 #scale_factor
    mlab.plot3d(field_line[:, 0]/sf, field_line[:, 1]/sf,
                field_line[:, 2]/sf)

#ax.set_xlim3d(-18000, 18000)
#ax.set_ylim3d(-18000, 18000)
#ax.set_zlim3d(-18000, 18000)
#plt.show()
mlab.show()
