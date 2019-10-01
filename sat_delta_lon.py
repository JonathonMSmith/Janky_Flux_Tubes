'''script to get lat lon alt for flux tubes connected to ground stations
'''
import csv
import pickle
import os
import numpy as np
import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from mayavi.sources.builtin_surface import BuiltinSurface
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Rectangle
import sami2py

def spherical(lats, lons, radii):
    '''inputs are arrays not scalars... although i guess it should work w/ both
       alts should be added to earth radius before this function
    '''
    #lat and lon to angle and azimuth
    pol = (90 - lats) * np.pi / 180
    az = lons * np.pi / 180
    X = radii * np.cos(az) * np.sin(pol)
    Y = radii * np.sin(az) * np.sin(pol)
    Z = radii * np.cos(pol)
    return X, Y, Z


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

    x, y, z = spherical(mag_eq[1], mag_eq[0], 6371)
    sf = 950 #scale_factor
    mlab.plot3d(x/sf, y/sf, z/sf, color=(0, 0, 0))

plot_3D = True
observatories = {'Haystack':(42.64, -71.45), 'McDonald':(30.67, -104.02),
                 'Arecibo':(18.3, -66.8), 'Zaquencipa':(5.6, -73.52),
                 'Jicamarca':(-11.95, -76.87), 'El Leoncito':(-31.8, -69.3),
                 'Mercedes':(-34.51, -59.4), 'SAAMER':(-53.79, -67.75),
                 'Mount John':(-43.99, 170.46), 'SAO':(-33.93, 18.48),
                 'Asiago':(45.87, 11.53)}

color_list = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 0, 0),
              (.5, .5, .5), (1, 0, 1), (0, 1, 1), (.6, .6, .2), (.2, .6, .6),
              (.6, .2, .6)]
flux_tube_dict = {}

if plot_3D:
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    max_1 = []
    max_2 = []
    R_cont = 6.5
    #make an earth
    sphere = mlab.points3d(0, 0, 0, scale_factor=13, color=(.67, .77, .93),
                           resolution=50, opacity=1, name="Earth")
    sphere = mlab.points3d(0, 0, 0, scale_factor=13*(6771/6371),
                           color=(.1, .1, .1), resolution=50,
                           opacity=.35, name="Orbit")
    Continents_src = BuiltinSurface(source="earth", name="Continents")
    continents = mlab.pipeline.surface(Continents_src, color=(0, 0, 0),
                                       extent=[-R_cont+.5, R_cont+.5, -R_cont,
                                               R_cont, -R_cont, R_cont])
    draw_mag_equator()
n = 0
for observatory in observatories:
    lon = observatories[observatory][1]
    lat = observatories[observatory][0]
    print(observatory+'::lat:'+str(lat)+'; lon:'+str(lon)+'; color:'+str(color_list[n]))
    flux_tube_dict[lon] = {}
    path = sami2py.utils.generate_path('petit', year=2012, lon=int(lon), day=0)
    if not os.path.exists(path):
        sami2py.run_model(tag='petit', day=0, year=2012, lon=int(lon), hrpr=.0,
                          dthr=0.05, hrinit=0.0, hrmax=.11, rmax=18000)

    S = sami2py.Model(tag='petit', day=0, year=2012, lon=int(lon))
    nf = np.shape(S.zalt)[1]
    l_400_1 = []
    l_400_2 = []
    l_200_1 = []
    l_200_2 = []
    for ft in range(nf):
        alts = S.zalt[:, ft]
        lons = S.glon[:, ft]
        lats = S.glat[:, ft]
        if not plot_3D:
            plt.plot(lats, alts)
        # find locations along flux tube within 10 degrees latitude of ground
        # based instruments below 251 km, or locations that the ground based
        # instruments can see.
        inds, = np.where((lats > lat-10) & (lats < lat+10) & (np.round(alts) < 251))
        # does the flux tube have an apex at or above 400 km? this is
        # necessary for the satellite to observe plasma on the flux tube
        has_400km, = np.where(alts >= 400)
        # if the flux tube doesn't connect with the ground station or the
        # satellite then continue
        if inds.size == 0 or has_400km.size == 0:
            continue
        #do some longitude shifting for calculations
        if np.max(lons) - np.min(lons) > 200:
            lons = [x+360 if x < 300 else x for x in lons]

        #make a dictionary for the flux tube to store all the stuff... why tho
        flux_tube_dict[lon][ft] = {}
        flux_tube_dict[lon][ft]['alts'] = alts
        flux_tube_dict[lon][ft]['lats'] = lats
        flux_tube_dict[lon][ft]['lons'] = lons

        #get the radial distance of the sat in... km?
        r = alts + 6371# / 2000

        #lat and lon to angle and azimuth
        theta = (90 - lats) * np.pi / 180
        phi = lons * np.pi / 180
        x = r * np.cos(phi) * np.sin(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(theta)

        #plot the flux tube
        if plot_3D:
            sf = 1000 #scale_factor
            mlab.plot3d(x/sf, y/sf, z/sf, color=color_list[n], line_width=1000)
        #get lons above 200 and 400 and then store the two values at the most
        #extreme latitudes
        lons_400 = [x for x, y in zip(lons, alts) if y > 400]
        lons_200 = [x for x, y in zip(lons, alts) if y > 200]
        l_400_1.append(lons_400[0])
        l_400_2.append(lons_400[-1])
        l_200_1.append(lons_200[0])
        l_200_2.append(lons_200[-1])
    n += 1
    dif_1 = []
    dif_2 = []
    for i in range(len(l_400_1)):
        dif_1.append(l_400_1[i] - l_200_2[i])
        dif_2.append(l_200_1[i] - l_400_2[i])
    max_1.append(np.max(np.abs(dif_1)))
    max_2.append(np.max(np.abs(dif_2)))
    print(np.max(np.abs(dif_1)))
    print(np.max(np.abs(dif_2)))

sat_delta_lon = {'observatories': observatories, 'flux_tubes': flux_tube_dict,
                 'delta_1': max_1, 'delta2': max_2}
pickle.dump(sat_delta_lon, open('sat_delta_lon.p', 'wb'))
if plot_3D:
    mlab.axes(extent=[0, 360, -80, 80, 0, 5])
    mlab.show()
else:
    plt.show()
