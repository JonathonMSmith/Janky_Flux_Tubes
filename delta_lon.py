'''Field line tracing for delta lons for petitsat acs
        an observatory in the Southern Hemisphere on a field line with a
        negative declination will be in darkness longer than the conjugate
        point on the dusk-side. 
        Summary:
        So. hem. + neg. dec. = dusk and positive delta glon
        So. hem. + pos. dec. = dawn and negative delta glon
        No. hem. + neg. dec. = dawn and negative delta glon
        No. hem. + pos. dec. = dusk and positive delta glon
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

def get_delta_lon(field_line, alt_1, alt_2):
    '''get the delta longitude between two altitudes on one field line
       Parameters:
       field_line: field line object from pysatMagVect
       alt_1: in this case a dummy altitude for the air-glow layer
       alt_2: the altitude at which to find the conjugate point, in this case
              roughly the ISS orbit altitude
    '''
    n = 0
    first_point = np.full([1, 3], np.nan)[0]
    critical_point = np.full([1, 3], np.nan)[0]
    for point in field_line:
        glat, glon, zalt = pmv.ecef_to_geocentric(point[0], point[1], point[2])
        if zalt >= alt_1 and n == 0:
            n += 1
            first_point = [glat, glon, zalt]
        if zalt >= alt_2 and np.abs(glat) < 52:
            critical_point = [glat, glon, zalt]
    delta_lon = first_point[1] - critical_point[1]
    if np.abs(delta_lon) > 180:
        delta_lon = (delta_lon / np.abs(delta_lon)) * (np.abs(delta_lon) - 360)
    print(critical_point[0])
    return delta_lon

def draw_mag_equator():
    '''Draws the magnetic equator using mayavi and some igrf grid output
    '''
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

def draw_picture():
    '''Draws a rough picture of a globe with continents, the magnetic equator,
        a transparent sphere at the 'orbit' altitude and the field lines
        connected to observatories
    '''
    R_cont = 1
    #make an earth
    erf = mlab.points3d(0, 0, 0, scale_factor=2, color=(.67, .77, .93),
                        resolution=50, opacity=1, name="Earth")
    orb = mlab.points3d(0, 0, 0, scale_factor=2*(6771/6371),
                        color=(.1, .1, .1), resolution=50,
                        opacity=.35, name="Orbit")
    Continents_src = BuiltinSurface(source="earth", name="Continents")
    continents = mlab.pipeline.surface(Continents_src, color=(0, 0, 0))
    draw_mag_equator()

    date = datetime.datetime(2015, 1, 1)
    for observatory in OBSERVATORIES:
        lon = OBSERVATORIES[observatory][1]
        lat = OBSERVATORIES[observatory][0]
        x, y, z = pmv.geocentric_to_ecef(lat, lon, 280)

        if lat < 0:
            direction = 1
        else:
            direction = -1

        field_line = pmv.field_line_trace([x, y, z], date=date, height=10,
                                          direction=direction)
        print(get_delta_lon(field_line, 300, 400))
        print(observatory)
        sf = 6371 #scale_factor
        mlab.plot3d(field_line[:, 0]/sf, field_line[:, 1]/sf,
                    field_line[:, 2]/sf)
    mlab.show()

def write_delta_info():
    '''Record some information about the observatories and their field lines.
        importantly the delta longitude between the conjugate orbit altitude.
        Writes all of this to a dump file. The dump file will need to be
        cleared manually between runs for the time being or it will keep
        appending lines over and again.
    '''
    date = datetime.datetime(2015, 1, 1)
    for observatory in OBSERVATORIES:
        lon = OBSERVATORIES[observatory][1]
        lat = OBSERVATORIES[observatory][0]
        x, y, z = pmv.geocentric_to_ecef(lat, lon, 280)

        if lat < 0:
            direction = 1
            hemisphere = 'Southern'
        else:
            direction = -1
            hemisphere = 'Northern'

        field_line = pmv.field_line_trace([x, y, z], date=date, height=10,
                                          direction=direction)
        delta_lon = get_delta_lon(field_line, 300, 400)

        if hemisphere == 'Northern' and delta_lon > 0:
            declination = 'positive'
            terminator = 'dusk'
        elif hemisphere == 'Northern' and delta_lon < 0:
            declination = 'negative'
            terminator = 'dawn'
        elif hemisphere == 'Southern' and delta_lon > 0:
            declination = 'negative'
            terminator = 'dusk'
        elif hemisphere == 'Southern' and delta_lon < 0:
            declination = 'positive'
            terminator = 'dawn'
        else:
            declination = 'nan'
            terminator = 'nan'
        output = ('Observatory:' + observatory +
                  '; Hemisphere:' + hemisphere +
                  '; Declination:' + declination +
                  '; Delta_Lon:' + str(np.ceil(delta_lon)) +
                  '; Terminator:' + terminator)
        dump_file = open('dump_file.txt', 'a')
        dump_file.write(output + '\n')
        dump_file.close()
