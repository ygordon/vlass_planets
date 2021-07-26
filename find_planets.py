###find which QL subtiles have planets in using Greg Sivakoff's code

from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body, get_moon
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np, argparse

from astropy.table import Table

########################################
#time of code check - not important for code to run
import time
start_time = time.time()
########################################


####################################################################################
####################################################################################
###params

#subtile_file = ('/Users/yjangordon/Documents/science/survey_data/VLASS/catalogues/'
#                + 'CIRADA_VLASS1QL_table3_subtile_info_v1a_correctdates.fits')
#comp_file = ('/Users/yjangordon/Documents/science/survey_data/VLASS/catalogues/'
#             + 'CIRADA_VLASS1QL_table1_components_v1a.fits')


####greg's code test
#test_time = Time("2014-09-22 23:22")
#test_coordinates = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')

jpos1 = SkyCoord('14h58m32.73s-15d50m45s')
jtim1 = Time('2018-01-01 14:13:41.925', scale='utc')

#jpos2 = SkyCoord('15h17m20.7s-17d02m18s')
#jtim2 = Time('2018-03-02 16:23:28.341', scale='utc')

j1subtile = 'J145838-153000'


eloc = EarthLocation.of_site('vla')
bodies = ('sun', 'mercury', 'venus', 'moon', 'mars', 'jupiter',
          'saturn', 'uranus', 'neptune', 'pluto')




####################################################################################
####################################################################################
##functions


#planet_warning = False
#with solar_system_ephemeris.set('de432s'):
#    for test_body in test_bodies:
#    body_coords = get_body(test_body, test_time, test_loc)
#    sep = test_coordinates.separation(body_coords)
#    if sep.arcsec < 60:
#        print("Warning: {} within {} arcseconds".format(test_body.capitalize(),np.ceil(sep.arcsec)))
#        planet_warning = True
#    break
#
#if planet_warning == False:
#    print ("All Clear")



def planet_find(coord, obstime, maxsep=30*u.arcmin):
    'use emphemerids to determine if a planet/sun/moon is potentially in field of obs'
    
#    t0 = time.time()

    eloc = EarthLocation.of_site('vla')
    bodies = ['sun', 'mercury', 'venus', 'moon', 'mars', 'jupiter',
              'saturn', 'uranus', 'neptune', 'pluto']

    planet_warning = False
    planet_warnings = {}
    with solar_system_ephemeris.set('de432s'):
        for body in bodies:
            body_coords = get_body(body, obstime, eloc)
            sep = body_coords.separation(coord) ###has to be body_coord.sep(coord) not coord.sep(body_coord), the latter for some reason gives a wrong result (probably the ordering and complexity of body_coords not being accounted for)
            if sep < maxsep:
                print(" ")
                print("Warning: {} within {} arcsec".format(body.capitalize(),np.ceil(sep.arcsec)))
                print(" ")
#                planet_warning = True
                planet_warnings[body] = SkyCoord(body_coords.ra, body_coords.dec)

    ###how long to run
#    t1 = time.time()
#    print('')
#    print('time = ', np.round(t1-t0, 3), 's')

    if len(planet_warnings.keys())>0:
        return planet_warnings


def search_subtiles_for_planets(data,
                                acol='OBSRA',
                                dcol='OBSDEC',
                                tcol='DATE-OBS',
                                ncol='Subtile',
                                search_radius=43*u.arcmin):
    'use 43 arcmin radius (subtiles are 60*60arcmin square) -- may lead to some false positives'
    
    coords = SkyCoord(ra=data[acol], dec=data[dcol])
    times = Time(data[tcol])
    names = data[ncol]
    
    bware = {}
    
    
    ###set up print statements to show progress through data
    t0 = time.time() #sets start time
    ni = len(data)/100
    ni = 10*(int(ni/10)) ##rounds to nearest 10
    if ni < 1:
        ni = 1 ###accounts for smaller data sets
    ti = np.arange(0, len(data), ni)

    
    for i in range(len(data)):
        co = coords[i]
        tm = times[i]
        nm = names[i]
        pf = planet_find(coord=co, obstime=tm, maxsep=search_radius)
        if pf is not None:
            bware[nm] = pf
        ###show progress
        if i in ti:
            print(i+1,'/',len(data), '   -   ', np.round(time.time()-t0, 2), 's')
    
    return bware


def planet_search_to_table(psdict):
    'create astropy table to output from dict returned by search_subtiles..'
    nkeys = list(psdict.keys())
    
    stiles = []
    planets = []
    plan_a = []
    plan_d = []
    
    for key in nkeys:
        pd = psdict[key]
        plans = list(pd.keys())
        coords = list(pd.values())
        for i in range(len(pd)):
            stiles.append(key)
            planets.append(plans[i].capitalize())
            plan_a.append(coords[i].ra.value)
            plan_d.append(coords[i].dec.value)

    outtab = {'Planet': planets,
              'RA': np.array(plan_a)*u.deg,
              'DEC': np.array(plan_d)*u.deg,
              'Subtile': stiles}
    
    return Table(outtab)


def parseargs():
    parser = argparse.ArgumentParser(description="command line arguments")
    parser.add_argument('obsfile', help='file containing observation info to check for planets')
    
    args = parser.parse_args()
    
    return args


####################################################################################
####################################################################################
###main

#subtiles = Table.read(subtile_file)
#
#test_st = subtiles[22011-5: 22011+5]
#
#
#tsearch = search_subtiles_for_planets(test_st)
#
#tout = planet_search_to_table(tsearch)


if __name__ == '__main__':
    args = parseargs()
    data = Table.read(args.obsfile)
    searchres = search_subtiles_for_planets(data)
    outdata = planet_search_to_table(searchres)
    outdata.write('planet_finding_results.fits', format='fits')
    print('END: ', ("--- %s seconds ---" % (time.time() - start_time)))




#########################################
#time taken for code to run
#print('END: ', ("--- %s seconds ---" % (time.time() - start_time)))


