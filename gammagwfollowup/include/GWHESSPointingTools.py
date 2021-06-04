#####################################################################
# Author: Monica Seglar-Arroyo
# Contributors: Halim Ashkar,  Fabian Schussler
# All the tools that are needed to follow-up a GW with an IACT (HESS)
# are described and implemented in the following.
#####################################################################
# Packages
# Packages

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from astropy.table import Table
from astropy import units as u
from astropy.utils.data import download_file
from astropy.utils import iers
import astropy.coordinates as co
from astropy.time import Time
from astropy.coordinates import (SkyCoord, EarthLocation, AltAz,
                                 get_sun, get_moon)
from astropy.io import fits
import datetime
import ephem
from mocpy import MOC
import numpy.ma as ma
from scipy.stats import norm
import time
import os

#iers_url_mirror='ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
#IERS_A_URL_MIRROR = 'https://datacenter.iers.org/data/9/finals2000A.all'
#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers.IERS_A_URL, cache=True))

# iers_file = './gammagwfollowup/finals2000A.all'
iers_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../dataset/finals2000A.all')
iers.IERS.iers_table = iers.IERS_A.open(iers_file)
#iers.IERS.iers_table = iers.IERS_A.open(download_file(iers_url_mirror, cache=True))

#Darkness criteria: max sun and moon altitude in degrees.

gSunDown = -18
#HorizonSun = '-17:43:48'
HorizonSun = '-18:00:00'
gMoonDown = -0.5
HorizonMoon = '-00:30:00'
gMoonGrey=50  #Altitude in degrees
gMoonPhase=60 #Phase in %
MoonSourceSeparation=30 #Separation in degrees
# The twilight definitions specify the position of the center of the Sun,
# while the rising and setting functions normally pay attention to the top edge of an object instead!



######################################################

##                      Classes                     ##

######################################################
class HESSObservatory:
    def __init__(self):
        self.Lat = -23.271778 * u.deg
        self.Lon = 16.50022 * u.deg
        self.Height = 1835 * u.m
        self.Location = EarthLocation(lat=self.Lat, lon=self.Lon,
                                      height=self.Height)

class CTASouthObservatory:

    def __init__(self):
        self.Name='South'
        self.Lat = -24.5 * u.deg
        self.Lon = -70.17 * u.deg
        self.Height = 2635 * u.m
        self.Location = EarthLocation(lat=self.Lat, lon=self.Lon,
                                      height=self.Height)


class CTANorthObservatory:

    def __init__(self):
        self.Name='North'
        self.Lat = 28.75 * u.deg
        self.Lon = 17.5 * u.deg
        self.Height = 2200 * u.m
        self.Location = EarthLocation(lat=self.Lat, lon=self.Lon,
                                      height=self.Height)


class Tools:
    '''
        class with different visibility check functions and other setting and rising of the sun and the moon functions
        '''

    @classmethod
    def IsDarkness(cls, obsTime,obsSite):
        sunAlt = Tools.SunAlt(obsTime,obsSite)
        moonAlt = Tools.MoonAlt(obsTime,obsSite)
        #print(obsTime,sunAlt,moonAlt)
        if sunAlt > gSunDown:
            return False
        if moonAlt > gMoonDown:
            return False

        return True

    @classmethod
    def IsGreyness(cls, obsTime,obsSite):
        # SUN altitude
        sunAlt = Tools.SunAlt(obsTime,obsSite)
        # MOON altitude
        moonAlt = Tools.MoonAlt(obsTime,obsSite)
        # MOON azimuth
        #moonAz = Tools.MoonAz(obsTime)
        # MOON phase
        moonPhase = Tools.MoonPhase(obsTime)

        if sunAlt > gSunDown:
            return False
        if moonAlt > gMoonGrey:
            return False
        if moonPhase > gMoonPhase and moonAlt > gMoonDown:
            return False
        return True

    @classmethod
    def MoonPhase(cls, obsTime):
        hess = HESSObservatory()
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(hess.Lon / u.deg)
        obs.lat = str(hess.Lat / u.deg)
        obs.elev = hess.Height / u.m
        obs.date = obsTime
        moon.compute(obs)

        #print("Phase of the moon = %s percent" % moon.phase)

        return moon.phase


    @classmethod
    def SunAlt(cls, obsTime,obsSite):
        sun = get_sun(Time(obsTime)).transform_to(AltAz(obstime=Time(obsTime),
                                                        location=obsSite.Location))
        #print(Time(obsTime),obsSite.Location)
        #print(get_sun(Time(obsTime)))
        #print(sun.alt/u.deg)
        return sun.alt / u.deg

    @classmethod
    def MoonAlt(cls, obsTime,obsSite):  # THE ERROR IS IN THIS FUNCTION!!!!
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        #print(obs)
        moon.compute(obs)
        #print('Altitude of the moon = ',moon.alt * 180. / np.pi)
        return moon.alt * 180. / np.pi

    @classmethod
    def MoonAz(cls, obsTime):
        hess = HESSObservatory()
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(hess.Lon / u.deg)
        obs.lat = str(hess.Lat / u.deg)
        obs.elev = hess.Height / u.m
        obs.date = obsTime
        moon.compute(obs)
        print('Azimuth of the moon = ',moon.az * 180. / np.pi)
        return moon.az * 180 / np.pi

    @classmethod
    def NextSunrise(cls, obsTime,obsSite):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonSun
        sun.compute(obs)
        nextSunrise = obs.next_rising(sun,use_center=True).datetime()
        return nextSunrise

    @classmethod
    def PreviousSunrise(cls, obsTime,obsSite):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonSun
        sun.compute(obs)
        previousSunrise = obs.previous_rising(sun,use_center=True).datetime()
        return previousSunrise

    @classmethod
    def NextSunset(cls, obsTime,obsSite):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonSun
        sun.compute(obs)
        nextSunset = obs.next_setting(sun,use_center=True).datetime()
        return nextSunset

    @classmethod
    def PreviousMoonrise(cls, obsTime,obsSite):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonMoon
        # print('NewHorizon=',obs.horizon)
        moon.compute()
        previousMoonrise = obs.previous_rising(moon,use_center=True).datetime()
        return previousMoonrise

    @classmethod
    def NextMoonrise(cls, obsTime,obsSite):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonMoon
        moon.compute()
        nextMoonrise = obs.next_rising(moon,use_center=True).datetime()
        return nextMoonrise

    @classmethod
    def PreviousSunset(cls, obsTime,obsSite):
        sun = ephem.Sun()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonSun
        sun.compute()
        previousSunset = obs.previous_setting(sun,use_center=True).datetime()
        return previousSunset

    @classmethod
    def PreviousMoonset(cls, obsTime,obsSite):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonMoon
        moon.compute()
        previousMoonset = obs.previous_setting(moon,use_center=True).datetime()
        return previousMoonset

    @classmethod
    def NextMoonset(cls, obsTime,obsSite):
        moon = ephem.Moon()
        obs = ephem.Observer()
        obs.lon = str(obsSite.Lon / u.deg)
        obs.lat = str(obsSite.Lat / u.deg)
        obs.elev = obsSite.Height / u.m
        obs.date = obsTime
        obs.horizon = HorizonMoon
        moon.compute()
        nextMoonset = obs.next_setting(moon).datetime()
        #print('NextMoonset',nextMoonset)
        previousMoonset = obs.previous_setting(moon,use_center=True).datetime()
        #print(previousMoonset)
        return nextMoonset

    @classmethod
    def TrustingDarknessSun(cls, obsTime,obsSite):
        DarkObsTime = obsTime
        referencetime=obsTime
        while (Tools.IsDarkness(DarkObsTime,obsSite) == False and ((DarkObsTime.hour >= referencetime.hour and DarkObsTime.day == referencetime.day) or (
                DarkObsTime.hour <= Tools.NextSunrise(referencetime,obsSite).hour and DarkObsTime.day == Tools.NextSunrise(
                referencetime,obsSite).day))):
            DarkObsTime = DarkObsTime + datetime.timedelta(minutes=1)
            #print(Tools.IsDarkness(DarkObsTime,obsSite))
        return DarkObsTime

    @classmethod
    def TrustingGreynessSun(cls, obsTime,obsSite):
        GreyObsTime = obsTime
        referencetime=obsTime
        while (Tools.IsGreyness(GreyObsTime,obsSite) == False and ((GreyObsTime.hour >= referencetime.hour and GreyObsTime.day == referencetime.day) or (
                GreyObsTime.hour <= Tools.NextSunrise(referencetime,obsSite).hour and GreyObsTime.day == Tools.NextSunrise(
                referencetime,obsSite).day))):
            GreyObsTime = GreyObsTime + datetime.timedelta(minutes=1)
            #print(Tools.IsDarkness(DarkObsTime,obsSite))
        return GreyObsTime

    @classmethod
    def TrustingDarknessMoon(cls, obsTime, referencetime,obsSite):
        DarkObsTime = obsTime
        # Make sure that its night
        if ((DarkObsTime.hour >= referencetime.hour and DarkObsTime.day == referencetime.day) or (
                DarkObsTime.hour <= Tools.NextSunrise(referencetime,obsSite).hour and DarkObsTime.day == Tools.NextSunrise(
                referencetime,obsSite).day)):
            while (Tools.IsDarkness(DarkObsTime,obsSite) == False):
                DarkObsTime = DarkObsTime + datetime.timedelta(minutes=1)
        return DarkObsTime

    @classmethod
    def TrustingGreynessMoon(cls, obsTime, referencetime,obsSite):
        GreyObsTime = obsTime
        # Make sure that its night
        if ((GreyObsTime.hour >= referencetime.hour and GreyObsTime.day == referencetime.day) or (
                GreyObsTime.hour <= Tools.NextSunrise(referencetime,obsSite).hour and GreyObsTime.day == Tools.NextSunrise(
                referencetime,obsSite).day)):
            while (Tools.IsGreyness(GreyObsTime,obsSite) == True):
                GreyObsTime = GreyObsTime + datetime.timedelta(minutes=1)
        if (Tools.IsGreyness(GreyObsTime,obsSite) == False):
            if ((GreyObsTime - obsTime)>=datetime.timedelta(minutes=10)):
                return True, GreyObsTime
            if ((GreyObsTime - obsTime)<datetime.timedelta(minutes=10)):
                return False, GreyObsTime


    @classmethod
    def UTCtoNamibia(cls, UTCtime):
        TimezonesDifference = datetime.timedelta(hours=2)
        NamibianTime = UTCtime + TimezonesDifference
        return NamibianTime

    @classmethod
    def NextObservationWindow(cls,time, obsSite):
        if(Tools.NextSunset(time, obsSite).hour >= time.hour >= Tools.PreviousSunrise(time,obsSite).hour and time.day == Tools.NextSunset(time, obsSite).day):
            time = Tools.NextSunset(time, obsSite)
            # print('Sunset', time)
            time = Tools.TrustingGreynessSun(time, obsSite)
            # print('Trusted', time)
        if(Tools.IsGreyness(time, obsSite) is True):
            return time
        elif ((Tools.IsGreyness(time, obsSite) is False)):
            time = Tools.TrustingGreynessSun(time, obsSite)
            #time=Tools.NextMoonset(time, obsSite)
            return time
        else:
            print('No window is found')
            return False

    @classmethod
    def CheckWindow(cls,time, obsSite):
        WindowDuration = datetime.timedelta(minutes=28)
        MinimalWindowDuration = datetime.timedelta(minutes=10)
        if (Tools.IsDarkness(time, obsSite) is True) and (Tools.IsDarkness(time + MinimalWindowDuration, obsSite) is True):
            Observe = True
        else:
            print('No window found')
            Observe = False

        return Observe
    @classmethod
    def GalacticPlaneBorder(cls,coords):
        lon=coords.galactic.l.value #x-coordinate
        lat=coords.galactic.b.value  #y-coordinate
        print(lon)
        print(lat)
        YouAreInside = False
        n=20
        if(lat<=10 and lat>=0 and lon<=130):
            n=lat-(1.0/13)*lon
        elif (lat <= 10 and lat >= 0 and lon>=240):
            n = lat - (1.0/12) * lon +20
        elif (lat >= -10 and lat <= 0 and lon <= 130):
            n = lat + (1.0/13) * lon
        elif (lat>= -10 and lat <= 0 and lon>=240):
            n = lat +(1.0/12) * lon - 20
        print(n)
        if np.absolute(n)<=10:
            YouAreInside=True
            #print('You got here')
        print(YouAreInside)
        return YouAreInside


######################################################

## Functions from BestCandidateon PGW ==> 2D

######################################################

def LoadHealpixMap(thisfilename):
    '''Download aLIGO HEALpix map and keep in cache
        RETURNS:
        --------

        tprob : array of p-values as a function of sky position
        tdistmu : array of distance estimate
        tdistsigma : array of error on distance estimates
        distnorm : array of distance normalisations
        detectors: which interferometers triggered
        event_id: ID of the event
        distmean: mean distance from the header
        disterr: error on distance from the header
        '''
    PrintFileName = "Loading LVC HEALPix map from file: " + thisfilename
    print(PrintFileName)
    fitsfile = fits.open(thisfilename)

    tevent_id = "Non specified"
    tdetectors = ""
    tdistmean = 0
    tdisterr = 0
    tdistmu = []
    tdistsigma = []
    tdistnorm = []

    if 'OBJECT' in fitsfile[1].header:
        tevent_id = fitsfile[1].header['OBJECT']
    else:
        tevent_id = "Non specified"

    if 'INSTRUME' in fitsfile[1].header:
        tdetectors = fitsfile[1].header['INSTRUME']
    else:
        tdetectors = "Non specified"

    if (fitsfile[1].header['TFIELDS'] == 4):
        tprob, tdistmu, tdistsigma, tdistnorm = hp.read_map(thisfilename, field=range(4))
        tdistmean = fitsfile[1].header['DISTMEAN']
        tdisterr = fitsfile[1].header['DISTSTD']
        print('Event has triggered ', tdetectors, ' => distance = {0:.2f}'.format(tdistmean),' +- {0:.2f}'.format(tdisterr),' Mpc')
    else:
        tprob = hp.read_map(thisfilename, field=range(1))
    # raise

    fitsfile.close()

    return tprob, tdistmu, tdistsigma, tdistnorm, tdetectors, tevent_id, tdistmean, tdisterr

def NightDarkObservation(time, obsSite, MaxNights,duration,minduration):
    '''
    Function that searches for an array of observation times that fulfilled darkness condition and window

    '''
    WindowDuration = datetime.timedelta(minutes=duration)
    MinimalWindowDuration = datetime.timedelta(minutes=minduration)
    AuxMax = 100
    NightDarkRuns = []
    isFirstNight = True
    # Loop for the nights of observation
    for i in range(0, MaxNights):
        # In case the alert arrives during the day, time is set to sunset
        if (Tools.NextSunset(time,obsSite).hour >= time.hour >= Tools.PreviousSunrise(
                time,obsSite).hour and time.day == Tools.NextSunset(time,obsSite).day):
            time = Tools.NextSunset(time,obsSite)
            #print('POSTTIME',time,'isdarks', Tools.IsDarkness(time,obsSite))
            time = Tools.TrustingDarknessSun(time,obsSite)
            #print('POSTPOSTTIME', time, 'isdarks', Tools.IsDarkness(time,obsSite))
            isFirstNight = False
        # in case the alert arrives during the night, the last observation night will be the fourth so MaxNights=4
        else:
            if (isFirstNight):
                #MaxNights = MaxNights+1
                MaxNights = MaxNights
                isFirstNight = False
        # Using Time 0 in order to define what would be a night through that goes from Time 0 to NextSunrise(Time 0)
        time0 = time
        for j in range(0, AuxMax):
            # Night condition fulfilled
            if ((time.hour >= time0.hour and time.day == time0.day) or (
                    time.hour <= Tools.NextSunrise(time0,obsSite).hour and time.day == Tools.NextSunrise(time0,obsSite).day)):
                # print('time.hour >= time0.hour and time.day == time0.day) or (time.hour <= Tools.NextSunrise(time0).hour and time.day == Tools.NextSunrise(time0).day')
                # print(time.hour ,'   ', time0.hour,'    ', time.day ,'   ',time0.day ,'   ',time.hour ,'   ',Tools.NextSunrise(time0).hour ,'   ', time.day ,'   ', Tools.NextSunrise(time0).day)
                # Dark and window conditions are fulfilled. Append time
                if (Tools.IsDarkness(time,obsSite) is True) and (Tools.IsDarkness(time + WindowDuration,obsSite) is True):
                    NightDarkRuns.append(time)
                    time = time + WindowDuration
                    # print('first IF', str(time),'isdarks',Tools.IsDarkness(time))
                # Window condition is not fulfilled. Look if possible window is bigger that MinimalWindowDuration and if so, append time
                if (Tools.IsDarkness(time,obsSite) is True) and (Tools.IsDarkness(time + WindowDuration,obsSite) is False):
                    # print('Second  IF', str(time),'isdarks',Tools.IsDarkness(time))
                    # REgardless of being the sunrise or the moonrise the reason of darkness(the end of the window) == False , check if minimum requirement for a wondow is fulfilled
                    if ((MinimalWindowDuration <= (Tools.NextMoonrise(time,obsSite) - time) <= WindowDuration) or MinimalWindowDuration <= (Tools.NextSunrise(time,obsSite) - time) <= WindowDuration):
                        NightDarkRuns.append(time)
                    # print('But as darkness is not fulfilled anymore.Now its', str(time),'We jump to Nextmoonset', str(Tools.NextMoonset(time)))
                    # print('==============================================')
                    # print('MOON')
                    # print("str(Tools.PreviousMoonrise(time))  str(Tools.PreviousMoonset(time))  str(Tools.NextMoonrise(time))   str(Tools.NextMoonset(time))")
                    # print(str(Tools.PreviousMoonrise(time)),'   ',str(Tools.PreviousMoonset(time)),'  ',str(Tools.NextMoonrise(time)),'   ', str(Tools.NextMoonset(time)))
                    # print('==============================================')
                    # print('SUN')
                    # print("str(Tools.PreviousSunrise(time))  str(Tools.PreviousSunset(time))  str(Tools.NextSunrise(time))   str(Tools.NextSunset(time))")
                    # print(str(Tools.PreviousSunrise(time)), '   ', str(Tools.PreviousSunset(time)), '  ',str(Tools.NextSunrise(time)), '   ', str(Tools.NextSunset(time)))
                    # print('==============================================')
                    time = Tools.NextMoonset(time,obsSite)
                    time = Tools.TrustingDarknessMoon(time,
                                                      time0,obsSite)  # This case means: Night and darkness given by nextmoonset but when using IsDarkness() gives False due to decimals
                # Dark is false
                if (Tools.IsDarkness(time,obsSite) is False) and ((time.hour >= time0.hour and time.day == time0.day) or (
                        time.hour <= Tools.NextSunrise(time0,obsSite).hour and time.day == Tools.NextSunrise(time0,obsSite).day)):
                    time = Tools.NextMoonset(time,obsSite)
                    time = Tools.TrustingDarknessMoon(time, time0,obsSite)
            # Night is over, break
            else:
                #print('NIGHT IS OVER/BREAK')
                # print('time.hour >= time0.hour and time.day == time0.day) or (time.hour <= Tools.NextSunrise(time0).hour and time.day == Tools.NextSunrise(time0).day')
                # print(time.hour ,'   ', time0.hour,'    ', time.day ,'   ',time0.day ,'   ',time.hour ,'   ',Tools.NextSunrise(time0).hour ,'   ', time.day ,'   ', Tools.NextSunrise(time0).day)
                break
    return NightDarkRuns


def NightDarkObservationwithGreyTime(time, obsSite, MaxNights,duration,minduration):
    '''
    Function that searches for an array of observation times that fulfilled darkness condition and window

    '''
    WindowDuration = datetime.timedelta(minutes=duration)
    MinimalWindowDuration = datetime.timedelta(minutes=minduration)
    AuxMax = 100
    NightDarkRuns = []
    isFirstNight = True
    # Loop for the nights of observation
    for i in range(0, MaxNights):
        # In case the alert arrives during the day, time is set to sunset
        if (Tools.NextSunset(time,obsSite).hour >= time.hour >= Tools.PreviousSunrise(
                time,obsSite).hour and time.day == Tools.NextSunset(time,obsSite).day):
            time = Tools.NextSunset(time,obsSite)
            #print('POSTTIME',time,'isdarks', Tools.IsDarkness(time,obsSite))
            time = Tools.TrustingGreynessSun(time,obsSite)
            #print('POSTPOSTTIME', time, 'isdarks', Tools.IsDarkness(time,obsSite))
            isFirstNight = False
        # in case the alert arrives during the night, the last observation night will be the fourth so MaxNights=4
        else:
            if (isFirstNight):
                #MaxNights = MaxNights+1
                MaxNights = MaxNights
                isFirstNight = False
        # Using Time 0 in order to define what would be a night through that goes from Time 0 to NextSunrise(Time 0)
        time0 = time
        for j in range(0, AuxMax):
            # Night condition fulfilled
            #print('Times', str(time))
            if ((time.hour >= time0.hour and time.day == time0.day) or (
                    time.hour <= Tools.NextSunrise(time0,obsSite).hour and time.day == Tools.NextSunrise(time0,obsSite).day)):
                # print('time.hour >= time0.hour and time.day == time0.day) or (time.hour <= Tools.NextSunrise(time0).hour and time.day == Tools.NextSunrise(time0).day')
                # print(time.hour ,'   ', time0.hour,'    ', time.day ,'   ',time0.day ,'   ',time.hour ,'   ',Tools.NextSunrise(time0).hour ,'   ', time.day ,'   ', Tools.NextSunrise(time0).day)
                # Dark and window conditions are fulfilled. Append time
                if (Tools.IsGreyness(time,obsSite) is True) and (Tools.IsGreyness(time + WindowDuration,obsSite) is True):
                    NightDarkRuns.append(time)
                    time = time + WindowDuration
                    #print('first IF', str(time),'isdarks',Tools.IsDarkness(time,obsSite),'is grey',Tools.IsGreyness(time,obsSite))
                # Window condition is not fulfilled. Look if possible window is bigger that MinimalWindowDuration and if so, append time
                if (Tools.IsGreyness(time,obsSite) is True) and (Tools.IsGreyness(time + WindowDuration,obsSite) is False):
                    #print('second IF', str(time),'isdarks',Tools.IsDarkness(time,obsSite),'is grey',Tools.IsGreyness(time,obsSite))
                    #print(Tools.IsGreyness(time + MinimalWindowDuration,obsSite))
                    windowbool, endwind = Tools.TrustingGreynessMoon(time,time0,obsSite)
                    #print(windowbool)
                    if windowbool:
                        NightDarkRuns.append(time)
                    time=endwind
                    #REgardless of being the sunrise or the moonrise the reason of darkness(the end of the window) == False , check if minimum requirement for a wondow is fulfilled
                    #if ((MinimalWindowDuration <= (Tools.NextMoonrise(time,obsSite) - time) <= WindowDuration) or MinimalWindowDuration <= (
                    #         Tools.NextSunrise(time,obsSite) - time) <= WindowDuration):
                    #    NightDarkRuns.append(time)
                    # This case means: Night and darkness given by nextmoonset but when using IsDarkness() gives False due to decimals
                    # print('But as darkness is not fulfilled anymore.Now its', str(time),'We jump to Nextmoonset', str(Tools.NextMoonset(time)))
                    # print('==============================================')
                    # print('MOON')
                    # print("str(Tools.PreviousMoonrise(time))  str(Tools.PreviousMoonset(time))  str(Tools.NextMoonrise(time))   str(Tools.NextMoonset(time))")
                    # print(str(Tools.PreviousMoonrise(time)),'   ',str(Tools.PreviousMoonset(time)),'  ',str(Tools.NextMoonrise(time)),'   ', str(Tools.NextMoonset(time)))
                    # print('==============================================')
                    # print('SUN')
                    # print("str(Tools.PreviousSunrise(time))  str(Tools.PreviousSunset(time))  str(Tools.NextSunrise(time))   str(Tools.NextSunset(time))")
                    # print(str(Tools.PreviousSunrise(time)), '   ', str(Tools.PreviousSunset(time)), '  ',str(Tools.NextSunrise(time)), '   ', str(Tools.NextSunset(time)))
                    # print('==============================================')

                    #time = Tools.NextMoonset(time,obsSite)
                # Dark is false
                if (Tools.IsGreyness(time,obsSite) is False) and ((time.hour >= time0.hour and time.day == time0.day) or (time.hour <= Tools.NextSunrise(time0,obsSite).hour and time.day == Tools.NextSunrise(time0,obsSite).day)):
                    #print('IsGreyness is',Tools.IsGreyness(time,obsSite))
                    ####time = Tools.NextMoonset(time,obsSite)
                    ####time = Tools.TrustingDarknessMoon(time, time0,obsSite)
                    break
            else:
                #print('NIGHT IS OVER/BREAK')
                # print('time.hour >= time0.hour and time.day == time0.day) or (time.hour <= Tools.NextSunrise(time0).hour and time.day == Tools.NextSunrise(time0).day')
                # print(time.hour ,'   ', time0.hour,'    ', time.day ,'   ',time0.day ,'   ',time.hour ,'   ',Tools.NextSunrise(time0).hour ,'   ', time.day ,'   ', Tools.NextSunrise(time0).day)
                break
    return NightDarkRuns


def ZenithAngleCut(prob, nside, time, MinProbCut,max_zenith,observatory,usegreytime):
    '''
    Mask in the pixels with zenith angle larger than 45
    '''
    #observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)
    frame = co.AltAz(obstime=time, location=observatory)
    pprob = prob

    mzenith = hp.ma(pprob)
    maskzenith = np.zeros(hp.nside2npix(nside), dtype=np.bool)

    pixel_theta, pixel_phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    ra = np.rad2deg(pixel_phi)
    dec = np.rad2deg(0.5 * np.pi - pixel_theta)
    targetCoord_map = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    altaz_map = targetCoord_map.transform_to(frame)
    maskzenith[altaz_map.alt.value < 90-max_zenith] = 1
    mzenith.mask = maskzenith
    #hp.mollview(mzenith)
    #plt.show()
    #plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_%g.png")

    yprob = ma.masked_array(pprob, mzenith.mask)
    #hp.mollview(yprob)
    #plt.savefig("/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/Zenithmask_prob_%g.png")

    #print('Integrated probability of the masked map', np.sum(yprob))

    if np.sum(yprob) < MinProbCut:
        ObsBool = False
    else:
        ObsBool = True

    if usegreytime and ObsBool:
        # Get Alt/Az of the Moon
        moonaltazs = get_moon(Time(time)).transform_to(AltAz(obstime=Time(time),
                                                        location=observatory))
        separations = altaz_map.separation(moonaltazs)
        mask_moonDistance = np.zeros(hp.nside2npix(nside), dtype=np.bool)
        mask_moonDistance[separations < MoonSourceSeparation * u.deg]=1
        mzenith = hp.ma(pprob)
        mzenith.mask = mask_moonDistance
        yprob = ma.masked_array(pprob, mzenith.mask)
        #hp.mollview(pprob)
        #hp.mollview(yprob)
        #plt.show()
        if np.sum(yprob) < MinProbCut:
            ObsBool = False
        else:
            ObsBool = True
        #print('Integrated probability of the masked map', np.sum(yprob))
        #hp.mollview(mzenith)
        #plt.show()
        # Get the mask that does a radius around, of 30 degs
        # Plot to check
        # Return a bool if there is any observable region

    return ObsBool, yprob

def ComputeProbability2D(prob,highres, radecs,ReducedNside,HRnside,MinProbCut, time,observatory,max_zenith, FOV, tname, ipixlist,ipixlistHR, counter,dirName, usegreytime,plot):
    '''
    Compute probability in 2D by taking the highest value pixel
    '''
    radius = FOV
    #P_GW = prob[ipix_disc].sum()
    frame = co.AltAz(obstime=time, location=observatory)
    thisaltaz = radecs.transform_to(frame)
    #pix_alt1 = thisaltaz.alt.value

    if usegreytime:
        moonaltazs = get_moon(Time(time)).transform_to(AltAz(obstime=Time(time),location=observatory))
        #Zenith and Moon angular distance mask
        pix_ra = radecs.ra.value[(thisaltaz.alt.value > 90-max_zenith)&(thisaltaz.separation(moonaltazs)>30* u.deg)]
        pix_dec = radecs.dec.value[(thisaltaz.alt.value > 90 - max_zenith)&(thisaltaz.separation(moonaltazs)>30* u.deg)]

    else:
        # Zenith angle mask
        pix_ra = radecs.ra.value[(thisaltaz.alt.value > 90-max_zenith)]
        pix_dec = radecs.dec.value[thisaltaz.alt.value > 90 - max_zenith]
        #pix_alt = pix_alt1[thisaltaz.alt.value > 90 - max_zenith]


    phipix = np.deg2rad(pix_ra)
    thetapix = 0.5 * np.pi - np.deg2rad(pix_dec)

    ipix = hp.ang2pix(ReducedNside,thetapix, phipix)

    dp_Pix_Fov = np.empty(len(pix_ra), dtype=object)

    cat_pix = Table([ipix,pix_ra, pix_dec,dp_Pix_Fov], names=('PIX','PIXRA', 'PIXDEC', 'PIXFOVPROB'))

    dp_dV_FOV = []

    xyzpix = hp.ang2vec(thetapix, phipix)

    for i in range(0, len(cat_pix)):
        ipix_discfull = hp.query_disc(HRnside, xyzpix[i], np.deg2rad(radius))
        maskComputeProb = [np.isin(ipix_discfull, ipixlistHR, invert=True)]
        dp_dV_FOV.append(highres[ipix_discfull[maskComputeProb]].sum())

    cat_pix['PIXFOVPROB'] = dp_dV_FOV

    # Mask already observed pixels
    #print('ipixlist',ipixlist)

    mask=[np.isin(cat_pix['PIX'], ipixlist,invert=True)]

    if all(np.isin(cat_pix['PIX'], ipixlist,invert=False)):
        maskcat_pix = cat_pix
    else:
        maskcat_pix=cat_pix[mask]

    # Sort table
    sortcat = maskcat_pix[np.flipud(np.argsort(maskcat_pix['PIXFOVPROB']))]
    # Chose highest

    targetCoord = co.SkyCoord(sortcat['PIXRA'][:1],sortcat['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))
    #print('targetCoord',targetCoord)

    P_GW = sortcat['PIXFOVPROB'][:1]

    # Include to the list of pixels already observed


    #ipix_discComputeProb = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
    #maskComputeProb=[np.isin(ipix_discComputeProb, ipixlist,invert=True)]
    #print('maskedP_GW',highres[ipix_discComputeProb[maskComputeProb]].sum())
    if(P_GW >= MinProbCut):
        phip = float(np.deg2rad(targetCoord.ra.deg))
        thetap = float(0.5 * np.pi - np.deg2rad(targetCoord.dec.deg))
        xyz = hp.ang2vec(thetap, phip)
        
        ipixlistHR.extend(hp.query_disc(HRnside, xyz, np.deg2rad(radius)))
        ipix_disc = hp.query_disc(ReducedNside, xyz, np.deg2rad(radius))
        ipixlist.extend(ipix_disc)

    ######################################

    # PLOT THE RESULTS
    if plot:
        path = dirName + '/EvolutionPlot'
        if not os.path.exists(path):
            os.mkdir(path, 493)
        #nside = 1024

        #hp.mollview(highres,title="With FoV circle")

        hp.gnomview(prob, xsize=500, ysize=500, rot=[targetCoord.ra.deg, targetCoord.dec.deg], reso=8.0)
        hp.graticule()
        #print('This skymap has nside equals to',hp.npix2nside(len(highres)))
        plt.savefig('%s/Pointing-prob_%g.png' % (path,counter))

        ipix_discplot = hp.query_disc(HRnside, xyz, np.deg2rad(radius))
        tt, pp = hp.pix2ang(HRnside, ipix_discplot)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)
        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
        #hp.visufunc.projplot(skycoord.ra, skycoord.dec, 'y.', lonlat=True, coord="C")
        #plt.show()
        #observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)
        
        hp.visufunc.projplot(sortcat['PIXRA'][:1],sortcat['PIXDEC'][:1], 'r.', lonlat=True, coord="C")
        MaxCoord = SkyCoord(sortcat['PIXRA'][:1],sortcat['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))
        separations = skycoord.separation(MaxCoord)
        tempmask = separations < (radius + 0.1 * radius) * u.deg
        tempmask2 = separations > (radius - 0.1  * radius) * u.deg
        hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'g.', lonlat=True, coord="C",linewidth=0.1)
        altcoord = np.empty(100)
        azcoord = np.random.rand(100) * 360
        #for i in range(0,1):
        #    altcoord.fill(90-(max_zenith-5*i))
        #    RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,location=observatory)
        #    RandomCoord_radec = RandomCoord.transform_to('fk5')
        #    hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        #plt.show()
        #plt.savefig('%s/Pointing-zencut_%g.png' % (path,counter))
        
    return P_GW, targetCoord, ipixlist,ipixlistHR

def SubstractPointings2D(tpointingFile,prob,nside,FOV,pixlist):
    radius=FOV
    print("Loading pointings from " + tpointingFile)
    rap, decP = np.genfromtxt(tpointingFile, usecols=(2, 3), dtype="str", skip_header=1,
                                             delimiter=' ',
                                             unpack=True)  # ra, dec in degrees
    coordinates=TransformRADec(rap, decP)
    P_GW = []
    for i in range(0,len(rap)):
        t = 0.5 * np.pi - coordinates[i].dec.rad
        p = coordinates[i].ra.rad
        xyz = hp.ang2vec(t, p)
        ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
        effectiveipix_disc = []
        #print('ipixdisc',ipix_disc)
        #print('pixlist',pixlist)
        for j in range(0, len(ipix_disc)):
            if not (ipix_disc[j] in pixlist):
                effectiveipix_disc.append(ipix_disc[j])
            pixlist.append(ipix_disc[j])
        P_GW.append(prob[effectiveipix_disc].sum())
        #print('effectiveipix_disc',effectiveipix_disc)
        #print('-----------')
        #print(prob[effectiveipix_disc].sum())
        #print(prob[ipix_disc].sum())
        #print(prob[pixlist].sum())
        #print('-----------')
    return pixlist,np.sum(P_GW)

######################################################
def TransformRADec(vra,vdec):
    if('h' in vra[0]):
        ra=[]
        dec=[]
        for i in range(0,len(vra)):
            coord=SkyCoord(vra[i].split('"')[1],vdec[i].split('"')[0],frame='fk5')
            print(coord)
            ra.append(coord.ra.deg)
            dec.append(coord.dec.deg)
    else:
        print(vra,vdec)
        ra = vra.astype(np.float)
        dec = vdec.astype(np.float)
        #np.float(vra)
        #dec = np.float(vdec)
    print(ra,dec)
    coordinates = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    return coordinates

## Extra functions from BestCandidateon PGal ==> 3D

######################################################

def LoadGalaxies(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    ra, dec, dist, bmag = np.genfromtxt(tgalFile, usecols=(0, 1, 2, 3), skip_header=1, unpack=True)  # ra, dec in degrees

    tcat = Table([ra, dec, dist, bmag], names=('RAJ2000', 'DEJ2000', 'Dist', 'Bmag'))
    return tcat

def LoadGalaxies_SteMgal(tgalFile):
    '''
    Load galaxy catalog as an Astropy Table
    '''

    print("Loading galaxy catalogue from " + tgalFile)

    ra, dec, dist, bmag, mgal = np.genfromtxt(tgalFile, usecols=(0, 1, 2, 3, 4), skip_header=1, unpack=True)  # ra, dec in degrees

    tcat = Table([ra, dec, dist, bmag, mgal], names=('RAJ2000', 'DEJ2000', 'Dist', 'Bmag', 'SteMgal'))
    return tcat

def CorrelateGalaxies_LVC(prob, distmu, distsigma, distnorm, cat, Info3D_available,MinimumProbCutForCatalogue):
    '''
    Correlates galaxies with GW 3D information following Going the Distance, then sort catalog by that value
    In case there is no 3 extra layers, the probability in 2D is asigned.
    '''
    ra=cat['RAJ2000']
    dec=cat['DEJ2000']
    dist=cat['Dist']

    # Translate RA,Dec of galaxies into theta,phi angles
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)


    # Get corresponding healpix pixel IDs
    npix = len(prob)
    nside = hp.npix2nside(npix)
    ipix = hp.ang2pix(nside, theta, phi)

    # Calculate probability in the space volumes

    pixarea = hp.nside2pixarea(nside)

    if(Info3D_available):
        dp_dV = prob[ipix] * distnorm[ipix] * norm(distmu[ipix], distsigma[ipix]).pdf(dist) / pixarea

    else:
        dp_dV = prob[ipix] / pixarea

    # Add dp_dV to catalogue
    cat['dp_dV'] = dp_dV
    
    # Select all values > 1% of the peak prob.

    total_dP_dV = dp_dV.sum()
    min_prob_cut = dp_dV > MinimumProbCutForCatalogue  * max(dp_dV)
    Gals = cat[min_prob_cut]
    # return array with list of Galaxies passing cuts, ordered by p-value

    tGals = Gals[np.flipud(np.argsort(Gals['dp_dV']))]
    #ascii.write(tGals, '/Users/hashkar/Desktop/GWfollowup/GW-Followup/tGals_noM.txt', names = ['RAJ2000','DEJ2000','Dist','Bmag','dp_dV'],overwrite=True)
    return tGals, total_dP_dV

def CorrelateGalaxies_LVC_SteMass(prob, distmu, distsigma, distmean, disterr, distnorm, cat, Info3D_available,MinimumProbCutForCatalogue):
    '''
    Correlates galaxies with GW 3D information following Going the Distance, then sort catalog by that value
    In case there is no 3 extra layers, the probability in 2D is asigned.
    '''
    beta = 1
    alpha = 0
    ra=cat['RAJ2000']
    dec=cat['DEJ2000']
    dist=cat['Dist']

    # Translate RA,Dec of galaxies into theta,phi angles

    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)


    # Get corresponding healpix pixel IDs
    npix = len(prob)
    nside = hp.npix2nside(npix)
    ipix = hp.ang2pix(nside, theta, phi)

    #get the pixles that are confined in the 90% area
    pixtab = Get90RegionPixGal(prob, 0.9, nside)


    #create new catalog with galaxies inside 90% region
    cat['pixtab'] = ipix
    pixtablist = np.in1d(ipix,pixtab)
    Gals = cat[pixtablist]

    #filter out galaxies outside the distance uncertaintity
    dist = Gals['Dist']
    distok = (dist < (distmean+2*disterr)) & (dist > (distmean-2*disterr))
    Gals = Gals[distok]

    #compute a new ipix
    ra=Gals['RAJ2000']
    dec=Gals['DEJ2000'] 
    dist = Gals['Dist']

    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)

    NewIpix  = hp.ang2pix(nside, theta, phi)


    # Calculate probability in the space volumes
    pixarea = hp.nside2pixarea(nside)

    if(Info3D_available):
        dp_dV_pos = prob[NewIpix] * distnorm[NewIpix] * norm(distmu[NewIpix], distsigma[NewIpix]).pdf(dist) / pixarea


    else:
        dp_dV_pos = prob[NewIpix] / pixarea

    # Add dp_dV to new catalogue
    Gals['dp_dV'] = dp_dV_pos

    # Select all values > 1% of the peak prob.

    total_dP_dV = dp_dV_pos.sum()
    min_prob_cut = dp_dV_pos > MinimumProbCutForCatalogue  * max(dp_dV_pos)
    Gals = Gals[min_prob_cut]

    if(Info3D_available):

        Mgal1 =  Gals['SteMgal']
        Pgal_pos = Gals['dp_dV']

        Mgal1 = np.nan_to_num(Mgal1)
        Mgal = 10**(Mgal1)
        Pgal_pos = np.nan_to_num(Pgal_pos)


        Gmass = Mgal/(np.sum(Mgal))
        alpha = (Pgal_pos).sum()/(Pgal_pos*Gmass).sum()
        dp_dV = (Pgal_pos)+(Pgal_pos*(alpha*beta*Gmass))

    Gals['dp_dV'] = dp_dV  

    total_dP_dV = dp_dV.sum()
    print(total_dP_dV)

    tGals = Gals[np.flipud(np.argsort(Gals['dp_dV']))]
    #ascii.write(tGals, '/Users/hashkar/Desktop/GWfollowup/GW-Followup/tGals.txt', names = ['RAJ2000','DEJ2000','Dist','Bmag','SteMgal', 'index' ,'dp_dV'],overwrite=True)

    return tGals, total_dP_dV

def VisibleAtTime(test_time, galaxies, maxz,observatory):
    '''Determine if prompt or afterglow follow-up is possible by knowing if there are galaxies with non-negligible probability of hosting the NSM in the FoV
     1) check if any galaxy is visible, if not --> AFTERGLOW

    2) loop over zenith angle and select subsets of galaxies

    3) stop if maximum p-value of this subset is smaller than 75% of the previous subset

    4) else: stricter cut on zenith and repeat

    5) take galaxy with highest p-value fulfilling both criteria as target

    RETURNS:
    --------
    bool `is_vis` : is visible now?
    np.ndarray `alt_az` : alt_az location of galaxies
    '''

    #print()
    #print("Check visibility at time {0}".format(test_time))

    # observatory time and location to look up visibility of objects

    #observatory = co.EarthLocation(lat=-23.271333 * u.deg,lon=16.5 * u.deg, height=1800 * u.m)

    frame = co.AltAz(obstime=test_time, location=observatory)
    #print('galaxies',galaxies)
    #print('galaxies',len(galaxies['RAJ2000']))
    radecs = co.SkyCoord(galaxies['RAJ2000'], galaxies['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    if (len(radecs) > 0):
        thisaltaz = radecs.transform_to(frame)

        # add altitude to topGals array
        # already sorted by descending probability value

        galaxies['Alt'] = thisaltaz.alt.value

        # check if any galaxy is visible at the moment
        thismask = thisaltaz.alt.value > 90 - maxz

        nGals = len(galaxies[thismask])

        #print('nGals',nGals)

        if (nGals == 0):
            #print("No galaxies visible within {0} deg zenith angle --> AFTERGLOW".format(maxz))

            return False, thisaltaz, galaxies
        else:
            #print("{0} galaxies are visible within {1} deg zenith angle ""--> Test for prompt follow up".format(nGals, maxz))

            return True, thisaltaz, galaxies
    else:
        thisaltaz = []
        return False, thisaltaz, galaxies

def FulfillsRequirement(theseGals, maxz,FOV,FulFillReq_Percentage,UsePix):
    '''
    Apply filter criteria to visible galaxy sample and compares them to get the best option of zenith angle

    '''

    #print("Check if galaxies with minimum p-value are among candidates...")

    # Initialise maximum p-value in map to 1

    maxp = 1
    mask = 0
    thisminz = 0

    alt = theseGals['Alt']

    for minz_aux in range(maxz, 5, -5):

        tmpmask = alt > 90 - (minz_aux)
        #print('TheseGals', theseGals[tmpmask])
        tmpGals = theseGals.copy()
        #print('len(tmpGals[tmpmask])', len(tmpGals[tmpmask]), 'without mask', len(tmpGals))
        # cut on zenith angle and select most probable galaxy

        if (len(tmpGals[tmpmask]) > 0):

            cur_maxp = tmpGals[tmpmask]['dp_dV'].max() / theseGals['dp_dV'].max()

            #print("{0} galaxies visible with zen < {1} deg - maximum p-value {2:0.3f}"

            #     .format(len(tmpGals[tmpmask]), minz_aux, cur_maxp))

           #print("Maximum probability of {0:0.3f} of global maximum of {1:3f}"

            #      .format(cur_maxp, theseGals['dp_dV'].max()))

            # define final mask

            if (maxz == minz_aux):
                maxp = cur_maxp
                mask = tmpmask
                thisminz = minz_aux

            if (cur_maxp > FulFillReq_Percentage * maxp):
                mask = tmpmask
                thisminz = minz_aux
            else:
                thisminz = minz_aux + 5
                break
    if UsePix:
        mask=alt >90-(thisminz+FOV)
    return mask, thisminz

def FulfillsRequirementGreyObservations(Ktime,theseGals,observatory):

    targetCoord = co.SkyCoord(theseGals['RAJ2000'], theseGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    frame = co.AltAz(obstime=Ktime, location=observatory)
    moonaltazs = get_moon(Time(Ktime)).transform_to(frame)

    altaz_map = targetCoord.transform_to(frame)
    separations = altaz_map.separation(moonaltazs)

    #Mask
    greymask=separations>MoonSourceSeparation*u.deg
    return greymask

def FulfillsRequirement_MinProb(thisGals_aux, maxz):
    ''' Same as FulfillsRequirements but at the end a supplementary cut is performed.
    This algorithm comes from the first version of the code but in the newer versions
    the algorithms separates this two options.

    '''

    print("Check if galaxies with minimum p-value are among candidates...")

    # Initialise maximum p-value in map to 1

    maxp = 1
    mask = 0

    for minz_aux in range(maxz, 5, -5):

        tmpmask = altaz.alt.value > 90 - minz_aux
        tmpGals = thisGals_aux.copy()

        # cut on zenith angle and select most probable galaxy

        if (len(tmpGals[tmpmask]) > 0):

            cur_maxp = tmpGals[tmpmask]['dp_dV'].max() / tGals['dp_dV'].max()

            print("{0} galaxies visible with zen < {1} deg - maximum p-value {2:0.3f}".format(len(tmpGals[tmpmask]),
                                                                                              minz_aux, cur_maxp))
            print("Maximum probability of {0:0.3f} of global maximum of {1:3f}".format(cur_maxp, tGals['dp_dV'].max()))

            if (maxz == minz_aux):
                maxp = cur_maxp
                mask = tmpmask
                minz = minz_aux

            if (cur_maxp > 0.75 * maxp):
                mask = tmpmask
                minz = minz_aux
            else:
                minz = minz_aux + 5
                break

    if (maxp < 0.02 * thisGals_aux['dp_dV'].max()):

        print("Probability too low, postpone observations --> AFTERGLOW")
        return False, mask, minz

    else:
        print('This minz= ', minz)
        return True, mask, minz


def Afterglow():
    print('Afterglow!')

def ComputeProbBCFOVSimple(prob,time, visiGals, allGals, tsum_dP_dV, nside, thisminz,max_zenith, FOV, tname,
                           tsavedcircle,dirName, doplot):
    '''Computes probability pgal and pgw in FoV and draws everything

    bool doplot when  = True is used to plot the maps

    RETURNS:

    --------

	P_Gal: Probability of galaxies within H.E.S.S. FoV in the LIGO signal region
	P_GW: Total probability within H.E.S.S. FoV of the Ligo signal.
	noncircleGal: Table of galaxies that are outside the circle(s) and inside the LIGO signal region


    '''

    targetCoord = co.SkyCoord(visiGals['RAJ2000'][:1], visiGals['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(visiGals['RAJ2000'], visiGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord3 = co.SkyCoord(allGals['RAJ2000'], allGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    dp_dVfinal = visiGals['dp_dV']
    # dp_dV = tGals['dp_dV']

    # Array of indices of pixels inside circle of HESS-I FoV

    radius = FOV

    t = 0.5 * np.pi - targetCoord[0].dec.rad

    p = targetCoord[0].ra.rad

    # print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg',t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg)

    xyz = hp.ang2vec(t, p)

    # print(xyz)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    # print('ipix_disc',ipix_disc)
    P_GW = prob[ipix_disc].sum()

    P_Gal = dp_dVfinal[targetCoord2.separation(targetCoord).deg < radius].sum() / tsum_dP_dV

    # print("Total probability within H.E.S.S. FoV: {0}".format(P_GW))

    # print("Probability of galaxies within H.E.S.S. FoV in the LIGO signal region:{0}".format(P_Gal))

    # all galaxies inside the current observation circle

    circleGal = visiGals[targetCoord2.separation(targetCoord).deg < radius]
    # print('Galaxies within the FoV: ', len(circleGal['RAJ2000']))

    # all galaxies outside the current observation circle, no visibility selection

    noncircleGal = allGals[targetCoord3.separation(targetCoord).deg > radius]

    if (doplot):
        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
        observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

        frame = co.AltAz(obstime=time, location=observatory)
        altaz_all = skycoord.transform_to(frame)
        tmask = altaz_all.alt.value > 90 - thisminz

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        path = dirName+ '/EvolutionPlot'
        if not os.path.exists(path):
            os.mkdir(path, 493)

        hp.mollview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (
            time.day, time.month, time.year, time.hour, time.minute))

        # hp.gnomview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (time.day, time.month, time.year, time.hour, time.minute),xsize = 4000,ysize=6000,rot = [90,-50],reso=0.8)
        hp.graticule()
        # plt.show()
        # plt.savefig("Figures/ExampleGW_%g.png" % (j))

        # draw all galaxies within zenith-angle cut
        # hp.visufunc.projscatter(finalGals['RAJ2000'], finalGals['DEJ2000'], lonlat=True, marker='*', color='g')
        # plt.savefig("Figures/ExampleGW_Galaxies_%g.png" % (j))

        # If I want to plot all gals, plot also the ones that are out of the circle
        # hp.visufunc.projscatter(noncircleGal['RAJ2000'], noncircleGal['DEJ2000'], lonlat=True, marker='*', color='g')

        # draw observation position, which is equivalent to galaxy with highest
        # probability
        hp.visufunc.projscatter(visiGals['RAJ2000'][:1], visiGals['DEJ2000'][:1], lonlat=True, marker='.', color='r',
                                linewidth=0.1)

        # draw circle of HESS-I FoV around best fit position

        # hp.visufunc.projscatter(allGals['RAJ2000'], allGals['DEJ2000'], lonlat=True, marker='.', color='g',linewidth=0.1)
        # plt.show()
        hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'k.', lonlat=True,
                             coord="C", linewidth=0.1)

        # hp.visufunc.projplot(tsavedcircle.ra, tsavedcircle.dec, 'r.', lonlat=True, coord="C",linewidth=0.1)

        # Draw H.E.S.S. visibility

        # altcoord= [np.random.randint(-90,90-thisminz) for _ in range(4000)]

        altcoord = np.empty(4000)

        altcoord.fill(90 - max_zenith)

        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                               location=observatory)

        RandomCoord_radec = RandomCoord.transform_to('fk5')

        hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C", linewidth=0.1)
        # plt.show()
        # Draw MinZ area

        print('Min Zenith= ', thisminz)

        altcoordmin = np.empty(4000)

        # altcoordmin.fill(90 - thisminz)
        azcoordmin = np.random.rand(4000) * 360

        RandomCoordmin = SkyCoord(azcoordmin, altcoordmin, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                                  location=observatory)

        RandomCoordmin_radec = RandomCoordmin.transform_to('fk5')

        # hp.visufunc.projplot(RandomCoordmin_radec.ra, RandomCoordmin_radec.dec, 'y.', lonlat=True, coord="C", marker='.', markersize = 8 )

        # plt.show()
        plt.savefig("%s/ExamplePointing_%g.png" % (path, len(ObservationTimearray)))

    return P_Gal, P_GW, talreadysumipixarray2

def ComputeProbBCFOV(prob,time, finalGals, visiGals, allGals, tsum_dP_dV, talreadysumipixarray, nside, thisminz,max_zenith, FOV, tname,
                     tsavedcircle, dirName,doplot):
    '''Computes probability Pgal and Pgw in FoV but it takes into account a list of pixels to avoid recounting already observed zones.
    Returns saved circle too (is it really needed? )
    bool doplot when  = True is used to plot the maps

    RETURNS:

    --------

	P_Gal: Probability of galaxies within H.E.S.S. FoV in the LIGO signal region
	P_GW: Total probability within H.E.S.S. FoV of the Ligo signal.
	noncircleGal: Table of galaxies that are outside the circle(s) and inside the LIGO signal region


    '''

    targetCoord = co.SkyCoord(finalGals['RAJ2000'][:1], finalGals['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(visiGals['RAJ2000'], visiGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord3 = co.SkyCoord(allGals['RAJ2000'], allGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    dp_dVfinal = visiGals['dp_dV']
    # dp_dV = tGals['dp_dV']

    # Array of indices of pixels inside circle of HESS-I FoV

    radius = FOV

    t = 0.5 * np.pi - targetCoord[0].dec.rad

    p = targetCoord[0].ra.rad

    # print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg',t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg)

    xyz = hp.ang2vec(t, p)

    # print(xyz)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    # print('ipix_disc',ipix_disc)

    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])
        talreadysumipixarray.append(ipix_disc[j])

    # print('talreadysumipixarray', talreadysumipixarray)

    P_GW = prob[effectiveipix_disc].sum()

    P_Gal = dp_dVfinal[targetCoord2.separation(targetCoord).deg < radius].sum() / tsum_dP_dV

    # print("Total probability within H.E.S.S. FoV: {0}".format(P_GW))

    # print("Probability of galaxies within H.E.S.S. FoV in the LIGO signal region:{0}".format(P_Gal))

    # all galaxies inside the current observation circle

    circleGal = visiGals[targetCoord2.separation(targetCoord).deg < radius]
    # print('Galaxies within the FoV: ', len(circleGal['RAJ2000']))

    # all galaxies outside the current observation circle, no visibility selection

    noncircleGal = allGals[targetCoord3.separation(targetCoord).deg > radius]

    if (doplot):
        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
        observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

        frame = co.AltAz(obstime=time, location=observatory)
        altaz_all = skycoord.transform_to(frame)
        tmask = altaz_all.alt.value > 90 - thisminz

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        path = dirName + '/EvolutionPlot'
        if not os.path.exists(path):
            os.mkdir(path, 493)

        #hp.mollview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (time.day, time.month, time.year, time.hour, time.minute))

        hp.gnomview(prob,xsize = 4000,ysize=6000,rot = [90,-50],reso=0.8)
        hp.graticule()
        # plt.show()
        # plt.savefig("Figures/ExampleGW_%g.png" % (j))

        # draw all galaxies within zenith-angle cut
        # hp.visufunc.projscatter(finalGals['RAJ2000'], finalGals['DEJ2000'], lonlat=True, marker='*', color='g')
        # plt.savefig("Figures/ExampleGW_Galaxies_%g.png" % (j))

        # If I want to plot all gals, plot also the ones that are out of the circle
        # hp.visufunc.projscatter(noncircleGal['RAJ2000'], noncircleGal['DEJ2000'], lonlat=True, marker='*', color='g')

        # draw observation position, which is equivalent to galaxy with highest
        # probability
        ###### hp.visufunc.projscatter(finalGals['RAJ2000'][:1], finalGals['DEJ2000'][:1], lonlat=True, marker='.', color='r',linewidth=0.1)

        # draw circle of HESS-I FoV around best fit position

        ######hp.visufunc.projscatter(allGals['RAJ2000'], allGals['DEJ2000'], lonlat=True, marker='.', color='g',linewidth=0.1)
        # plt.show()
        hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'r.', lonlat=True,
                             coord="C", linewidth=0.1)

        # hp.visufunc.projplot(tsavedcircle.ra, tsavedcircle.dec, 'r.', lonlat=True, coord="C",linewidth=0.1)

        # Draw H.E.S.S. visibility

        # altcoord= [np.random.randint(-90,90-thisminz) for _ in range(4000)]

        altcoord = np.empty(4000)

        altcoord.fill(90 - max_zenith)

        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                               location=observatory)

        RandomCoord_radec = RandomCoord.transform_to('fk5')

        hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C", linewidth=0.1)
        # plt.show()
        # Draw MinZ area

        print('Min Zenith= ', thisminz)

        altcoordmin = np.empty(4000)

        altcoordmin.fill(90 - thisminz)
        azcoordmin = np.random.rand(4000) * 360

        RandomCoordmin = SkyCoord(azcoordmin, altcoordmin, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                                  location=observatory)

        RandomCoordmin_radec = RandomCoordmin.transform_to('fk5')

        hp.visufunc.projplot(RandomCoordmin_radec.ra, RandomCoordmin_radec.dec, 'y.', lonlat=True, coord="C", marker='.', markersize = 8 )

        # plt.show()
        tsavedcircle = skycoord[tempmask & tempmask2]
        plt.savefig("%s/ExamplePointing.png" % (path))

    return P_Gal, P_GW, noncircleGal, talreadysumipixarray, tsavedcircle

def SimpleGWprob(prob,finalGals,talreadysumipixarray,FOV,nside):
    '''Computes probability Pgw in FoV but it takes into account a list of pixels to avoid recounting already observed zones.
    bool doplot when  = True is used to plot the maps
    '''
    targetCoord = co.SkyCoord(finalGals['RAJ2000'][:1], finalGals['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))
    radius = FOV
    t = 0.5 * np.pi - targetCoord[0].dec.rad
    p = targetCoord[0].ra.rad
    xyz = hp.ang2vec(t, p)
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])

    # print('talreadysumipixarray', talreadysumipixarray)

    probability = prob[effectiveipix_disc].sum()
    return probability

def SubstractPointings(tpointingFile,galaxies,talreadysumipixarray,tsum_dP_dV,FOV,prob,nside):

    #targetCoord = co.SkyCoord(galaxies['RAJ2000'], galaxies['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    #Read PointingsFile

    print("Loading pointings from " + tpointingFile)
    rap, decP, = np.genfromtxt(tpointingFile, usecols=(3,4), dtype="str", skip_header=1,
                                             delimiter=' ',
                                             unpack=True)  # ra, dec in degrees

    coordinates=TransformRADec(rap, decP)
    ra=coordinates.ra.deg
    dec=coordinates.dec.deg

    PGW= []
    PGAL=[]

    updatedGalaxies = galaxies

    for i in range(0,len(ra)):
        #print('ra{1}', ra[i])
        updatedGalaxies,pgwcircle,pgalcircle,talreadysumipixarray =SubstractGalaxiesCircle(updatedGalaxies,ra[i], dec[i],talreadysumipixarray,tsum_dP_dV,FOV,prob,nside)
        PGW.append(pgwcircle)
        PGAL.append(pgalcircle)
    return ra,dec,updatedGalaxies, PGW, PGAL,talreadysumipixarray

def SubstractGalaxiesCircle(galaux, ra,dec ,talreadysumipixarray,tsum_dP_dV,FOV,prob,nside):
    radius = FOV
    coordinates = co.SkyCoord(ra, dec, frame='fk5', unit=(u.deg, u.deg))
    #print('coordinatesPointing',coordinates)
    targetCoord = co.SkyCoord(galaux['RAJ2000'], galaux['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))
    dp_dVfinal = galaux['dp_dV']
    #print('coordinatesPointing[0].dec.rad',coordinates.dec.rad)

    t = 0.5 * np.pi - coordinates.dec.rad
    p = coordinates.ra.rad
    xyz = hp.ang2vec(t, p)
    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))
    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])
        talreadysumipixarray.append(ipix_disc[j])

    P_GW = prob[effectiveipix_disc].sum()
    P_Gal = dp_dVfinal[targetCoord.separation(coordinates).deg < radius].sum() / tsum_dP_dV
    #Define galaxies that are outside the FoV

    newgalaxies = galaux[targetCoord.separation(coordinates).deg > radius]

    return newgalaxies, P_GW,P_Gal,talreadysumipixarray
#####################################################

## Extra functions from PGalinFoVOptimised ==> 3D ##

#####################################################

def PGalinFOV(prob,cat,galpix,FOV, totaldPdV,nside,UsePix):
    '''
        Computes probability Pgal in FoV
    '''
    if UsePix :
        targetCoord = co.SkyCoord(galpix['PIXRA'], galpix['PIXDEC'],frame='fk5', unit=(u.deg, u.deg))
    else:
        targetCoord = co.SkyCoord(galpix['RAJ2000'], galpix['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(cat['RAJ2000'], cat['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    #dp_dVfinal = finalGals['dp_dV']
    dp_dV = cat['dp_dV']

    # Array of indices of pixels inside circle of HESS-I FoV

    radius = FOV

    t = 0.5 * np.pi - targetCoord.dec.rad

    p = targetCoord.ra.rad

    #print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord.ra.deg, targetCoord.dec.deg)

    xyz = hp.ang2vec(t, p)

    #print(xyz)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))

    P_GW = prob[ipix_disc].sum()

    Pgal_inFoV = dp_dV[targetCoord2.separation(targetCoord).deg <= radius].sum() / totaldPdV

    return Pgal_inFoV


def ModifyCatalogue(prob,cat, FOV, totaldPdV,nside):
    '''
     Computes the integrated Pgal in FoV for a list of calues using Pgal in FoV and sorts the catalog
     using that quantity as a criteria
    '''
    #lengthSG=0.02*len(cat)
    lengthSG = 1000
    SelectedGals = cat[:lengthSG]
    dp_dV_FOV = []
    #print('len(cat[RAJ2000])', len(cat['RAJ2000']))
    #print('len(SelectedGals[RAJ2000])',len(SelectedGals['RAJ2000']))
    for l in range(0, len(cat['dp_dV'])):
        if(l<len(SelectedGals['dp_dV'])):
            dp_dV_FOV.append(PGalinFOV(prob,cat, SelectedGals[l], FOV, totaldPdV,nside,UsePix=False))
        else:
            dp_dV_FOV.append(0)

    cat['dp_dV_FOV'] = dp_dV_FOV

    tcat = cat[np.flipud(np.argsort(cat['dp_dV_FOV']))]

    return tcat


def ComputeProbPGALIntegrateFoV(prob,time, centerPoint,UsePix, visiGals, allGalsaftercuts, tsum_dP_dV, talreadysumipixarray, nside,
                                thisminz,max_zenith,FOV,counter, tname,dirName, doplot):
    '''
        Same as ComputeProbBCFOV but it does not return circle coordinates.
    '''

    if UsePix:
        targetCoord = co.SkyCoord(centerPoint['PIXRA'][:1], centerPoint['PIXDEC'][:1], frame='fk5', unit=(u.deg, u.deg))

    else:
        targetCoord = co.SkyCoord(centerPoint['RAJ2000'][:1], centerPoint['DEJ2000'][:1], frame='fk5', unit=(u.deg, u.deg))

    targetCoord2 = co.SkyCoord(visiGals['RAJ2000'], visiGals['DEJ2000'], frame='fk5', unit=(u.deg, u.deg))

    targetCoord3 = co.SkyCoord(allGalsaftercuts['RAJ2000'], allGalsaftercuts['DEJ2000'], frame='fk5',
                               unit=(u.deg, u.deg))

    dp_dVfinal = visiGals['dp_dV']

    # Array of indices of pixels inside circle of HESS-I FoV

    radius = FOV
    #print('FOV',FOV)
    t = 0.5 * np.pi - targetCoord[0].dec.rad

    p = targetCoord[0].ra.rad

    #print('t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg', t, p, targetCoord[0].ra.deg, targetCoord[0].dec.deg)

    xyz = hp.ang2vec(t, p)

    # translate pixel indices to coordinates

    ipix_disc = hp.query_disc(nside, xyz, np.deg2rad(radius))


    effectiveipix_disc = []

    for j in range(0, len(ipix_disc)):
        if not (ipix_disc[j] in talreadysumipixarray):
            effectiveipix_disc.append(ipix_disc[j])
        talreadysumipixarray.append(ipix_disc[j])

    P_GW = prob[effectiveipix_disc].sum()

    P_Gal = dp_dVfinal[targetCoord2.separation(targetCoord).deg < radius].sum() / tsum_dP_dV
    #TotalGal_FOV = len(dp_dVfinal[targetCoord2.separation(targetCoord).deg < radius])

    #print("Total probability within H.E.S.S. FoV: {0}".format(P_GW))

    #print("Probability of galaxies within H.E.S.S. FoV in the LIGO signal region:{0}".format(P_Gal))

    # all galaxies inside the current observation circle

    circleGal = visiGals[targetCoord2.separation(targetCoord).deg < radius]
    #print('Galaxies within the FoV: ', len(circleGal['RAJ2000']))

    # all galaxies outside the current observation circle, no visibility selection

    noncircleGal = allGalsaftercuts[targetCoord3.separation(targetCoord).deg > radius]

    if (doplot):

        path = dirName + '/EvolutionPlot'
        if not os.path.exists(path):
            os.mkdir(path, 493)

        tt, pp = hp.pix2ang(nside, ipix_disc)
        ra2 = np.rad2deg(pp)
        dec2 = np.rad2deg(0.5 * np.pi - tt)

        skycoord = co.SkyCoord(ra2, dec2, frame='fk5', unit=(u.deg, u.deg))
        observatory = co.EarthLocation(lat=-23.271333 * u.deg, lon=16.5 * u.deg, height=1800 * u.m)

        frame = co.AltAz(obstime=time, location=observatory)
        altaz_all = skycoord.transform_to(frame)
        tmask = altaz_all.alt.value > 90 - thisminz

        separations = skycoord.separation(targetCoord)
        tempmask = separations < (radius + 0.05 * radius) * u.deg
        tempmask2 = separations > (radius - 0.05 * radius) * u.deg

        #path = os.path.dirname(os.path.realpath(__file__)) + tname
        #if not os.path.exists(path):
        #    os.mkdir(path, 493)

        #hp.mollview(prob, title="GW prob map (Ecliptic)   %s/%s/%s %s:%s UTC" % (time.day, time.month, time.year, time.hour, time.minute), xsize=2000)
        hp.gnomview(prob, xsize=500, ysize=500, rot=[targetCoord.ra.deg, targetCoord.dec.deg], reso=5.0)

        hp.graticule()
        #plt.savefig("%s/ExampleGW_%g.png" % (tname,j))

        # draw all galaxies within zenith-angle cut
        ##### hp.visufunc.projscatter(allGalsaftercuts['RAJ2000'], allGalsaftercuts['DEJ2000'], lonlat=True, marker='.',color='g', linewidth=0.1)
        #plt.savefig("%s/ExampleGW_Galaxies_%g.png" % (tname,j))

        # If I want to plot all gals, plot also the ones that are out of the circle
        # hp.visufunc.projscatter(noncircleGal['RAJ2000'], noncircleGal['DEJ2000'], lonlat=True, marker='*', color='g')

        # draw observation position, which is equivalent to galaxy with highest
        # probability
        # hp.visufunc.projscatter(finalGals['RAJ2000'][:1], finalGals['DEJ2000'][:1], lonlat=True, marker='.', color='r',linewidth=0.1)

        # draw circle of HESS-I FoV around best fit position

        hp.visufunc.projplot(skycoord[tempmask & tempmask2].ra, skycoord[tempmask & tempmask2].dec, 'r.', lonlat=True,
                             coord="C")

        # Draw H.E.S.S. visibility

        # altcoord= [np.random.randint(-90,90-thisminz) for _ in range(4000)]

        altcoord = np.empty(4000)

        altcoord.fill(90 - max_zenith)

        azcoord = np.random.rand(4000) * 360

        RandomCoord = SkyCoord(azcoord, altcoord, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                               location=observatory)

        RandomCoord_radec = RandomCoord.transform_to('fk5')

        hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        #MOON


        #hp.visufunc.projplot(RandomCoord_radec.ra, RandomCoord_radec.dec, 'b.', lonlat=True, coord="C")
        # Draw MinZ area

        print('Min Zenith= ', thisminz)

        altcoordmin = np.empty(4000)

        altcoordmin.fill(90 - thisminz)

        azcoordmin = np.random.rand(4000) * 360

        RandomCoordmin = SkyCoord(azcoordmin, altcoordmin, frame='altaz', unit=(u.deg, u.deg), obstime=time,
                                  location=observatory)

        RandomCoordmin_radec = RandomCoordmin.transform_to('fk5')

        # hp.visufunc.projplot(RandomCoordmin_radec.ra, RandomCoordmin_radec.dec, 'y.', lonlat=True, coord="C")

        # plt.show()
        plt.savefig("%s/ExamplePointing%g.png" % (path,counter))

    return P_Gal, P_GW, noncircleGal, talreadysumipixarray

def MOC_confidence_region(infile, percentage, short_name=' ', save2File=False):

    # reading skymap
    hpx = hp.read_map(infile, verbose=False)
    npix = len(hpx)
    nside = hp.npix2nside(npix)

    sort = sorted(hpx, reverse=True)
    cumsum = np.cumsum(sort)
    index, value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))

    # finding ipix indices confined in a given percentage
    index_hpx = range(0, len(hpx))
    hpx_index = np.c_[hpx, index_hpx]

    sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
    value_contour = sort_2array[0:index]

    j = 1
    table_ipix_contour = []

    for i in range(0, len(value_contour)):
        ipix_contour = int(value_contour[i][j])
        table_ipix_contour.append(ipix_contour)

    # from index to polar coordinates
    theta, phi = hp.pix2ang(nside, table_ipix_contour)
    # converting these to right ascension and declination in degrees
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
    from astropy.table import Table
    contour_ipix = Table([ra, dec], names=('RA[deg]', 'DEC[deg]'), meta={'ipix': 'ipix table'})

    # setting MOC order
    from math import log
    moc_order = int(log(nside, 2))

    # creating a MOC map from the contour_ipix table
    moc = MOC.from_table(contour_ipix, 'RA[deg]', 'DEC[deg]', moc_order)

    # writing MOC file in fits
    if (save2File):
        moc.write(short_name + '_MOC_' + str(percentage), format='fits')
    return moc


def MOC_confidence_region2D(hpx, percentage, short_name=' ', save2File=False):
    """
        Multi-Order coverage map (MOC) of sky area enclosed within a contour plot
        at a given confidence level.

        Input:
        infile: healpix format
        LVC probability sky map
        percentage: float
        probability percentage of the enclosed area
        short_name: str
        output file name

        Output: fits format
        MOC map named "short_name"_"percentage"

        Remark: for json format change the statement
        "moc.write(short_name+'_MOC_'+str(percentage), format='fits' )" -->
        "moc.write(short_name+'_MOC_'+str(percentage), format='json' )"
        """

    # reading skymap
    npix = len(hpx)
    nside = hp.npix2nside(npix)

    sort = sorted(hpx, reverse=True)
    cumsum = np.cumsum(sort)
    index, value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))

    # finding ipix indices confined in a given percentage
    index_hpx = range(0, len(hpx))
    hpx_index = np.c_[hpx, index_hpx]

    sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
    value_contour = sort_2array[0:index]

    j = 1
    table_ipix_contour = []

    for i in range(0, len(value_contour)):
        ipix_contour = int(value_contour[i][j])
        table_ipix_contour.append(ipix_contour)
    
    # from index to polar coordinates
    theta, phi = hp.pix2ang(nside, table_ipix_contour)
    # converting these to right ascension and declination in degrees
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
    from astropy.table import Table
    contour_ipix = Table([ra, dec], names=('RA[deg]', 'DEC[deg]'), meta={'ipix': 'ipix table'})

    # setting MOC order
    from math import log
    moc_order = int(log(nside, 2))

    # creating a MOC map from the contour_ipix table
    moc = MOC.from_table(contour_ipix, 'RA[deg]', 'DEC[deg]', moc_order)

    # writing MOC file in fits
    if (save2File):
        moc.write(short_name + '_MOC_' + str(percentage), format='fits')
    return moc


def randomDate(start, end, prop):
    return strTimeProp(start, end, '%Y-%m-%d %H:%M:%S', prop)


def strTimeProp(start, end, format, prop):
    """Get a time at a proportion of a range of two formatted times.

        start and end should be strings specifying times formated in the
        given format (strftime-style), giving an interval [start, end].
        prop specifies how a proportion of the interval to be taken after
        start.  The returned time will be in the specified format.
        """

    stime = time.mktime(time.strptime(start, format))
    etime = time.mktime(time.strptime(end, format))
    # print etime

    ptime = stime + prop * (etime - stime)
    # print time.strftime(format, time.localtime(ptime))
    # print time.strftime(format, time.gmtime(ptime))

    return time.strftime(format, time.localtime(ptime))

######################################################

## Extra functions from the use of Center of Pixels ##
## as center of pointings, in a 3D treatment

######################################################

def ModifyCataloguePIX(pix_ra1, pix_dec1, test_time, maxz, prob,cat, FOV, totaldPdV,nside, NewNside, minz):
    
    #To do:
    #for a faster time:
    #subtract the summed pixels
    #chose only visbile pixels / done
    #adjust the number of interations and pixel
    
    #####################
    #pprob = hp.pixelfunc.ud_grade(prob, 64)#power = -2
    #pixel_theta, pixel_phi = hp.pix2ang((hp.npix2nside(len(pprob))), np.arange(len(pprob)))
    
    #pix_ra1 = np.rad2deg(pixel_phi)
    #pix_dec1 = np.rad2deg(0.5 * np.pi - pixel_theta)
    ##################
    
    #Cuts on azimuth angle  (in the probability region)
    observatory = co.EarthLocation(lat=-23.271333 * u.deg,lon=16.5 * u.deg, height=1800 * u.m)
    
    
    frame = co.AltAz(obstime=test_time, location=observatory)
    
    radecs = co.SkyCoord(pix_ra1, pix_dec1, frame='fk5', unit=(u.deg, u.deg))
    thisaltaz = radecs.transform_to(frame)
    
    #pix_alt1 = thisaltaz.alt.value
    
    pix_ra = radecs.ra.value[thisaltaz.alt.value > 90 - (minz)]
    pix_dec = radecs.dec.value[thisaltaz.alt.value > 90 - (minz)]
    #pix_alt = pix_alt1[thisaltaz.alt.value > 90 - (minz)]
    
    dp_Pix_Fov = np.empty(len(pix_ra), dtype=object)
    
    cat_pix = Table([pix_ra, pix_dec, dp_Pix_Fov ], names=('PIXRA', 'PIXDEC', 'PIXFOVPROB'))

    
    ##############################################################
    #Possible:  select the pixels that only have a prob > certain value
    #       To do so attribute for each pix its prob: start with nside initial, put in table prob, reduce resolution, make cut...
    #Note : maybe highest PROBFOV pix is not visible ? ? ? check fullfills requirements
    ###############################################################
    
    dp_dV_FOV = []

    #iteration on chosen pixel to calculate the probability on their field of view using galaxies
    for l in range(0, len(cat_pix)):
        dp_dV_FOV.append(PGalinFOV(prob,cat,cat_pix[l], FOV, totaldPdV, nside,UsePix=True))


    cat_pix['PIXFOVPROB'] = dp_dV_FOV

    ttcat = cat_pix[np.flipud(np.argsort(cat_pix['PIXFOVPROB']))]
    return ttcat



def Get90RegionPixReduced(hpxx, percentage, Nnside):

    nside = 512  # size of map used for contour determination
    hpx = hp.ud_grade(hpxx, nside, power=-2)

    sort = sorted(hpx, reverse=True)
    cumsum = np.cumsum(sort)
    index, value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))
        
    # finding ipix indices confined in a given percentage
    index_hpx = range(0, len(hpx))
    hpx_index = np.c_[hpx, index_hpx]
        
    sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
    value_contour = sort_2array[0:index]
        
    j = 1
    table_ipix_contour = []
            
    for i in range(0, len(value_contour)):
        ipix_contour = int(value_contour[i][j])
        table_ipix_contour.append(ipix_contour)
    #area = len(table_ipix_contour)*hp.nside2pixarea(nside, True)
    # from index to polar coordinates
    theta1, phi1 = hp.pix2ang(nside, table_ipix_contour)
    area = len(table_ipix_contour)*hp.nside2pixarea(nside, True)
                
    # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
    #contour_ipix = Table([ra, dec], names=('RA[deg]', 'DEC[deg]'), meta={'ipix': 'ipix table'})

    #reducing resolution to et a faser execution
    R_ipix = hp.ang2pix(Nnside, theta1, phi1) #list of pixel indices in the new map
    R_ipix = list(set(R_ipix)) #Removing/keeping 1 duplicate from list)
                
    # from index to polar coordinates
    theta, phi = hp.pix2ang(Nnside, R_ipix)
                        
    # converting these to right ascension and declination in degrees
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
                        
    return ra, dec, area


def Get90RegionPixGal(hpxx, percentage, Nside):

    nside = Nside  # size of map used for contour determination
    #hpx = hp.ud_grade(hpxx, nside, power=-2)
    hpx = hpxx
    sort = sorted(hpx, reverse=True)
    cumsum = np.cumsum(sort)
    index, value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))
        
    # finding ipix indices confined in a given percentage
    index_hpx = range(0, len(hpx))
    hpx_index = np.c_[hpx, index_hpx]
        
    sort_2array = sorted(hpx_index, key=lambda x: x[0], reverse=True)
    value_contour = sort_2array[0:index]
        
    j = 1
    table_ipix_contour = []
            
    for i in range(0, len(value_contour)):
        ipix_contour = int(value_contour[i][j])
        table_ipix_contour.append(ipix_contour)
    return table_ipix_contour


