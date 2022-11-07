#!/home/antoniar/python3.6_env/bin/python

#################### WARNING: Turn off sending triggers when testing!!!!!!!!!!!!!

import numpy as np
import scipy as sp
import os
from datetime import datetime
from datetime import timedelta
from datetime import time
import logging
import sys
import six
import voeventparse
import pandas
import requests
import pytz
from lxml import etree as ET
import lofar_maintenance

# SGR Imports
import voeventdb.remote as vr
import voeventdb.remote.apiv1 as apiv1
from voeventdb.remote.apiv1 import FilterKeys

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


################# Runing the 4PiSky Broker #################
# This script is called when the broker receives a VOEvent from the 4PiSky Boker.
# See the example here: https://github.com/4pisky/fourpiskytools/blob/master/examples/listen_for_voevents.sh
# Edit line 14 of this script to HANDLER=./sgrb_swift.py


# Some variables that I need everywhere, so to make life easier...
global ObsMax, ObsMin, AltCut, MaxDwell, Name, email, phoneNumber, Affiliation, ProjectCodeSGRBs, ProjectCodeXRFs, GRB_trig_dur, template, username,LOFARlocation, grblist, delaymins

ObsMax= 120. #time in minutes
ObsMin= 30.
AltCut=20. #minimum elevation for observations in degrees
MaxDwell=8. # maximum dwell time in minutes
MinDwell=4. # minimum dwell time (while testing)
LOFARlocation = EarthLocation(lat=52.9088*u.deg,lon=6.8689*u.deg,height=0.*u.m)
CalObsT = 10.
delaymins=10.

Name = "Name"
username = 'username'
email = "email@adress"
phoneNumber = "tel"
Affiliation = "ASTRON"
ProjectCodeSGRBs = 'LC'
ProjectCodeXRFs = 'LC'
ProjectCodeSGRs = 'LC'
GRB_trig_dur = 2.
Rate_Signif = 10.
template = 'sgrb_template.xml'
grblist = 'triggeredGRBs_TrigIDs.txt'
sgrlist = 'triggeredSGRs_TrigIDs.txt'

#List of SGRs to trigger on, and longer delay
sgr_name_list = ['sgr1935','SGR1900+14','SGR1806-20',
                 'SGR_0501+4516','SGR0418+5729',
                 '4U_0142+61','SGR1801-23',
                 'sgr1818','sgr1845',
                 'sgr1834']

sgr_coord_list = [[293.73199,21.89673],[286.8041667,9.3261111],[272.1638333,-20.4109722],
                  [75.283,45.2752778],[64.6060,57.4890],
                  [26.5925417,61.7510556],[270.2454167,-22.9466667],
                  [274.5154,-16.1255],[281.2278333,-2.9480833],
                 [278.717,-8.766]]
delaymins_sgr=300.

#################
# This code is edited from Tim's example here: https://github.com/4pisky/fourpiskytools/blob/master/examples/process_voevent_from_stdin_2.py
# Some of these can be used directly from the 4pisky tutorials but not all.

logging.basicConfig(filename='sgrb_swift.log',level=logging.INFO)
logger = logging.getLogger('notifier')
logger.handlers.append(logging.StreamHandler(sys.stdout))

def main():
    if six.PY2:
        stdin = sys.stdin.read()
    else:
        stdin = sys.stdin.buffer.read()

    v = voeventparse.loads(stdin)

    handle_voevent(v)
    return 0

def handle_voevent(v):

    if is_grb(v):
        logger.info(format(v.attrib['ivorn']))
        logger.info('is GRB')
        RA, Dec, time, parameters, ivorn = handle_grb(v) # You may need to edit this for non-GRBs and non-Swift

        f= open(grblist,"a+")
        if str(parameters[None]['TrigID']['value']) in f.read():
            logger.info("GRB already processed " + parameters[None]['TrigID']['value'])
            f.close()
            return
        else:
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+','+format(time)+','+str(RA)+','+str(Dec)+','+format(parameters[None]['Integ_Time']['value'])+','+format(parameters[None]['Rate_Signif']['value']))

            f.write(str(parameters[None]['TrigID']['value'])+','+format(parameters[None]['Integ_Time']['value'])+','+format(parameters[None]['Rate_Signif']['value'])+','+str(time)+'\n')
            f.close()
            logger.info('Written data to: '+grblist)
            return

    if is_slew_to_grb(v):
        logger.info(format(v.attrib['ivorn']))
        logger.info('Swift is slewing to source')
        RA, Dec, parameters, ivorn = handle_grb2(v)
        logger.info(parameters[None]['TrigID']['value'])
        f= open(grblist,"r")
        logger.info('reading '+grblist)
        lines=f.readlines()

        tmp=0
        for line in lines:
            #logger.info(line)
            params=line.split(',')
            if params[0]==str(parameters[None]['TrigID']['value']):
                    tmp=1
                    integtime=float(params[1])
                    ratesignif=float(params[2])
                    tmpdate = str(params[3])
                    tmpdate = tmpdate.split('+')
                    tmpdate = tmpdate[0]
                    trigtime=datetime.strptime(tmpdate, '%Y-%m-%d %H:%M:%S.%f').replace(tzinfo=pytz.utc)
        f.close()

        if tmp==0:
            logger.info("Not a triggered GRB")
            return

        now = datetime.utcnow().replace(tzinfo=pytz.utc)
        logger.info(now)
        logger.info(trigtime)
        logger.info(now - trigtime)

        if now-trigtime < timedelta(minutes=delaymins):
            filterburst(RA, Dec, now, parameters, v, integtime, ratesignif)
            logger.info(delaymins, float(parameters[None]['Wait_Time']['value'])/60. < delaymins)
            if float(parameters[None]['Wait_Time']['value'])/60. < delaymins:
                logger.info('Slew time check passed')
            else:
                logger.info('Delayed slew,  >'+str(delaymins)+' min after GRB')
        else:
                logger.info('Delayed slew,  >'+str(delaymins)+' min after GRB')

    if is_sgr_2(v):
        logger.info(format(v.attrib['ivorn']))
        logger.info('is SGR burst')
        RA, Dec, time, parameters, ivorn = handle_grb(v)

        now = datetime.utcnow().replace(tzinfo=pytz.utc)
        logger.info(now)
        logger.info(time)
        logger.info(now - time)

        if now-time < timedelta(minutes=delaymins_sgr):
            logger.info('Passed time check')
            filterSGRs(RA, Dec, now, parameters, v)

        else:
            logger.info('Failed time check')


    elif is_swift_pointing(v):
        handle_pointing(v)

    elif is_ping_packet(v):
        handle_ping_packet(v)

    else:
        handle_other(v)



def is_grb(v):
    ivorn = v.attrib['ivorn']
    if ivorn.find("ivo://nasa.gsfc.gcn/SWIFT#BAT_GRB_Pos") == 0:
        return True
    return False

def is_slew_to_grb(v):
    ivorn = v.attrib['ivorn']
    if ivorn.startswith("ivo://nasa.gsfc.gcn/SWIFT#SC_Slew"):
        return True
    return False

def is_sgr(v):

    ivorn = v.attrib['ivorn']
    if ivorn.find("ivo://nasa.gsfc.gcn/SWIFT#BAT_Trans_Pos") == 0:
        for sgr_name in sgr_name_list:
            if sgr_name in str(v.Why.Inference.Name):
                return True
            # Can add in coordinate-based triggering in future here
            else:
                return False
    return False


def is_sgr_2(v):
    err = 0.1
    ivorn = v.attrib['ivorn']
    if ivorn.find("ivo://nasa.gsfc.gcn/SWIFT#BAT_GRB_Pos") == 0:
        coords = voeventparse.get_event_position(v)
        ra=coords.ra
        dec=coords.dec

        for i, sgr_coord in enumerate(sgr_coord_list):
            if sgr_coord[0]-err < ra < sgr_coord[0]+err:
                if sgr_coord[1]-err < dec < sgr_coord[1]+err:
                    return True
                else:
                    return False
            else:
                return False

    elif ivorn.find("ivo://nasa.gsfc.gcn/SWIFT#BAT_Known") == 0:
        coords = voeventparse.get_event_position(v)
        ra=coords.ra
        dec=coords.dec

        for i, sgr_coord in enumerate(sgr_coord_list):
            if sgr_coord[0]-err < ra < sgr_coord[0]+err:
                if sgr_coord[1]-err < dec < sgr_coord[1]+err:
                    return True
                else:
                    return False
            else:
                return False
    else:
        return False



def is_swift_pointing(v):
    ivorn = v.attrib['ivorn']
    if ivorn.startswith("ivo://nasa.gsfc.gcn/SWIFT#Point_Dir_"):
        return True
    return False

def is_ping_packet(v):
    ivorn = v.attrib['ivorn']
    role = v.attrib['role']
    ivorn_tail = ivorn.rsplit('/', 1)[-1]
    stream = ivorn_tail.split('#')[0]
    if ( role == voeventparse.definitions.roles.test and
         stream == 'TEST'
        ):
        return True
    return False

def handle_grb(v):
    ivorn = v.attrib['ivorn']
    coords = voeventparse.get_event_position(v)
    ra=coords.ra
    dec=coords.dec
    parameters = voeventparse.convenience.pull_params(v)
    time = voeventparse.convenience.get_event_time_as_utc(v, index=0) # this is a datetime.datetime object
    return ra, dec, time, parameters, ivorn

def handle_grb2(v):
    logger.info('in handle grb')
    ivorn = v.attrib['ivorn']
    logger.info("VOEvent received. IVORN: "+ivorn)
    coords = voeventparse.get_event_position(v)
    ra=coords.ra
    dec=coords.dec
    logger.info(coords)
    parameters = voeventparse.convenience.pull_params(v)
    return ra, dec, parameters, ivorn


def handle_pointing(v):
    ivorn = v.attrib['ivorn']
    coords = voeventparse.get_event_position(v)
    text = "Swift repointing, coords are {}".format(coords)

def handle_ping_packet(v):
   #logger.info("Packet received matches 'ping packet' filter.")
   handle_other(v)

def handle_other(v):
    ivorn = v.attrib['ivorn']
    #logger.info("VOEvent received. IVORN: "+ivorn)


################# Filter the possible short GRBs and send alert #################

def filterburst(RA, Dec, now, parameters, v, integtime, ratesignif):
    logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'In XRF filtering code')
    if ratesignif > Rate_Signif:
        filterXRFs(RA, Dec, now, parameters,v) # changed to use the time now rather than the trigger time
    elif integtime < GRB_trig_dur:
        filterSGRBs(RA, Dec, now, parameters,v)
    else:
        logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Failed trigger criteria')
    return

def filterXRFs(RA, Dec, time, parameters,v):
    # Trig_dur needs to be less than the input, readFromEvent needs to be checked.
    # calculate time of event in mjds
    logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Rate significance check passed. '+str(parameters[None]['Rate_Signif']['value']))
    # Source needs to stay above horizon
    index = horizon(time,RA, Dec)
    logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Index from altitude check: '+str(index))
    if index != 0:
        # we need a calibrator source
        calibrator = find_cal(time,RA, Dec,index)
        if calibrator['Calibrators'] != 'None':
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Calibrator found: '+str(calibrator['Calibrators'])+','+str(calibrator['CalRA'])+','+str(calibrator['CalDec']))
            xmlname = writeXMLfiles(format(parameters[None]['TrigID']['value']),time,RA, Dec,calibrator,index,ProjectCodeSGRBs)
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'XML file written to: '+xmlname)
            sendXMLtoLOFAR(xmlname)
        else:
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Failed to find a calibrator')
    else:
        logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Not above horizon within time constraints')

def filterSGRBs(RA, Dec, time, parameters,v):
    # Trig_dur needs to be less than the input, readFromEvent needs to be checked.
    # calculate time of event in mjds
    logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Integration time check passed. '+str(parameters[None]['Integ_Time']['value'])+','+str(GRB_trig_dur))

    logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+str(time)+','+str(RA)+','+str(Dec))

    # Source needs to stay above horizon
    index = horizon(time,RA, Dec)
    logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Index from altitude check: '+str(index))
    if index != 0:
        # we need a calibrator source
        calibrator = find_cal(time,RA, Dec,index)
        if calibrator['Calibrators'] != 'None':
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Calibrator found: '+str(calibrator['Calibrators'])+','+str(calibrator['CalRA'])+','+str(calibrator['CalDec']))
            xmlname = writeXMLfiles(format(parameters[None]['TrigID']['value']),time,RA, Dec,calibrator,index,ProjectCodeSGRBs)
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'XML file written to: '+xmlname)
            sendXMLtoLOFAR(xmlname)
        else:
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Failed to find a calibrator')
    else:
        logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Not above horizon within time constraints')

def filterSGRs(RA, Dec, time, parameters,v):
    # Needs a time-check, and also a check for last known LOFAR trigger
    # Look for swift detected bursts in past 2 weeks
    my_filters = {FilterKeys.role: 'observation',
    FilterKeys.cone: '['+str(RA)+','+str(Dec)+',0.1]',
    FilterKeys.authored_since: time - timedelta(days=7),
    FilterKeys.authored_until: time,
    }

    rel_sgr_ivorns = apiv1.list_ivorn(filters=my_filters,order=apiv1.OrderValues.author_datetime_desc,n_max=10)
    if len(rel_sgr_ivorns) > 0:
        logger.info('Found'+str(rel_sgr_ivorns)+' relevant IVORNs.')

        index = horizon(time,RA, Dec)
        logger.info('TrigID: '+format(parameters[None]['TrigID']['value'])+', '+'Index from altitude check: '+str(index))

        if index != 0:
             # we need a calibrator source
            calibrator = find_cal(time,RA, Dec,index)
            if calibrator['Calibrators'] != 'None':
                logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', '+'Calibrator found: '+str(calibrator['Calibrators'])+','+str(calibrator['CalRA'])+','+str(calibrator['CalDec']))
                xmlname = writeXMLfiles(format(parameters[None]['TrigID']['value']).replace("-", "0"),time,RA, Dec,calibrator,index,ProjectCodeSGRs)
                logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', '+'XML file written to: '+xmlname)
                sendXMLtoLOFAR(xmlname)
            else:
                logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', '+'Failed to find a calibrator')
        else:
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', '+'Not above horizon within time constraints')
    else:
        logger.info('Not enough relevant ivorns.')
        return
################# Check the elevation of the source #################

def horizon(time,RA,Dec):
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', in horizon function')
            az0=calcAltAz(time+timedelta(minutes=MaxDwell),RA,Dec) # altitude at max dwell time
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', az0')
            az0b=calcAltAz(time+timedelta(minutes=MaxDwell),RA,Dec) # altitude at min dwell time
            az2=calcAltAz(time+timedelta(minutes=(MaxDwell+ObsMax)),RA,Dec) # altitude at maximum end time
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', az2')
            az1=calcAltAz(time+timedelta(minutes=(MaxDwell+ObsMin)),RA,Dec) # altitude at minimum end time
            logger.info('TrigID: '+format(parameters[None]['TrigID']['value']).replace("-", "0")+', az1')

            if az0<AltCut: # the source remains below minimum elevation by the max dwell time
                return 0
            if az0b>AltCut: # source is above min elevation at start of observation
                if az1<AltCut: # the source sets before the minimum observation duration so don't observe
                    return 0
                if az2>AltCut: # the source is above the minimum for the whole observation, observe for full observation
                    return 1
                if az1>AltCut and az2<AltCut: # the source drops between the max and min time, so only observe for the minimum
                    return 2
            if az0>AltCut and az0b<AltCut: # source rises above min elevation during dwell time, wait till dwell time to start observations
                if az1<AltCut: # the source sets before the minimum observation duration so don't observe
                    return 0
                if az2>AltCut: # the source is above the minimum for the whole observation, so observe for full duration
                    return 3
                if az1>AltCut and az2<AltCut: # the source drops between the max and min time, so only observe for the minimum
                    return 4

    # if 1 is returned, ask for the full time immediately, no constraints
    # if 2 is returned just ask for minimum (until smarter method,
    # this also helps if there is a delay) and no constraints.
    # the max dwell time is considered here - if source would set during
    # the latest possible observation, we will not observe
    # need to keep the dwell time short for now.
    # if 0 is returned, don't observe
    # 3 returned, start full observation after dwell time
    # 4 returned, ask for minimum duration after the dwell time



def find_cal(time,RA, Dec, index):
    calibrators=pandas.read_csv('calibrators.csv',sep=',',header=0)
    separation=648000.
    for index2, cal in calibrators.iterrows():
        c1 = SkyCoord(RA,Dec,unit='deg')
        c2 = SkyCoord(float(cal.ra),float(cal.dec),unit='deg')
        septmp = c1.separation(c2)
        septmp = septmp.arcsecond
        startT = time+timedelta(minutes=MinDwell)
        if index == 1:
            timeStart = startT+timedelta(minutes=ObsMax+MaxDwell+2.)
            timeEnd = timeStart+timedelta(minutes=CalObsT)
        elif index == 2:
            timeStart = startT+timedelta(minutes=ObsMin+MaxDwell+2.)
            timeEnd = timeStart+timedelta(minutes=CalObsT)
        elif index == 3:
            timeStart = startT+timedelta(minutes=ObsMax+MaxDwell+2.)
            timeEnd = timeStart+timedelta(minutes=CalObsT)
        elif index == 4:
            timeStart = startT+timedelta(minutes=ObsMin+MaxDwell+2.)
            timeEnd = timeStart+timedelta(minutes=CalObsT)
        else:
            return {          'Calibrators':'None',
                              'CalSep':0,
                              'CalRA':0,
                              'CalDec':0}

        alt_start = calcAltAz(timeStart,RA,Dec) # check altitude of calibrator at start of observation
        alt_end = calcAltAz(timeEnd,RA,Dec) # check altitude of calibrator at end of observation
        if septmp < separation and alt_start > AltCut and alt_end > AltCut:
            # if calibrator is above min elevation for duration and is closer than previous calibrator, it becomes the optimum calibrator
            separation = septmp
            optcal = cal.src
            optra = cal.ra
            optdec = cal.dec
    if optcal:
        return {'Calibrators':optcal,
                              'CalSep':(separation/(60.*60.)),
                              'CalRA':optra,
                              'CalDec':optdec}
    else:
        return {'Calibrators':'None',
                              'CalSep':0,
                              'CalRA':0,
                              'CalDec':0}





################# Write the XML for LOFAR #################
def generateXML(GRB,RA,Dec,CalRA,CalDec,start,maxDur,minDur,end,calname,startCal,endCal,stations,ProjectCode):
    # this is particularly ugly so I suspect there is a better way? Also specific for the sgrb_template.xml

    tree = ET.parse(template)
    root = tree.getroot()

    if len(stations) != 0:
        for station in stations:
            for child in root.iter('specification'):
                for child2 in child.findall('activity'):
                    for child3 in child2.findall('observation'):
                        for child4 in child3.findall('stationSelectionSpecification'):
                            for child5 in child4.findall('stationSelection'):
                                for child6 in child5.findall('stations'):
                                    for child7 in child6.findall('station'):
                                        for child8 in child7.findall('name'):
                                            if child8.text==station:
                                                child6.remove(child7)
                                                child7.remove(child8)
                                                logger.info(station+' is removed due to maintenance')

    for child in root.iter('userName'):
        child.text = username
    for child in root.findall('contactInformation'):
        for child2 in child.findall('name'):
            child2.text = Name
        for child2 in child.findall('email'):
            child2.text = email
        for child2 in child.findall('phoneNumber'):
            child2.text = phoneNumber
        for child2 in child.findall('affiliation'):
            child2.text = Affiliation
    for child in root.findall('projectReference'):
        for child2 in child.findall('ProjectCode'):
            child2.text = ProjectCode
    for child in root.findall('specification'):
        for child2 in child.findall('projectReference'):
            for child3 in child2.findall('ProjectCode'):
                child3.text = ProjectCode
        for child2 in child.iter('container'):
            for child3 in child2.findall('folder'):
                for child4 in child3.findall('name'):
                    if child4.text == "TestOpened_1":
                        child4.text = 'Trig'+str(GRB)
        for child2 in child.findall('activity'):
            for child3 in child2.findall('observation'):
                for child4 in child3.findall('timeWindowSpecification'):
                    for child5 in child4.findall('minStartTime'):
                        child5.text = start
                    for child5 in child4.findall('maxEndTime'):
                        child5.text = end
                    for child5 in child4.findall('duration'):
                        for child6 in child5.findall('duration'):
                            if child6.text == "PT7200S":
                                child6.text="PT"+str(int(maxDur)*60)+"S"
                for child4 in child3.findall('name'):
                    if child4.text == "3C295/1/TO":
                        child4.text=calname+'/1/TO'
                        for child6 in child3.findall('timeWindowSpecification'):
                            for child5 in child6.findall('minStartTime'):
                                child5.text = startCal
                            for child5 in child6.findall('maxEndTime'):
                                child5.text = endCal
                                for child5 in child4.findall('duration'):
                                    for child6 in child5.findall('duration'):
                                        if child6.text == "PT120S":
                                            child6.text="PT"+str(int(CalObsT)*60)+"S"
                for child4 in child3.findall('description'):
                    if child4.text == "3C295/1/TO (Target Observation)":
                        child4.text=calname+'/1/TO (Target Observation)'
        for child2 in child.findall('activity'):
            for child3 in child2.findall('measurement'):
                for child4 in child3.iter('name'):
                    if child4.text == "Target":
                        for child5 in child3.findall('ra'):
                            child5.text = str(RA)
                        for child5 in child3.findall('dec'):
                            child5.text = str(Dec)
                    if child4.text == "Calibrator":
                        for child5 in child3.findall('ra'):
                            child5.text = str(CalRA)
                        for child5 in child3.findall('dec'):
                            child5.text = str(CalDec)

    tree.write(str(GRB)+'.xml',xml_declaration=True)
    return str(GRB)+'.xml'

def writeXMLfiles(GRB,time,RA, Dec,calibrator,index,ProjectCode):
    if index != 0:
        starttime = time+timedelta(minutes=MinDwell)
        if index == 1:
            start=starttime.strftime("%Y-%m-%dT%H:%M:%S")
            maxDur = ObsMax
            endtime =  starttime+timedelta(minutes=maxDur)
            end = endtime.strftime("%Y-%m-%dT%H:%M:%S")
        elif index == 2:
            start=starttime.strftime("%Y-%m-%dT%H:%M:%S")
            maxDur = ObsMin
            endtime =  starttime+timedelta(minutes=maxDur)
            end = endtime.strftime("%Y-%m-%dT%H:%M:%S")
        elif index == 3:
            starttime=starttime+timedelta(seconds=MaxDwell*60.)
            start=starttime.strftime("%Y-%m-%dT%H:%M:%S")
            maxDur = ObsMax
            endtime =  starttime+timedelta(minutes=maxDur)
            end = endtime.strftime("%Y-%m-%dT%H:%M:%S")
        elif index == 4:
            starttime=starttime+timedelta(seconds=MaxDwell*60.)
            start=starttime.strftime("%Y-%m-%dT%H:%M:%S")
            maxDur = ObsMin
            endtime =  starttime+timedelta(minutes=maxDur)
            end = endtime.strftime("%Y-%m-%dT%H:%M:%S")

        startCal = endtime+timedelta(minutes=2.)
        startCal = startCal.strftime("%Y-%m-%dT%H:%M:%S")
        endCaltmp = endtime+timedelta(minutes=20.)
        endCal = endCaltmp.strftime("%Y-%m-%dT%H:%M:%S")
        stations = check_LOFAR_maintenance(starttime,endCaltmp)
        xmlname = generateXML(GRB,RA,Dec,calibrator['CalRA'],calibrator['CalDec'],start,maxDur,ObsMin,end,calibrator['Calibrators'],startCal,endCal,stations,ProjectCode)
        return xmlname

#test generateXML
#generateXML(123.4,56.7,10.2,45.3,"2018-11-23T15:21:44","7200","720")

def sendXMLtoLOFAR(xmlname): # code received from Auke last year, needs updating...
    os.system("echo curl --insecure --data-binary @"+xmlname+" --netrc 'https://proxy.lofar.eu/rtrest/triggers/?format=json'")
# Uncomment line below to send xml to LOFAR
    os.system("curl --insecure --data-binary @"+xmlname+" --netrc 'https://proxy.lofar.eu/rtrest/triggers/?format=json'")

def check_LOFAR_maintenance(start,end):
    start=start.strftime("%Y-%m-%d %H:%M:%S")
    end=end.strftime("%Y-%m-%d %H:%M:%S")
    stations=lofar_maintenance.getMaintenance(start,end)
    return stations

def calcAltAz(time,RA,Dec):
    time = Time(time.strftime("%Y-%m-%d %H:%M:%S"))
    position = SkyCoord(RA,Dec,unit='deg')
    AltAzPos = position.transform_to(AltAz(obstime=time,location=LOFARlocation))
    return AltAzPos.alt.deg

if __name__ == '__main__':
    sys.exit(main())


#v=voeventparse.loads('<?xml version=\'1.0\' encoding=\'UTF-8\'?>\n<voe:VOEvent xmlns:voe="http://www.ivoa.net/xml/VOEvent/v2.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ivorn="ivo://nasa.gsfc.gcn/SWIFT#BAT_Known_Pos_-1637829720-907" role="observation" version="2.0" xsi:schemaLocation="http://www.ivoa.net/xml/VOEvent/v2.0  http://www.ivoa.net/xml/VOEvent/VOEvent-v2.0.xsd"><Who><AuthorIVORN>ivo://nasa.gsfc.tan/gcn</AuthorIVORN><Author><shortName>VO-GCN</shortName><contactName>Scott Barthelmy</contactName><contactPhone>+1-301-286-3106</contactPhone><contactEmail>scott.barthelmy@nasa.gov</contactEmail></Author><Date>2021-09-16T02:39:55</Date><Description>This VOEvent message was created with GCN VOE version: 15.08 25jun21</Description></Who><What><Param name="Packet_Type" value="141"/><Param name="Pkt_Ser_Num" value="6600"/><Param name="TrigID" value="-1637829720" ucd="meta.id"/><Param name="Event_TJD" value="19472" unit="days" ucd="time"/><Param name="Event_SOD" value="70678.00" unit="sec" ucd="time"/><Param name="Event_Inten" value="112" unit="cts" ucd="phot.count;em.gamma.soft"/><Param name="Integ_Time" value="0.016" unit="sec" ucd="time.interval"/><Param name="Phi" value="-58.35" unit="deg" ucd="pos.az.azi"/><Param name="Theta" value="30.90" unit="deg" ucd="pos.az.zd"/><Param name="Soln_Status" value="0x808"/><Param name="Misc_flags" value="0x40000000"/><Param name="Image_Signif" value="14.00" unit="sigma" ucd="stat.snr"/><Param name="Cat_Num" value="10762" ucd="meta.id"/><Group name="Solution_Status"><Param name="Point_Source" value="false"/><Param name="VERY_Lo_Image_Signif" value="true"/><Param name="Target_in_Flt_Catalog" value="true"/><Param name="Target_in_Gnd_Catalog" value="false"/><Param name="Near_Bright_Star" value="false"/><Param name="Spatial_Prox_Match" value="false"/><Param name="Temporal_Prox_Match" value="false"/><Param name="Test_Submission" value="false"/></Group><Group name="Misc_Flags"><Param name="Values_Out_of_Range" value="false"/><Param name="Near_Bright_Star" value="false"/><Param name="Err_Circle_in_Galaxy" value="false"/><Param name="Galaxy_in_Err_Circle" value="false"/></Group><Param name="Coords_Type" value="1" unit="dn"/><Param name="Coords_String" value="source_object"/><Group name="Obs_Support_Info"><Description>The Sun and Moon values are valid at the time the VOEvent XML message was created.</Description><Param name="Sun_RA" value="174.00" unit="deg" ucd="pos.eq.ra"/><Param name="Sun_Dec" value="2.59" unit="deg" ucd="pos.eq.dec"/><Param name="Sun_Distance" value="116.49" unit="deg" ucd="pos.angDistance"/><Param name="Sun_Hr_Angle" value="-8.00" unit="hr"/><Param name="Moon_RA" value="295.26" unit="deg" ucd="pos.eq.ra"/><Param name="Moon_Dec" value="-25.27" unit="deg" ucd="pos.eq.dec"/><Param name="MOON_Distance" value="47.23" unit="deg" ucd="pos.angDistance"/><Param name="Moon_Illum" value="74.62" unit="%" ucd="arith.ratio"/><Param name="Galactic_Long" value="57.25" unit="deg" ucd="pos.galactic.lon"/><Param name="Galactic_Lat" value="0.81" unit="deg" ucd="pos.galactic.lat"/><Param name="Ecliptic_Long" value="300.63" unit="deg" ucd="pos.ecliptic.lon"/><Param name="Ecliptic_Lat" value="42.84" unit="deg" ucd="pos.ecliptic.lat"/></Group><Description>Type=141: The sub-sub-threshold Swift-BAT trigger position notice.</Description></What><WhereWhen><ObsDataLocation><ObservatoryLocation id="GEOLUN"/><ObservationLocation><AstroCoordSystem id="UTC-FK5-GEO"/><AstroCoords coord_system_id="UTC-FK5-GEO"><Time unit="s"><TimeInstant><ISOTime>2021-09-15T19:37:58.00</ISOTime></TimeInstant></Time><Position2D unit="deg"><Name1>RA</Name1><Name2>Dec</Name2><Value2><C1>293.7450</C1><C2>21.8909</C2></Value2><Error2Radius>0.0666</Error2Radius></Position2D></AstroCoords></ObservationLocation></ObsDataLocation><Description>The RA,Dec coordinates are of the type: source_object.</Description></WhereWhen><How><Description>Swift Satellite, BAT Instrument</Description><Reference uri="http://gcn.gsfc.nasa.gov/swift.html" type="url"/></How><Why importance="0.90"><Inference probability="0.98"><Name>sgr1935p2154</Name><Concept>process.variation.burst;em.gamma</Concept></Inference></Why><Description>\n  </Description></voe:VOEvent>')
#handle_voevent(v)
