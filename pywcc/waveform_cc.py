#!/anaconda/bin/python
"""
Copyright (c) 2017, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Dennise Templeton, templeton4@llnl.gov.
LLNL-CODE-737623. All rights reserved.
    
This file is part of PyWCC. For details, see
https://github.com/templetond/pywcc.
    
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
(as published by the Free Software Foundation) version 2.1
dated February 1999.
    
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
terms and conditions of the GNU General Public License for more
details.
    
You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free
Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
MA 02111-1307 USA
    
Please also see the NOTICE and LICENSE files for our notice and the LGPL.
    
PyWCC 0.1.1
Functions for single-component or multiple-component network waveform cross
correlation of seismic data. Cross-correlates waveform templates with 
continuous seismic waveform data, associates the detections, identifies the 
template with the highest cross correlation coefficent, and outputs a
catalog of detections above an absolute cross correlation threshold value.
"""

from obspy import read
from obspy import UTCDateTime

import os
import sys
import glob
import datetime
import warnings
import argparse
import numpy as np
from math import floor
from distutils.util import strtobool
from scipy.signal import fftconvolve
from scipy.ndimage.filters import uniform_filter


def WvfmSimilarity(st_data, st_templates_list, Tshift):
    """
    Compares template waveforms to longer continuous seismic data streams.
    Returns the cross correlation coefficient, and an estimate of origin time 
    for new detection. Both template and data streams should have identical 
    trace ids to be compared. Accepts one-component or multiple-component
    streams. All components will be shifted together.
    
    :type st_data: obspy.core.stream.Stream
    :param st_data: A Stream object containing the continuous data
    :type st_templates: List of :class:'~obspy.core.stream.Stream'
    :param st_templates: List of template streams
    :param Tshift: time between origin time and template start time
    :returns: correlation coefficient between template and data streams in 
        steps of one sample, and an estimate of origin time
    """
    
    # Remove data streams with zero amplitude traces
    for tr_ in st_data:
        if np.sum(tr_.data == 0):
            msg = 'WARNING: Data trace is all zeros, removing. '
            warnings.warn(msg)
            st_data.remove(tr_)

    # Calculate time values
    time_beg = UTCDateTime(st_data[0].stats.starttime)
    time_mod = time_beg + Tshift

    # Extract traces from template list
    for st_templates in st_templates_list:
        # List of traceIDs in template stream
        ids = [tr.id for tr in st_templates]
        
        # Remove template traceID if no identical traceID in data stream
        for id_ in reversed(ids):
            if not st_data.select(id=id_):
                print 'WARNING: Data trace missing. Removing template: ', id_
                ids.remove(id_)

        cc = None
        # Iterate over each validated template traceID stream
        for id_ in ids:
            
            tr1 = st_data.select(id=id_)[0]
            tr2 = st_templates.select(id=id_)[0]
            if len(tr1) > len(tr2):
                data_short = tr2.data
                data_long = tr1.data
            else:
                data_short = tr1.data
                data_long = tr2.data

            # Correlation coefficient determined from convolution
            corr_top = fftconvolve((data_long - np.mean(data_long)), \
                                   (data_short - np.mean(data_short))[::-1], \
                                   mode='valid')
            lds = len(data_short)
            ldl = len(data_long)
            ld = ldl-lds+1
            val_a = uniform_filter(data_long, lds, mode='constant')
            val_b = uniform_filter(data_long*data_long, lds, mode='constant')
            sind = int(floor((lds+1)/2))
            std_dl = (np.sqrt(val_b - val_a*val_a))[lds-sind:-sind+1]
            std_ds = data_short.std()
            corr_bot = std_dl * std_ds * lds
            cor_coef = corr_top / corr_bot

            try:
                cc += cor_coef
            except TypeError:
                cc = cor_coef
            except ValueError:
                cc = None
                break

        if cc is None:
            msg = 'WARNING: No valid data-template pairs. \
                  Correlation has not run.'
            warnings.warn(msg)

        # Estimated Origin Time of Potential Detection
        st_delta = st_data[0].stats.delta
        orig_t = [(time_mod + (st_delta * ind)) for ind in xrange(0,len(cc),1)]
    
    #  Mean correlation coefficient over all streams
    if ids:
        coef = np.array(cc) / len(ids)
        return coef, orig_t
    else:
        return 0, 0


def WSwrapper(pre_ot, post_ot, bp_low, bp_high, dec, cc_thr, sta_l, comp_l,
             data_dir='./ContinuousData', out_dir='./OutputFiles',
             templ_dir='./TemplateData', templ_file='./TemplateCatalog',
             verbose=False):
    """
    :param pre_ot: Seconds to cut template waveform prior to \
        reference origin time
    :param post_ot: Seconds to cut template waveform after \
        reference origin time
    :param bp_low: Low corner pass band frequency
    :param bp_high: High corner pass band frequency
    :param dec: Decimation rate
    :param sta_l: List of station name strings
    :param comp_l: List of component name strings
    :param cc_thr: Cross Correlation threshold for output catalog
    :type cc_thr: str
    :type data_dir: str
    :param data_dir: Continuous data directory path. Defaults to run directory
    :type out_dir: str
    :param out_dir: Path to output directory. Defaults to run directory.
    :type templ_dir: str
    :param templ_dir: Path to template directory. Defaults to run directory
    :type templ_file: str
    :param templ_file: Path to template file. Defaults to run directory. \
        Template file will have form: \
        TEMPLATE_NAME YYYY/MM/DD HH:MM:SS.MMM LAT LON DEPTH
    :param verbose: If True, will print program messages to screen
        
    """

    # Issue input file warnings
    try:
        strtobool(verbose)
    except ValueError:
       sys.exit('WARNING: Non-boolean type specified for verbose argument '
                'in input file. Change to True / False.')

    try:
        pre_ot = float(pre_ot)
        post_ot = float(post_ot)
        bp_low = float(bp_low)
        bp_high = float(bp_high)
        cc_thr = float(cc_thr)
    except ValueError:
        sys.exit('WARNING: Invalid value(s) in input file. Arguments pre_ot, '
                 'post_ot, bp_low, bp_high, dec, and cc_thr should be integers '
                 'or floats.')

    try:
        dec = int(dec)
    except ValueError:
        sys.exit('WARNING: Invalid value of dec in input file. '
                 'Argument should be an integer.')

    # Absolute file path from relative, if necessary
    templ_file_full = os.path.abspath(templ_file)
    templ_dir_full = os.path.abspath(templ_dir)
    data_dir_full = os.path.abspath(data_dir)
    out_dir_full = os.path.abspath(out_dir)
    
    # Create output directory, if necessary
    if not os.path.exists(out_dir_full):
        os.makedirs(out_dir_full)

    # If verbose mode: Raw output file name and directory
    if verbose == 'True':
        ymds = os.path.basename(data_dir_full)
        trig_day_file_nm = 'RawOutput_' + ymds
        out_path_nm = os.path.join(out_dir_full,trig_day_file_nm)
        try:
            os.remove(out_path_nm)
        except OSError:
            pass

    # Read template file
    try:
        g = open(templ_file_full,'r')
    except IOError:
        print 'WARNING: Could not read template catalog file: ', templ_file_full
        sys.exit()

    with g:
        tf = g.readlines()
    g.close()

    # Get continuous data
    if verbose == 'True':
        print datetime.datetime.now(), 'Getting continuous data: ', ymds
    try:
        os.chdir(data_dir_full)
    except OSError:
        print 'WARNING: No such data directory: ', data_dir_full
        sys.exit()

    for wvfm_cont in glob.glob('*'):
        # Read in data and save in 'st_c'
        if not 'st_c' in locals():
            try:
                st_c = read(wvfm_cont)
            except TypeError:
                pass
        else:
            try:
                st_c += read(wvfm_cont)
            except TypeError:
                pass

    # If no data, give errors
    if not 'st_c' in locals():
        print 'WARNING: No data files in directory specified: ', data_dir_full
        sys.exit()
    elif len(st_c) == 0:
        print 'WARNING: Empty data files in directory specified: ',data_dir_full
        sys.exit()

    # Process data
    st_c.merge(method=1,fill_value='interpolate',interpolation_samples=200)
    st_c.decimate(dec, no_filter=False, strict_length=False)
    utc_beg_c = UTCDateTime(st_c[0].stats.starttime)
    utc_end_c = UTCDateTime(st_c[0].stats.endtime)

    st_c.detrend(type = 'linear')
    st_c.detrend(type = 'demean')
    st_c.filter('bandpass',freqmin=bp_low,freqmax=bp_high,zerophase=True)
    st_c.trim(starttime=utc_beg_c,endtime=utc_end_c,pad=True,fill_value=0)
    
    sta_l = sta_l.split()
    comp_l = comp_l.split()
    raw_lines = []
    raw_lines_cc = []
    # Iterate over each template event
    for line in tf:
        
        # Date-time of template event
        lfields = line.split()
    
        det = lfields[0]
        date_sep = lfields[1].split('/')
        time_sep = lfields[2].split(':')

        year_s = date_sep[0]
        mon_s = date_sep[1]
        day_s = date_sep[2]
        hr_s = time_sep[0]
        min_s = time_sep[1]
        sec_f = float(time_sep[2])
        sec_fields = time_sep[2].split('.')
        sec_is = sec_fields[0]

        year_n = int(year_s)
        mon_n = int(mon_s)
        day_n = int(day_s)
        hr_n = int(hr_s)
        min_n = int(min_s)
        sec_n = int(sec_is)
        msec_n = int((sec_f - sec_n) * 1000000)

        # UTC seconds pre/post zero detection time to cut template
        sec_beg_t=datetime.datetime(year_n,mon_n,day_n,hr_n,min_n,sec_n,msec_n)\
                 - datetime.timedelta(seconds=abs(pre_ot))
        sec_end_t=datetime.datetime(year_n,mon_n,day_n,hr_n,min_n,sec_n,msec_n)\
                 + datetime.timedelta(seconds=abs(post_ot))
        utc_beg_t=UTCDateTime(sec_beg_t)
        utc_end_t=UTCDateTime(sec_end_t)

        # Template data directory
        ymdhms = ''.join([year_s,mon_s,day_s,hr_s,min_s,sec_is])
        sac_dir = os.path.join(templ_dir_full, ymdhms)

        try:
            os.chdir(sac_dir)
        except OSError:
            print 'WARNING: No such template directory: ', sac_dir
            sys.exit()

        frst_itr = 1
        # Iterate over each station in station list

        for sta_nm in sta_l:

            # Iterate over each component in component list
            os.chdir(sac_dir)
            for comp_nm in comp_l:

                # Station-Component search string
                sta_comp_s = '.'.join(['*', sta_nm, comp_nm, '*'])

                # Remove data from previous iteration
                if 'st' in locals():
                    del st
       
                # Iterate over each template station-component
                os.chdir(sac_dir)
                for sac_file in glob.glob(sta_comp_s):
        
                    # Read in data and save in 'st'
                    if not 'st' in locals():
                        try:
                            st = read(sac_file)
                        except TypeError:
                            pass
                    else:
                        try:
                            st += read(sac_file)
                        except TypeError:
                            pass

                # Process template data, if no errors
                if not 'st' in locals():
                    print 'NOTICE: Template file ', sta_comp_s, \
                          ' missing from directory: ', sac_dir
                elif len(st) == 0:
                    print 'NOTICE: Template file ', sta_comp_s, \
                          ' empty in directory : ', sac_dir
                else:
                    st.sort(['starttime'])
                    st.merge(method=1,fill_value='interpolate', \
                             interpolation_samples=200)
                    st.decimate(dec, no_filter=False, strict_length=False)
                    st.detrend(type = 'linear')
                    st.detrend(type = 'demean')
                    st.filter('bandpass',freqmin=bp_low,freqmax=bp_high, \
                             zerophase=True)
                    st.trim(starttime=utc_beg_t,endtime=utc_end_t,pad=True, \
                            fill_value=0)
                    st.taper(0.05,type='hann',side='both')
                    
                    if frst_itr == 1:
                        mul_str_b = st
                        frst_itr = 0
                    else:
                        mul_str_b += st

        if verbose == 'True':
            print datetime.datetime.now(), 'Correlating template event: ', det

        # Compare templates to data, if no errors
        if not 'st' in locals():
            sys.exit('WARNING: Zero template data files in directories.')
        elif len(st) == 0:
            sys.exit('WARNING: Empty template data files in directories.')
        else:
            trig_b, origin_time = \
                   WvfmSimilarity(st_c,[mul_str_b], abs(pre_ot))
            
            # Combine results above threshold
            raw_lines_cc.extend([(u'{} {} {:.2f}\n'.format( \
                                det, origin_time[k], trig_b[k])) \
                                for k in (xrange(0,len(trig_b),1)) \
                                if trig_b[k] >= cc_thr])

    # Write out raw output file
    if verbose == 'True' and len(raw_lines_cc) != 0:
        with open(out_path_nm, 'a') as fp:
            fp.writelines(raw_lines_cc)

    # Associate raw detections into events, if any
    if len(raw_lines_cc) == 0:
        print 'NOTICE: No detections for data in: ', data_dir_full
    else:
        AssociateDetections(out_dir_full, raw_lines_cc, tf, ymds, cc_thr)


def AssociateDetections(out_dir_a, raw_dets, t_file, ymds_, cc_thr_):
    """
    :param: out_dir_a: location of output directory
    :param: raw_dets: raw detections
    :param: t_file: template information
    :param: ymds_: date string for output file name
    :param: cc_thr_: cross-correlation coefficient threshold value
    """
    
    # Template catalog dictionary with latitude, longitude, depth
    template_dict = {}
    for tline in t_file:
        tfields = tline.split()
        tnum = tfields[0]
        tlat = tfields[3]
        tlon = tfields[4]
        tdpth = tfields[5]
        template_dict[tnum] = [tlat, tlon, tdpth]
    
    # Initialize detection dictionary
    bogus_val = '1900-01-01T00:00:00.000000Z'
    det_dict = {bogus_val: [0, 0, 0, 0, 0]}

    # Iterate over each raw detection
    for dline in raw_dets:
        dfields = dline.split()
        det_num = dfields[0]
        det_ot = dfields[1]
        det_stat = float(dfields[2])

        # Detector Association Time (in seconds)
        DELTA_TIME = 0.5
        
        add_det = 1
        # For detections within deltaT, keep detection with largest DetStat
        for t_key in det_dict.keys():
            
            time_diff = UTCDateTime(t_key) - UTCDateTime(det_ot)
            if abs(time_diff) <= DELTA_TIME and det_dict[t_key][3] < det_stat:
                del det_dict[t_key]
                det_dict[det_ot] = [template_dict[det_num][0], \
                                    template_dict[det_num][1], \
                                    template_dict[det_num][2],det_stat,det_num]
                add_det = 0
            elif abs(time_diff) <= DELTA_TIME and det_dict[t_key][3]>=det_stat:
                add_det = 0

        # Add detection if not within deltaT of existing detection
        if add_det == 1:
            det_dict[det_ot] = [template_dict[det_num][0], \
                                template_dict[det_num][1], \
                                template_dict[det_num][2], det_stat, det_num]

        # Write out results
        assoc_file_nm = 'OutputFile_' + ymds_ + '_' + str(cc_thr_)
        out_path_nm_a = os.path.join(out_dir_a,assoc_file_nm)
        with open(out_path_nm_a,'w') as g:
            for key, val in sorted(det_dict.items()):
                if val[3] >= cc_thr_:
                    key_dt = datetime.datetime.strptime(key, \
                             '%Y-%m-%dT%H:%M:%S.%fZ').strftime( \
                             '%Y/%m/%d %H:%M:%S.%f')
                    g.write( \
                        u'{:23} {:>8} {:>9} {:>6} {:>54} {:>6}:{:.2f}\n'.format( \
                        key_dt[0:23], val[0], val[1], val[2], 2, val[4], val[3]))


if __name__ == '__main__':
    
    argp = argparse.ArgumentParser()
    argp.add_argument('input_file', type=str)
    args = argp.parse_args()
    
    kwargs = {}
    with open(args.input_file,'r') as infile:
        for params in infile:
            a = params.split()
            try:
                kwargs[a[0]] = ' '.join(a[1:])
            except IndexError:
                sys.exit('WARNING: Problem with input file. Check syntax.')

    if len(sys.argv) == 2:
        WSwrapper(**kwargs)
    else:
        print 'WARNING: Incomplete list of arguments in input file. Stopping.'
        sys.exit('usage: python WaveformCC.py InputFile')

