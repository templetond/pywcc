# pywcc
Python tool for seismic waveform cross-correlation

Description
===========

PyWCC version 0.1.1 is a tool to compute seismic waveform cross-correlation on single- 
or multiple-component seismic data across a network of seismic sensors. 
PyWCC correlates waveform data templates with continuous seismic data, associates 
the resulting detections, identifies the template with the highest cross-correlation
coefficient, and outputs a catalog of detections above a user-defined absolute 
cross-correlation threshold value. A note on the correlation: PyWCC computes 
cross-correlations in the frequency-domain for increased computational speed. 


Installation
============

PyWCC runs on Python 2.7 and currently has the following major package dependencies: 
numpy, obspy, python, scipy

These dependencies can be easily installed via Anaconda.


Usage
=====

Once dependencies are installed and PyWCC is downloaded, PyWCC can be run from the 
command line via: 

>> /anaconda/bin/python /PATH/TO/FILE/WaveformCC.py /PATH/TO/FILE/InputFile

or from a Matlab prompt via:

>> [status, out] = system(‘/anaconda/bin/python /PATH/TO/FILE/WaveformCC.py /PATH/TO/FILE/InputFile’)

where a zero value of the status variable indicates the successful completion of 
the command and the out variable contains any status messages printed to screen.


I/O Files
=========

Format of InputFile:
pre_ot     Seconds prior to template zero-time (origin-time) to include in template
post_ot    Seconds following template zero-time (origin-time) to include in template
bp_low     Bandpass filter low corner
bp_high    Bandpass filter high corner
dec        Decimation factor for sample rate. Only integer decimation factor allowed.
cc_thr     Absolute cross-Correlation threshold
sta_l      List of stations to include in analysis, separated by a space
comp_l     List of components to include in analysis, separated by a space
data_dir   Path to continuous seismic data directory. Default=‘./ContinuousData’
out_dir    Path to output directory. Default=‘./OutputFiles’
temp_dir   Template catalog file, inclusive of PATH. Default=‘./TemplateCatalog’
templ_file Path to template data directory. Default=‘./TemplateData’
verbose    Verbose option True/[False]. If True will output raw data file as well
           as final data file and will also print status messages to screen.

Template catalog file has form:
TEMPLATE_NAME YYYY/MM/DD HH:MM:SS.MMM LAT LON DEPTH
YYYY year
MM month
DD day
HH hour
MM minute
SS seconds
MMM decimal seconds
LAT latitude of template event
LON longitude of template event
DEPTH depth of template event 

PyWCC output catalog file has form:
YYYY/MM/DD HH:MM:SS.MMMMMM LAT LON DEPTH CC DET
LAT latitude of the template DET, taken from the Template Catalog
LON longitude of the template DET, taken from the Template Catalog
DEPTH depth of the template DET, taken from the Template Catalog
CC cross-correlation coefficient between the template DET and the new detection
DET template with the highest correlation.

Raw template waveform data should include data before and after the template start
and end times to accommodate possible filter effects. The template waveform data 
should have accurate header information: network, station, location, and channel 
codes. These trace identifiers should be exactly the same as those in the 
continuous data headers. 


Example Test
============

An example test using seismic data from the Bradys Hot Springs geothermal field
is included in the directory ./example. Waveform data, metadata, or data products for 
this study were accessed through the Induced Seismicity Data Website at the Lawrence 
Berkeley National Laboratory which is supported by the U.S. DOE Office of Geothermal Technology.


Acknowledgment
==============

This work was supported by the U.S. Department of Energy (DOE), Office of Energy 
Efficiency and Renewable Energy (EERE), Geothermal Technologies Office (GTO). This 
work was performed under the auspices of the U.S. Department of Energy by Lawrence 
Livermore National Laboratory under Contract DE-AC52-07NA27344.
