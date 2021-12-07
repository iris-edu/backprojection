import os
import sys

from subprocess import Popen, PIPE

from PIL import Image

"""
    Description:

    A Python file that contains all BackProjection data product parameters. You may modify this file to customize 
    plot and animation production. All parameter definitions in this file must follow Python rules. Each 
    parameter group in this file is commented for clarification.

    Copyright and License:

    This software Copyright (c) 2021 IRIS (Incorporated Research Institutions for Seismology).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at
    your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/.

    History:
        2021-11-01 Manoch: v.2021.305 r2 (Python 3) data product release.
"""

# Import the aftershocks libraries.
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# The run settings.
log_to_screen = False
verbose = False

timing = True
timing_threshold = 0.01

# Always set to 1 as this is used for testing ONLY (reduces resolution but speeds up computation for testing).
grid_factor = 1

# Directories.
src_dir = os.path.join(parent_dir, 'src')
image_dir = os.path.join(parent_dir, 'image')
scratch_dir = os.path.join(parent_dir, 'scratch')
video_dir = os.path.join(parent_dir, 'video')
log_dir = os.path.join(parent_dir, 'log')
xml_dir = os.path.join(parent_dir, 'xml')
param_dir = os.path.join(parent_dir, 'param')
lib_dir = os.path.join(parent_dir, 'lib')
assets_dir = os.path.join(parent_dir, 'assets')
metadata_dir = os.path.join(parent_dir, 'metadata')
data_dir = os.path.join(parent_dir, 'data')

ffmpeg = 'ffmpeg'
try:
    process = Popen([ffmpeg, "-version"], stdout=PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()
except Exception as ex:
    print(f"[ERR] ffmpeg application (ffmpeg) not found!\n{ex}")
    sys.exit(2)

# qtfaststart enables streaming and pseudo-streaming of QuickTime and
# MP4 files by moving metadata and offset information to the beginning of the file.
# set to None to disable.
qtfaststart = '/opt/dmc/anaconda/64bit/bin/qtfaststart'
qtfaststart = None

# Time integration resampling, This is to speed up processing. Time integration will be performed at
# the resampling sample interval. Set to None to prevent resampling. The output sample rate will only exactly match
# the selected decimation rate if the original to final rate ratio is factorable by 2,3,5 and 7. Otherwise,
# the closest factorable rate will be chosen.
decimate = 1

# The global trenches file.
# For available linestyles see: https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
global_trenches_file = 'global_trenches_AV.txt'
trench_linestyle = (0, (5, 10))

# Plot and animations.
figure_size = (9.5, 11)
dpi = 150
video_size = (8, 10.8)
video_dpi = 150
frames_per_second = 6
tile_frames_per_second = 8

# Fill color for the beam power plot.
beam_power_fill_color = 'silver'

# Colormap to display local maxima by time.
# https://matplotlib.org/stable/tutorials/colors/colormaps.html
# To use time_cmap, set time_colors to an empty list.
time_colors = ['#2f4f4f', '#8b4513', '#006400', '#bdb76b', '#4b0082', '#ff0000', '#00ced1', '#ffa500', '#ffff00',
               '#00ff00', '#00fa9a', '#0000ff', '#ff00ff', '#6495ed', '#ff1493', '#ffc0cb']
time_color = list()
time_alpha = 0.5
time_cmap = 'jet'

# Logo for the plot.
logo_file = ''
logo_image = os.path.join(assets_dir, logo_file)
if os.path.isfile(logo_image):
    image = Image.open(logo_image)
    logo_width, logo_height = image.size
elif logo_file:
    print(f'[WARN] logo file {logo_image} not found')
    logo_width = 0
    logo_height = 70
else:
    logo_width = 0
    logo_height = 70

"""Code from the question the OffsetImage is given an argument zoom=0.9. 
   This means that each pixel of the original image takes 0.9/0.72=1.25 pixels on the screen. 
   Hence 5 pixels of the original image needs to be squeezed into 4 pixels on the screen. 
   see: https://stackoverflow.com/questions/48639369/does-adding-images-in-pyplot-lowers-their-resolution"""
logo_zoom = 72.0 / dpi
logo_alpha = 1

logo_x = logo_width * logo_zoom * 0.75
logo_y = logo_height * logo_alpha * 0.75

# Padding in pixels between logo and text.
logo_padding = 3

# Logo location as axes pixels
"""'figure points'      points from the lower left corner of the figure
   'figure pixels'      pixels from the lower left corner of the figure
   'figure fraction'    (0, 0) is lower left of figure and (1, 1) is upper right
   'axes points'        points from lower left corner of axes
   'axes pixels'        pixels from lower left corner of axes
   'axes fraction'      (0, 0) is lower left of axes and (1, 1) is upper right
   'data'       use the axes data coordinate system"""
logo_coords = 'axes fraction'
logo_location = (0.01, 0.01)

# Label for plot timestamp.
doi = '10.17611/dp/bp.1'
production_label = f'dp.backprojection doi:{doi}'
production_label_font_size = 8

font_size = dict()
font_size['video'] = {'label': 12, 'legend': 10, 'title': 18, 'time': 16, 'network': 16, 'insufficient_data': 18}
font_size['image'] = {'label': 11, 'legend': 10, 'title': 14, 'time': 16, 'network': 16, 'insufficient_data': 18}

fedcatalog_service_url = 'http://service.iris.edu/irisws/fedcatalog/1/query?'

# Set sta_too_close_km <= 0 to disable sparsifying.
# Virtual or reference networks [NetLat, fNetLon, fdisttooclosekm, statweightdist, cccmin1, cccmin2, statweightazi]
vn_name = 'NA'
virtual_networks = {'GSN': {'lat': -90.0, 'lon': 0.0, 'corners': (-90.0, -180.0, 90.0, 180.0), 'name': 'GSN Stations',
                    'color': 'red', 'marker': '^', 'ccc_min': (0.5, 0.55), 'network': '_GSN', 'xcorr_min': 0.45,
                            'max_sta_count': 300, 'sta_too_close_km': 0.0, 'sta_too_close_deg': 0.0,
                            'std_max': 0.3,
                            'sta_weight_dist': 500.0, 'sta_weight_azim': -1.0, 'sparse_patch_count': 4},

                    'AU': {'lat': -24.3, 'lon': 134.4, 'corners': (-50.0, 90.0, 10.0, 175.0), 'name': 'Australia',
                           'color': 'red', 'marker': '^', 'ccc_min': (0.6, 0.65), 'network': '*', 'xcorr_min': 0.55,
                           'max_sta_count': 50, 'sta_too_close_km': 111.19, 'sta_too_close_deg': 2.0,
                           'std_max': 0.3,
                           'sta_weight_dist': 250.0, 'sta_weight_azim': -1.0, 'sparse_patch_count': 4},

                    'NA': {'lat': 38.9, 'lon': -98.4, 'corners': (0.0, -150.0, 50.0, 0.0), 'name': 'North America',
                           'color': 'red', 'marker': '^', 'ccc_min': (0.6, 0.65), 'network': '*', 'xcorr_min': 0.55,
                           'max_sta_count': 50, 'sta_too_close_km': 111.19, 'sta_too_close_deg': 2.0,
                           'std_max': 0.3,
                           'sta_weight_dist': 250.0, 'sta_weight_azim': -1.0, 'sparse_patch_count': 4},

                    'EU': {'lat': 51.6, 'lon': 20.6, 'corners': (32.0, -10.0, 70.0, 60.0), 'name': 'Europe',
                           'color': 'red', 'marker': '^', 'ccc_min': (0.6, 0.65), 'network': '*', 'xcorr_min': 0.55,
                           'max_sta_count': 50, 'sta_too_close_km': 111.19, 'sta_too_close_deg': 2.0,
                           'std_max': 0.3,
                           'sta_weight_dist': 250.0, 'sta_weight_azim': -1.0, 'sparse_patch_count': 4}
                    }
earthquakes = {'color': 'blue', 'marker': '*'}


# Earth radius, km.
earth_radius = 6378.1

sta_too_close_deg_inc = 0.25
#sta_too_close_deg_inc = {'default': 0.25, 'condition': '>', 'ranges': {500: 0.5, 900: 1.0}}
sta_too_close_deg_init = {'default': 0.25, 'condition': '>', 'ranges': {250: 0.25, 500: 0.5, 750: 0.75, 1000: 1.0}}

vn_check_dist = False
vn_azimuth = 25
vn_min_radius = 0

# For these virtual networks do not perform azimuth check.
vn_azimuth_exception_list = ['GSN']

# Waveform request time window in seconds before and after the event time:
request_time_before = 600.0
request_time_after = {'default': 1300, 'condition': '>', 'ranges': {8.6: 1600.0}}

#channel_order = {'BHZ': 1, 'HHZ': 2}
channel_order = {'BHZ': 1}
request_channel = ','.join(list(channel_order.keys()))

dc_to_exclude = []

# STD QC of waveforms. STD of the trace within std_window seconds before the event time is calculated. If computed
# STD is more than std_max, trace is rejected.
std_check = True
std_window = request_time_before * 0.8
std_offset = 10

# Merge traces with gaps and fill with zero.
merge_gaps = False

# To avoid making one single large request for data to a data center, it is better to make multiple requests.
# The parameter _chunck_count_ in the parameter file determines the number of
# stations per request (chunk) that will be sent to each data center. This number should be adjusted based on the
# number of station-channels involved and the length of each request to avoid memory issues or data center timeouts.
chunk_count = 10

# Filter [filt1,filt2,filt3,filt4]
trace_filter = {'low': {'corners': (0.05, 0.25), 'label': '0.05 to 0.25 Hz'},
                'high': {'corners': (0.25, 1.0), 'label': '0.25 to 1.0Hz'}
                }

# Filter to use for the back projection animation. The "high" filter is used for all except GSN
bp_filter = {'EU': 'high', 'NA': 'high', 'AU': 'high', 'GSN': 'low'}

# fnthroot1; fnthroot2=1./fnthroot1
# Stacking root (integer 1: linear, 3: cube root stacking) [root, fnthroot1=root, fnthroot2=1./fnthroot1].
stacking_root = 2

# Beam averaging window length (seconds).
beam_average_seconds = 5

# ixcorryesno
do_xcorr = False

# Use absolute values for MCCC?
xcorr_abs = False

# The minimum number of stations that we must have.
min_num_sta = 1

# Length (sec) of synthetic trapezoid function to replace data.
# 0 to use data.  Less than 1.0 will give you a 1-sec wide triangle.
trapezoid_length = 0

# The taper fraction of the cosine taper is applied to the waveform data in the time domain before deconvolution.
taper_fraction = 0.05

"""" Output type after deconvolution.
DISP"
displacement, output unit is meters
"VEL"
velocity, output unit is meters/second
"ACC"
acceleration, output unit is meters/second**2"""
output_type = 'VEL'

# A bandpass filter to apply in frequency domain to the data before deconvolution.
pre_filter = [0.018, 0.02, 2, 2.2]
prefilter_label = f'pre-filter {pre_filter[1]} to {pre_filter[2]} Hz'

# Peak amplitude marker on videos. For available markers see:
# https://matplotlib.org/stable/api/markers_api.html
# peak_marker_lw indicates that peak_marker will be plotted when amplitude drops below this value.
peak_marker = 'P'
peak_marker_size = 40
peak_marker_color = 'maroon'
peak_marker_lw = 2
peak_marker_max = 0.05

# Maximum allowable distance between max peak and the event location (km). Set to a large number to deactivate.
peak_offset_max_km = 1000

# Phase to use (P, PP, PKIKP, S) #738.
phase = 'P'
pre_phase_seconds = 5.0
post_phase_seconds = 120.0

stf_pre_phase_seconds = 10.0
stf_post_phase_seconds = 40.0
stf_search_seconds = 45

# Cross-correlation shift in seconds.
xcorr_window = [5.0, 15.0, 22.0, 30.0, 40.0]
xcorr_shift_default = 5.0

# Trace sampling.
trace_sampling_frequency = 20.0
trace_sampling_interval = 1.0 / trace_sampling_frequency

# Min and max distance from earthquake [delmin,delmax].
eq_min_radius = 30
eq_max_radius = 97.0
vn_max_radius = 120.0

# Distance circles to draw (degrees).
distance_circles = (30.0, 60.0, 95.0)
distance_circle_labels = [f'{int(d)}°' for d in distance_circles]
map_width = 28000000 * 0.8

# Basemap distance scale length in km.
scalebar_km = 100

# Find the local maxima peak and plot these within peak_marker_base_factor (0.5 = 50%) of the global maximum.
peak_marker_base_factor = 0.5

qc_max_peak_factor = 0.7
qc_max_peak_count = 5
qc_vn_max_distance = 25

# Maximum peak amplitude marker size.
peak_marker_size_max = 600
peak_marker_size = 60

# Azi_min Azimuth_max (deg going clockwise 0-360) [fazimin,fazimax].
azimuth_min = 0.0
azimuth_max = 360.0

# SNR window.
snr_window = {'pre-p': [-50, -10], 'p': [-5, 35]}

# Pre-P / P SNR on the unfiltered trace.
min_snr = 2.0

# Pre-P / P SNR on the filtered trace (set to None to disable).
min_filtered_snr = 2.0

# Should use filtered trace tp align MCCC [ifiltalign=0].
filtered_trace_align = False

# Seconds before and after P-arrivals trace must exist.
seconds_before_p = 45.0
seconds_after_p = 30.0

""" All Variable settings must be given in a dictionary structure (see lib.case):
        default value
        condition to impose >, >=, ==, <, <=, !=
        range: value pairs
"""
# STF trace length settings.
stf_length_dict = {'default': 60.0, 'condition': '>', 'ranges': {7.0: 100.0, 7.5: 120.0, 8.0: 170.0}}

# Cross-correlation window length based on magnitude.
mag_window_length_dict = {'default': 15.0, 'condition': '>', 'ranges': {7.0: 17.5, 7.9: 20.0}}

# BP initial time in seconds [initialtime].
bp_initial_seconds = 0.0

# Seconds before the event time when BP starts.
bp_time_offset = 30.0

# BP time increment [fitincrement].
bp_time_inc_dict = {'default': 0.25, 'condition': '>', 'ranges': {8.49: 0.5}}

# BP total time, time increment, and time offset [ittotal].
bp_time_total_dict = {'default': 90.0, 'condition': '>', 'ranges': {6.99: 130.0, 7.49: 180.0, 7.99: 210.0, 8.49: 280.0,
                                                                    8.99: 630.0, 9.09: 630.0, 9.19: 630.0}}

# The total number of seconds of the Hanning averaging taper [tavg]
hann_t_avg_dict = {'default': 10.0, 'condition': '>', 'ranges': {6.99: 10.0, 7.99: 15.0, 8.49: 20.0, 8.99: 30.0}}

# Latitude grid increment [glatinc, gLatHalf]. Grid latitude starts gLatHalf before the eq latitude'
grid_decimal_places = 2
grid_latitude_inc_dict = {'default': 0.1, 'condition': '>', 'ranges': {8.51: 0.12, 8.99: 0.13, 9.09: 0.15}}
grid_latitude_half_dict = {'default': 2.0, 'condition': '>', 'ranges': {6.99: 2.5, 7.49: 3.0, 7.99: 3.5, 8.49: 4.0,
                                                                        8.99: 5.0, 9.09: 6.0, 9.19: 10.0}}
# Synthetic triangle width in seconds.
triangle_width = 5.0

# Travel time model.
travel_time_model = 'iasp91'

# Create a summary plot.
create_summary = True

# Minimum/maximum trace time (seconds) for the summary plots
max_summary_trace_time = 90
min_summary_trace_time = -25

# Y position of the summary plot trace labels (max 1).
label_y_position = 0.70

# What to create.
create_animation = True

# Arrow characters.
arrows = (u'$\u2191$', u'$\u2193$')

pcolormesh_grid_factor = 5

# Resolution of boundary database to use. Can be c (crude), l (low), i (intermediate), h (high), f (full) or None.
# If None, no boundary data will be read in. Higher-res datasets are much slower to draw. None will be the fastest
# with no boundaries drawn.
basemap_resolution = 'h'
basemap_countries = False
basemap_states = False

# Azimuthal Equidistant Projection, shortest route from the center of the map to any other point is a straight line.
basemap_projection = 'aeqd'

# Options for the basemap's continent colors.
fill_continents = True
fill_continents_color = '#D3D3D3'
fill_continents_alpha = 0.2
sta_map_color = 'gray'

bp_colors = [(1.0, 1.0, 1.0),
             (0.95, 1.0, 1.0),
             (0.93, 1.0, 1.0),
             (0.9, 1., 1.),
             (0.85, 1.0000, 1.0000),
             (0.8, 1.0000, 1.0000),
             (0.7, 1.0000, 1.0000),
             (0.6, 1.0000, 1.0000),
             (0.5, 1.0000, 1.0000),
             (0.4, 1.0000, 1.0000),
             (0.3, 1.0000, 1.0000),
             (0.2, 1.0000, 1.0000),
             (0.1, 1.0000, 1.0000),
             (0, 1.0000, 1.0000),
             (0.0769, 1.0000, 0.9231),
             (0.1538, 1.0000, 0.8462),
             (0.2308, 1.0000, 0.7692),
             (0.3077, 1.0000, 0.6923),
             (0.3846, 1.0000, 0.6154),
             (0.4615, 1.0000, 0.5385),
             (0.5385, 1.0000, 0.4615),
             (0.6154, 1.0000, 0.3846),
             (0.6923, 1.0000, 0.3077),
             (0.7692, 1.0000, 0.2308),
             (0.8462, 1.0000, 0.1538),
             (0.9231, 1.0000, 0.0769),
             (1.0000, 1.0000, 0),
             (1.0000, 0.9231, 0),
             (1.0000, 0.8462, 0),
             (1.0000, 0.7692, 0),
             (1.0000, 0.6923, 0),
             (1.0000, 0.6154, 0),
             (1.0000, 0.5385, 0),
             (1.0000, 0.4615, 0),
             (1.0000, 0.3846, 0),
             (1.0000, 0.3077, 0),
             (1.0000, 0.2308, 0),
             (1.0000, 0.1538, 0),
             (1.0000, 0.0769, 0),
             (1.0000, 0, 0),
             (0.9231, 0, 0),
             (0.8462, 0, 0),
             (0.7692, 0, 0)]


