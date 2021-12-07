#!/usr/bin/env python
import sys
import os
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.transforms import Bbox

import warnings

import shutil

from obspy.clients.fdsn import Client
from obspy.geodetics import degrees2kilometers
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy import UTCDateTime
from obspy.core.stream import Stream

from datetime import datetime

from time import time
import math

import subprocess

import getopt

import glob
import numpy as np

from PIL import Image

import matplotlib
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
matplotlib.use('Agg')

from scipy.interpolate import griddata
from scipy.signal import argrelextrema
import pickle

# Import the back projection parameters and libraries.
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
lib_dir = os.path.join(parent_dir, 'lib')
param_dir = os.path.join(parent_dir, 'param')

sys.path.append(param_dir)
sys.path.append(lib_dir)
import back_projection_param as param

import back_projection_lib as lib
from back_projection_lib import print_message, ObjDict


# Ignore user warnings.
warnings.filterwarnings("ignore")

"""
    Description:

    This is the Python 3 code behind the IRIS DMC's BackProjection Data product (http://ds.iris.edu/spud/backprojection)
    and it can producing the individual plots and animations that are part of the BackProjection  product 
    (http://ds.iris.edu/spud/backprojection).

    The code can be configured via its parameter file "back_projection_param.py" and via the command line arguments. 
    Currently parameters are optimized for use with four virtual networks defined in the parameter file. 

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
        2011-10-01 Alex Hutko: r1, development and initial product release (Fortran).

"""

# Script info.
script_version = 'v.2021.305'
script = sys.argv[0]
script = os.path.basename(script)

# Parameters.
up_arrow, down_arrow = (u'$\u2191$', u'$\u2193$')

fedcatalog_service_url = param.fedcatalog_service_url

logo_image = param.logo_image
logo_zoom = param.logo_zoom
logo_alpha = param.logo_alpha
logo_coords = param.logo_coords
logo_location = param.logo_location
logo_width = param.logo_width
logo_height = param.logo_height
logo_padding = param.logo_padding

basemap_countries = param.basemap_countries
basemap_states = param.basemap_states
basemap_projection = param.basemap_projection

pcolormesh_grid_factor = param.pcolormesh_grid_factor

production_label_font_size = param.production_label_font_size

bp_filter = None

vn = param.virtual_networks
vn_list = list(vn.keys())
vn_min_radius = param.vn_min_radius
vn_max_radius = param.vn_max_radius
vn_azimuth = param.vn_azimuth

eq_min_radius = param.eq_min_radius
eq_max_radius = param.eq_max_radius

dc_to_exclude = param.dc_to_exclude

grid_factor = param.grid_factor

distance_circles = param.distance_circles
distance_circle_labels = param.distance_circle_labels

timing = param.timing

decimate = param.decimate
pre_phase_seconds = param.pre_phase_seconds
post_phase_seconds = param.pre_phase_seconds

trace_sampling_frequency = param.trace_sampling_frequency
trace_sampling_interval = param.trace_sampling_interval

chunk_count = param.chunk_count

earthquakes = param.earthquakes

xcorr_shift = param.xcorr_shift_default

create_synthetics = param.create_animation
create_summary = param.create_summary
create_animation = param.create_animation

image_dir = lib.mkdir(param.image_dir)
video_dir = lib.mkdir(param.video_dir)
scratch_dir = lib.mkdir(param.scratch_dir)
log_dir = lib.mkdir(param.log_dir)
metadata_dir = lib.mkdir(param.metadata_dir)
data_dir = lib.mkdir(param.data_dir)

font_size = param.font_size

log_to_screen = param.log_to_screen

# Travel time cache to speed things up.
tt_cache = dict()

# Distance between pair of stations.
intra_station_dist = dict()

# Station inventory list.
inventory_list = dict()

coastline_skipped = False
verbose = param.verbose


def usage():
    """The usage message.
    """
    new_line = '\n'
    print(f'{new_line}{new_line}{script} ({script_version}):', file=sys.stdout, flush=True)

    print(f'{new_line}{new_line}This is the Python 3 code behind the IRIS DMC\'s BackProjection data product (BP):\n'
          f'http://ds.iris.edu/ds/products/backprojection/\n'
          f'\nand can producing the individual plots and animations that are part of the PackProjection product:\n'
          f'http://ds.iris.edu/spud/backprojection{new_line}', file=sys.stdout, flush=True)

    print(f'The code can be configured via its parameter file "back_projection_param.py" or via the '
          f'command line arguments. {new_line}{new_line}'
          f'Currently parameters are optimized for use with four preconfigured virtual networks:',
          file=sys.stdout, flush=True)

    print(f'\t\tvirtual network\t\tname', file=sys.stdout, flush=True)
    print(f'\t\t{15 * "="}\t\t{11 * "="}', file=sys.stdout, flush=True)
    _vn = param.virtual_networks
    _ind = sorted(_vn.keys())
    for _i in _ind:
        print(f'\t\t{_i}\t\t\t{_vn[_i]["name"]}', file=sys.stdout, flush=True)

    print(f'{new_line}Virtual networks could be modified by changing the '
          f'"virtual_networks" parameter in the parameter file.'
          f'{new_line}{new_line}command line options:{new_line}'
          f'\t-h --help\t\tthis message{new_line}'
          f'\t-v --verbose \t\t[default: {verbose}] turns the verbose mode on{new_line}'
          f'\t-a --anim\t\t[default: {create_animation}] create animations [True/False]{new_line}'
          f'\t-l --logscreen\t\t[default: {param.log_to_screen}] send the log messages to screen{new_line}'
          f'\t-n --vnet\t\t[required] virtual network code (see above) {new_line}'
          f'\t-m --emag\t\t[required] event magnitude {new_line}'
          f'\t-s --summary\t\t[default: {create_summary}] create summary plot [True/False]{new_line}'
          f'\t-t --etime\t\t[required] the event time as "YYYY-MM-DDTHH:MM:SS" {new_line}'
          f'\t-x --elon\t\t[required] the event longitude {new_line}'
          f'\t-y --elat\t\t[required] the event latitude {new_line}'
          f'\t-z --edepth\t\t[required] the event depth (km) {new_line}'
          f'\t-d --decimate\t\t[default: {decimate}] the desired animations sample rate in seconds'
          f'{new_line}\t\t\t\tThe output sample rate will only exactly match the selected'
          f'{new_line}\t\t\t\tdecimation rate if the ratio of original to final rate is a whole number'
          f'{new_line}\t\t\t\tOtherwise, the closest factorable rate will be chosen. {new_line}'
          f'\t-f --factor\t\t[default: {grid_factor}] this parameter could be used for testing. '
          f'The grid spacing is {new_line}'
          f'\t\t\t\tmultiplied by this factor (-f 5 or -f 10 are reasonable choices) {new_line}'
          f'\t\t\t\tand as a result the resolution is reduced and computation take place faster.{new_line}'
          f'\t\t\t\tYou could use this option to '
          f'test your parameters before a production run.{new_line}{new_line}'
          f'Examples:{new_line}'
          f'\tlower resolution (faster, good for tuning parameters):{new_line}'
          f'\t\tpython src/{script} -m 7.8 -y 55.030 -x -158.522 -z 28 -n AU -t 2020-07-22T06:12:44 -l -f 5'
          f'{new_line}{new_line}\thigher resolution (slower, for final product):{new_line}'
          f'\t\tpython src/{script} -m 7.8 -y 55.030 -x -158.522 -z 28 -n AU -t 2020-07-22T06:12:44 -l'
          , file=sys.stdout, flush=True)
    print('\n\n', file=sys.stdout, flush=True)


def full_extent(ax, pad=0.0):
    """Get the full extent of an axes, including axes labels, tick labels, and
    titles. FroM;
    https://stackoverflow.com/questions/4325733/save-a-subplot-in-matplotlib/4328608"""
    # For text objects, we need to draw the figure first, otherwise the extents
    # are undefined.
    ax.figure.canvas.draw()
    items = ax.get_xticklabels() + ax.get_yticklabels()
    items += [ax, ax.title]
    bbox = Bbox.union([item.get_window_extent() for item in items])
    return bbox.expanded(1.0 + pad, 1.0 + pad)


def plot_stack(_stack):
    """This a debug function to plot a list"""
    for _grid_key in _stack:
            stack_list = list()
            for _t_key in _stack[_grid_key]:
                stack_list.append(_stack[_grid_key][_t_key])
            plt.plot(stack_list)
            plt.show()
            return


def make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, search_radius, anim_tag=''):
    """Create an insufficient data image."""

    # Since Basemap takes a while to draw, save it as a pickle file and reload it in the loop.
    # File for use by Python's object serialization module pickle. Added time tag to make it unique to this run.
    pickle_file = os.path.join(param.scratch_dir, f'{vn_name}_{int(datetime.now().timestamp())}.pickle')
    save_pickle = True

    # What type of media is generated?
    media = 'video'

    coastline_skipped = False
    title = f"{eq_datetime.strftime('%Y-%m-%d %H:%M:%S')} M{eq_magnitude} Z={eq_depth}km\n" \
            f"{vn_name} VNet, {param.trace_filter[param.bp_filter[vn_name]]['label']}"
    if anim_tag == 'syn':
        anim_tag = f'BP_{anim_tag}'
        title = f'{title} (ARF, synthetics)'
    else:
        anim_tag = f'BP'

    t7 = time()

    # Get grid locations.
    latitude, longitude = lib.set_grid(eq_lat, eq_lon, eq_magnitude)

    lon_0 = eq_lon
    lat_0 = eq_lat

    # Video frame layout.
    subplot_columns = 1
    subplot_tall_rows = 1
    subplot_short_rows = 1
    tall_to_short_height = 3

    file_tag = '_'.join([anim_tag, vn_name, lib.file_name_tag(eq_datetime)])
    screen_file_tag = '_'.join([anim_tag, 'screen', vn_name, lib.file_name_tag(eq_datetime)])

    # Remove image files if they exist from previous runs.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, file_tag)}*.png')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print_message('ERR', f'Failed to remove\n {_er}', flush=True, log=log_file)

    delimiter = ' '
    # Production date time stamp.
    production_date = lib.version_timestamp(script_version, search_radius, delimiter=delimiter)

    # The moving dot should show color of amplitude.
    normalize = colors.Normalize(vmin=0, vmax=1)

    fig = plt.figure(figsize=param.video_size, facecolor='white')
    ax = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                          (0, 0), rowspan=tall_to_short_height, colspan=1)

    lat_min =latitude['start']
    lat_max = latitude['end']
    lon_min = longitude['start']
    lon_max = longitude['end']

    # Create the basemap
    width = lat_max - lat_min
    width = degrees2kilometers(width) * 1000.0
    print_message('INFO', f'Basemap(width={width}, height={width}, projection={basemap_projection}, lat_0={lat_0}, '
                          f'lon_0={lon_0}, resolution={param.basemap_resolution})', flush=True, log=log_file)
    bm = Basemap(width=width, height=width, projection=basemap_projection,
                 lat_0=lat_0, lon_0=lon_0, resolution=param.basemap_resolution)

    # Avoid areas without coastlines.
    try:
        bm.drawcoastlines(color=param.fill_continents_color)
    except Exception as ex:
        if not coastline_skipped:
            print_message('WARN', f'Skipped drawcoastlines:\n{ex}', flush=True, log=log_file)
            coastline_skipped = True
        pass

    if basemap_countries:
        bm.drawcountries(color=param.fill_continents_color)
        if basemap_states:
            bm.drawstates(color=param.fill_continents_color)

    # labels = [left,right,top,bottom].
    bm.drawparallels(np.arange(int(lat_min), int(lat_max), 1), labels=[1, 0, 0, 0],
                     fontsize=font_size[media]['label'], linewidth=0.0,
                     labelstyle='+/-', fmt='%0.0f')
    bm.drawmeridians(np.arange(int(lon_min), int(lon_max), 2), labels=[0, 0, 0, 1],
                     rotation=0, fontsize=font_size[media]['label'],
                     linewidth=0.0,
                     labelstyle='+/-', fmt='%0.0f')

    trench_x, trench_y = lib.read_global_trenches(bmap=bm)

    # Earthquake location as map units.
    xpt, ypt = bm(eq_lon, eq_lat)

    # Mark the earthquake location.
    bm.plot([xpt], [ypt], marker=earthquakes['marker'], markerfacecolor='black',
            markeredgecolor=earthquakes['color'],
            markersize=15, label='event')

    # plt.ylabel('Latitude', labelpad=param.ylabel_pad, fontsize=font_size[media]['label'])
    # plt.xlabel('Longitude', labelpad=param.xlabel_pad, fontsize=font_size[media]['label'])

    bm.plot(trench_x, trench_y, color='black', linestyle=param.trench_linestyle, linewidth=0.5)

    # Insufficient Data in the middle.
    plt.text(0.5, 0.8, f'Insufficient Data', fontsize=font_size[media]['insufficient_data'],
             horizontalalignment='center',
             verticalalignment='center',
             backgroundcolor='white', color='maroon', weight='bold', transform=ax.transAxes)

    plt.title(title, fontsize=font_size[media]['title'])

    # Get info on this subplot so we can align the one below it.
    map_bbox = ax.get_position()

    # Plot the beam power.
    ax0 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                           (3, 0), rowspan=1, colspan=1)
    ax0.axes.yaxis.set_ticklabels([])

    # Get info on this subplot so we can align the one above it.
    # We always want to adopt the map's width since it is dynamic.
    power_bbox = ax0.get_position()
    ax0.set_position([map_bbox.x0, power_bbox.y0, map_bbox.width, power_bbox.height])

    # Plot the logo.
    logo_x0 = map_bbox.x0 * param.video_dpi
    if os.path.isfile(param.logo_image):
        image = np.array(Image.open(param.logo_image))
        im = OffsetImage(image, zoom=logo_zoom, alpha=logo_alpha, zorder=1000)
        b_box = AnnotationBbox(im, xy=logo_location, xycoords=logo_coords, box_alignment=(0.0, 0.0),
                               boxcoords='offset pixels', frameon=False)
        ax0.add_artist(b_box)

    xytext = (logo_padding,  - 1.5 * logo_height)

    plt.annotate(production_date, xy=xytext, xycoords='axes pixels',
                 xytext=(0, 0),
                 textcoords='offset pixels', horizontalalignment='left',
                 fontsize=production_label_font_size,
                 verticalalignment='center')
    plt.savefig(os.path.join(scratch_dir, f'{file_tag}_'
                                          f'{0:06d}.png'), bbox_inches='tight', pad_inches=0.25,
                dpi=param.video_dpi, facecolor='white')
    plt.close()

    print_message('INFO', f'Creating the video:', flush=True, log=log_file)
    # Apply -vf pad=ceil(iw/2)*2:ceil(ih/2)*2 filter to avoid eight not divisible by 2 (1644x1491)
    # error without rescaling.
    # The -r before the input means the video will play at that number of the original images per second.
    # Have to define -r for both input and output to avoid dropped frames due to different default input and
    # and the desired output rate.
    command = f'{param.ffmpeg} -loop 1 -r {param.frames_per_second}  ' \
              f'-i {os.path.join(scratch_dir, file_tag)}_%06d.png ' \
              f'-c:v libx264 -pix_fmt yuv420p -crf 23 -t {bp_t_total} ' \
              f'-r {param.frames_per_second} ' \
              f'-vf pad=ceil(iw/2)*2:ceil(ih/2)*2 ' \
              f'-y {os.path.join(video_dir, file_tag)}.mp4'.split()
    print_message('INFO', f'Creating the video: {command}', flush=True, log=log_file)
    subprocess.call(command)

    if param.qtfaststart is not None:
        command = f'{param.qtfaststart} {os.path.join(video_dir, file_tag)}.mp4 ' \
                  f'{os.path.join(video_dir, file_tag)}_q.mp4'.split()
        print_message('INFO', f'QTfaststart the video: {command}', flush=True, log=log_file)
        subprocess.call(command)

    # Remove image files if they exist from previous runs.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, file_tag)}*.png')
        for file_index, this_file in enumerate(files_to_remove):
            if file_index == 0:
                shutil.move(this_file, f'{os.path.join(video_dir, screen_file_tag)}.png')
            else:
                os.remove(this_file)
    except Exception as _er:
        print_message('ERR', f'Failed to remove\n {_er}', flush=True, log=log_file)

    return


def make_animation(anim_stack_start, anim_stack_end, anim_stack_amp, anim_stack_amp_loc, anim_stack_max,
                   anim_stack_max_loc, anim_global_max, search_radius, anim_tag='', grid_factor=1):
    """Create animation from the stacked traces."""

    # What type of media is generated?
    media = 'video'

    image_count = 0

    coastline_skipped = False

    # Since Basemap takes a while to draw, save it as a pickle file and reload it in the loop.
    # File for use by Python's object serialization module pickle. Added time tag to make it unique to this run.
    pickle_file = os.path.join(param.scratch_dir, f'{vn_name}_{int(datetime.now().timestamp())}.pickle')
    save_pickle = True

    title = f"{eq_datetime.strftime('%Y-%m-%d %H:%M:%S')} M{eq_magnitude} Z={eq_depth}km\n" \
            f"{vn_name} VNet, {param.trace_filter[param.bp_filter[vn_name]]['label']}"
    if anim_tag == 'syn':
        anim_tag = f'BP_{anim_tag}'
        title = f'{title} (ARF, synthetics)'
    else:
        anim_tag = f'BP'

    t7 = time()

    # Get grid locations.
    latitude, longitude = lib.set_grid(eq_lat, eq_lon, eq_magnitude, grid_factor=grid_factor)

    lon_0 = eq_lon
    lat_0 = eq_lat

    # Make the custom color map.
    bp_cmap = lib.make_cmap(param.bp_colors, bit=False, log=log_file)

    # Video frame layout.
    subplot_columns = 1
    subplot_tall_rows = 1
    subplot_short_rows = 1
    tall_to_short_height = 3

    file_tag = '_'.join([anim_tag, vn_name, lib.file_name_tag(eq_datetime)])

    # Remove image files if they exist from previous runs.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, file_tag)}*.png')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print_message('ERR', f'Failed to remove\n {_er}', flush=True, log=log_file)

    time_key_list = list(anim_stack_amp.keys())

    delimiter = ' '
    # Production date time stamp.
    production_date = lib.version_timestamp(script_version, search_radius, delimiter=delimiter)

    # The moving dot should show color of amplitude.
    normalize = colors.Normalize(vmin=0, vmax=1)
    max_peak = -1
    for time_key in time_key_list:
        # Limit the animation to the actual BP time limits.
        if float(time_key) < anim_stack_start or float(time_key) > anim_stack_end:
            continue

        fig = plt.figure(figsize=param.video_size, facecolor='white')
        ax = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                              (0, 0), rowspan=tall_to_short_height, colspan=1)

        lat = list()
        lon = list()
        val = list()
        for _index, _val in enumerate(anim_stack_amp[time_key]):
            lat_key, lon_key = anim_stack_amp_loc[time_key][_index]
            lat.append(float(lat_key))
            lon.append(float(lon_key))
            val.append(_val / anim_global_max)

        # Must be np arrays for grid
        lat_list = np.array(lat)
        lon_list = np.array(lon)
        value_list = np.array(val)

        # Find the min and max of coordinates.
        lon_min = lon_list.min()
        lat_min = lat_list.min()
        lon_max = lon_list.max()
        lat_max = lat_list.max()

        # Create the basemap and save it to reuse
        if save_pickle:
            width = lat_max - lat_min
            width = degrees2kilometers(width) * 1000.0
            bm = Basemap(width=width, height=width, projection=basemap_projection,
                         lat_0=lat_0, lon_0=lon_0, resolution=param.basemap_resolution)

            # Avoid areas without coastlines.
            try:
                bm.drawcoastlines(color=param.fill_continents_color)
            except Exception as ex:
                if not coastline_skipped:
                    print_message('WARN', f'Skipped drawcoastlines:\n{ex}', flush=True, log=log_file)
                    coastline_skipped = True
                pass

            if basemap_countries:
                bm.drawcountries(color=param.fill_continents_color)
                if basemap_states:
                    bm.drawstates(color=param.fill_continents_color)

            pickle.dump(bm, open(pickle_file, 'wb'), -1)
            save_pickle = False

        else:
            # Read pickle back in and plot it again (should be much faster).
            bm = pickle.load(open(pickle_file, 'rb'))

        # labels = [left,right,top,bottom].
        bm.drawparallels(np.arange(int(lat_min), int(lat_max), 1), labels=[1, 0, 0, 0],
                         fontsize=font_size[media]['label'],
                         linewidth=0.0,
                         labelstyle='+/-', fmt='%0.0f')
        bm.drawmeridians(np.arange(int(lon_min), int(lon_max), 2), labels=[0, 0, 0, 1], rotation=0,
                         fontsize=font_size[media]['label'],
                         linewidth=0.0,
                         labelstyle='+/-', fmt='%0.0f')

        # Avoid areas without coastlines.
        try:
            bm.drawcoastlines(color=param.fill_continents_color)
        except Exception as ex:
            if not coastline_skipped:
                print_message('WARN', f'Skipped drawcoastlines:\n{ex}', flush=True, log=log_file)
                coastline_skipped = True
            pass

        if basemap_countries:
            bm.drawcountries(color=param.fill_continents_color)
            if basemap_states:
                bm.drawstates(color=param.fill_continents_color)

        # Read global trenches coordinates.
        trench_x, trench_y = lib.read_global_trenches(bmap=bm)

        # Earthquake location in map units.
        xpt, ypt = bm(eq_lon, eq_lat)

        # Mark the earthquake location.
        bm.plot([xpt], [ypt], marker=earthquakes['marker'], markerfacecolor=earthquakes['color'],
                markeredgecolor='white',
                markersize=15, label='event')

        # plt.ylabel('Latitude', labelpad=param.ylabel_pad, fontsize=font_size[media]['label'])
        # plt.xlabel('Longitude', labelpad=param.xlabel_pad, fontsize=font_size[media]['label'])

        # Plot a vertical distance scale.
        scale_deg = kilometer2degrees(param.scalebar_km)
        lat0 = lat_min + (lat_max - lat_min) / 4.0
        lon0 = lon_max - (lon_max - lon_min) / 10.0
        lat1, lon1 = lib.get_location(lat0, lon0, 0, scale_deg)
        xs, ys = bm((lon0, lon0), (lat0, lat1))

        # Use the same _X to ensure scale is vertical.
        bm.plot([xs[0], xs[0]], ys, color='black', linewidth=2)
        plt.text(xs[0], ys[1], f'  {param.scalebar_km} km', fontsize=font_size[media]['legend'],
                 horizontalalignment='center', rotation=90,
                 verticalalignment='bottom',
                 path_effects=[pe.withStroke(linewidth=2, foreground='white')])

        # Now let's grid the data. Find the number of grid points in each direction.
        lon_num = pcolormesh_grid_factor * int(((lon_max - lon_min) / longitude['inc']) + 1)
        lat_num = pcolormesh_grid_factor * int(((lat_max - lat_min) / latitude['inc']) + 1)

        # Create a uniform mesh for contouring. First transfer lon, lat to map units (x, y).
        x_old, y_old = bm(lon_list, lat_list)
        x_new = np.linspace(min(x_old), max(x_old), lon_num)
        y_new = np.linspace(min(y_old), max(y_old), lat_num)

        # Basic mesh in map's x, y.
        grid_x, grid_y = np.meshgrid(x_new, y_new)

        try:
            # Interpolate at the new grid points.
            # Method : {'linear', 'nearest', 'cubic'}, optional.
            grid_v = griddata((x_old, y_old), value_list, (grid_x, grid_y), method='cubic')
        except Exception as _er:
            print_message('ERR', f'Griding failed: {_er}', log=log_file)
            sys.exit(1)

        # Create a pseudocolor plot.
        bm.pcolormesh(grid_x, grid_y, grid_v, cmap=bp_cmap, alpha=0.7, shading='auto', linewidths=0,
                      vmin=0, vmax=1)

        # Plot a + at the peak location.
        if anim_stack_max[time_key] / anim_global_max < param.peak_marker_max:
            _x, _y = bm(anim_stack_max_loc[time_key][1], anim_stack_max_loc[time_key][0])
            bm.scatter(_x, _y, s=param.peak_marker_size, marker=param.peak_marker, facecolor=param.peak_marker_color,
                       edgecolors=param.peak_marker_color, linewidths=param.peak_marker_lw,
                       linestyle='None', zorder=1000)

        # Plot the trench.
        bm.plot(trench_x, trench_y, color='black', linestyle=param.trench_linestyle, linewidth=0.5)

        # Frame time on the upper right
        plt.text(0.95, 0.95, f'{float(time_key):0.1f} sec ({anim_stack_max_loc[time_key][0]:0.2f}, '
                             f'{anim_stack_max_loc[time_key][1]:0.2f})', fontsize=font_size[media]['time'],
                 horizontalalignment='right', verticalalignment='center',
                 backgroundcolor='white', transform=ax.transAxes)

        plt.title(title, fontsize=font_size[media]['title'])

        # Get info on this subplot so we can align the one below it.
        map_bbox = ax.get_position()

        # Plot the beam power.
        ax0 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                               (3, 0), rowspan=1, colspan=1)
        times = list(anim_stack_max.keys())
        times = list(np.array(times, dtype=float))
        values = list(anim_stack_max.values())
        max_value = max(values)
        values = np.array(values) / max_value
        values = list(values)
        ax0.fill_between(times, 0, values, facecolor=param.beam_power_fill_color)
        ax0.set_xlim(stack_start, stack_end)
        ax0.set_ylim(bottom=0.0)
        ax0.set_xlabel('time relative to origin (sec.)')
        ax0.set_ylabel('beam power')
        ax0.axes.yaxis.set_ticklabels([])

        # Get info on this subplot so we can align the one above it.
        # We always want to adopt the map's width since it is dynamic.
        power_bbox = ax0.get_position()
        ax0.set_position([map_bbox.x0, power_bbox.y0, map_bbox.width, power_bbox.height])

        # Plot the logo.
        logo_x0 = map_bbox.x0 * param.video_dpi
        if os.path.isfile(param.logo_image):
            image = np.array(Image.open(param.logo_image))
            im = OffsetImage(image, zoom=logo_zoom, alpha=logo_alpha, zorder=1000)
            b_box = AnnotationBbox(im, xy=logo_location, xycoords=logo_coords, box_alignment=(0.0, 0.0),
                                   boxcoords='offset pixels', frameon=False)
            ax0.add_artist(b_box)

        xytext = (logo_padding,  - 1.5 * logo_height)

        plt.vlines(float(time_key), 0.0, 1.0)
        v = float(anim_stack_max[time_key]) / max_value
        if v > max_peak:
            max_time_frame_file = f'{os.path.join(scratch_dir, file_tag)}_{image_count + 1:06d}.png'
            max_peak = v
        plt.scatter([float(time_key)], [v], cmap=bp_cmap,
                    c=[v], s=40, marker='o', edgecolors='k', linewidths=0.3, norm=normalize,
                    linestyle='None', zorder=1000)

        image_count += 1

        plt.annotate(production_date, xy=xytext, xycoords='axes pixels',
                     xytext=(0, 0),
                     textcoords='offset pixels', horizontalalignment='left',
                     fontsize=production_label_font_size,
                     verticalalignment='center')
        plt.savefig(os.path.join(scratch_dir, f'{file_tag}_'
                                              f'{image_count:06d}.png'), bbox_inches='tight', pad_inches=0.25,
                    dpi=param.video_dpi, facecolor='white')
        if timing:
            t7 = lib.time_it(t7, f'Image for time {time_key}', log=log_file)

        plt.close()

    # create the video.
    try:
        print_message('INFO', f'Creating the video:', flush=True, log=log_file)
        # Apply -vf pad=ceil(iw/2)*2:ceil(ih/2)*2 filter to avoid eight not divisible by 2 (1644x1491)
        # error without rescaling.
        # The -r before the input means the video will play at that number of the original images per second.
        command = f'{param.ffmpeg} -r {param.frames_per_second}  ' \
                  f'-i {os.path.join(scratch_dir, file_tag)}_%06d.png ' \
                  f'-c:v libx264 -pix_fmt yuv420p -crf 23 -t {bp_t_total} ' \
                  f'-r {param.frames_per_second} ' \
                  f'-vf pad=ceil(iw/2)*2:ceil(ih/2)*2 ' \
                  f'-y {os.path.join(video_dir, file_tag)}.mp4'.split()
        print_message('INFO', f'Creating the video: {command}', flush=True, log=log_file)
        subprocess.call(command)
        shutil.move(max_time_frame_file, f'{os.path.join(video_dir, file_tag)}.png')

        if param.qtfaststart is not None:
            command = f'{param.qtfaststart} {os.path.join(video_dir, file_tag)}.mp4 ' \
                      f'{os.path.join(video_dir, file_tag)}_q.mp4'.split()
            print_message('INFO', f'QTfaststart the video: {command}', flush=True, log=log_file)
            subprocess.call(command)
    except Exception as _er:
        print_message('ERR', f'Command {command} failed\n {_er}', flush=True, log=log_file)

    # Remove the pickle file.
    try:
        print_message('INFO', f'Removing the pickle file {pickle_file}', log=log_file)
        os.remove(pickle_file)
    except Exception as _er:
        print_message('ERR', f'Failed to remove\n {_er}', flush=True, log=log_file)

    # Remove image files if they exist from previous runs.
    try:
        files_to_remove = glob.glob(f'{os.path.join(scratch_dir, file_tag)}*.png')
        for this_file in files_to_remove:
            os.remove(this_file)
    except Exception as _er:
        print_message('ERR', f'Failed to remove\n {_er}', flush=True, log=log_file)

    return


"""Main code"""
"""
eq_date_time = '2020-07-22T06:12:44'
eq_datetime = UTCDateTime(eq_date_time)
eq_magnitude = 7.8
eq_lat = 55.030
eq_lon = -158.522
eq_depth = 10.0
vn_name = 'AU'
"""
event_id, vn_name, event_mag, event_date_time, event_lon, event_lat, event_depth = \
    (None, None, None, None, None, None, None)

# Get the input parameters.
try:
    options, remainder = getopt.getopt(sys.argv[1:], 'hlva:d:e:f:n:m:s:t:x:y:z:',
                                       ['anim', 'help', 'summary', 'verbose', 'eid=', 'factor=', 'vnet=',
                                        'emag=', 'decimate=', 'etime=', 'elon=', 'elat=', 'edepth=', 'logscreen='])
    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-a', '--anim'):
            if arg.strip().lower() == 'true':
                create_synthetics = True
                create_animation = True
            else:
                create_synthetics = False
                create_animation = False
        elif opt in ('-l', '--logscreen'):
            log_to_screen = True
        elif opt in ('-s', '--summary'):
            if arg.strip().lower() == 'true':
                create_summary = True
            else:
                create_summary = False
        elif opt in ('-v', '--verbose'):
            verbose = True
        elif opt in ('-e', '--eid'):
            event_id = arg.strip()
        elif opt in ('-f', '--factor'):
            try:
                grid_factor = float(arg.strip())
            except Exception as ex:
                usage()
                print_message('ERR', f'Invalid factor {arg.strip()}\n'
                              f'{ex}')
                sys.exit(2)
        elif opt in ('-n', '--vnet'):
            vn_name = arg.strip()
        elif opt in ('-m', '--emag'):
            try:
                event_mag = float(arg.strip())
            except Exception as ex:
                usage()
                print_message('ERR', f'Invalid magnitude {arg.strip()}\n'
                              f'{ex}')
                sys.exit(2)
        # The output sample rate will only exactly match the selected decimation rate if the ratio of original
        # to final rate is factorable by 2,3,5 and 7. Otherwise, the closest factorable rate will be chosen.
        elif opt in ('-d', '--decimate'):
            try:
                decimate = float(arg.strip())
            except Exception as ex:
                usage()
                print_message('ERR', f'Invalid sampling {arg.strip()}\n'
                              f'{ex}')
                sys.exit(2)
        elif opt in ('-t', '--etime'):
            event_date_time = arg.strip()
            event_date_time = event_date_time.replace(' ', 'T')
        elif opt in ('-x', '--elon'):
            try:
                event_lon = float(arg.strip())
            except Exception as ex:
                usage()
                print_message('ERR', f'Invalid longitude {arg.strip()}\n'
                              f'{ex}')
                sys.exit(2)
        elif opt in ('-y', '--elat'):
            try:
                event_lat = float(arg.strip())
            except Exception as ex:
                usage()
                print_message('ERR', f'Invalid latitude {arg.strip()}\n'
                              f'{ex}')
                sys.exit(2)
        elif opt in ('-z', '--edepth'):
            try:
                event_depth = arg.strip()
            except Exception as ex:
                usage()
                print_message('ERR', f'Invalid depth {arg.strip()}\n'
                              f'{ex}')
                sys.exit(2)
        else:
            print_message('WARN', f'option {opt} not recognized and will be ignored!')
except getopt.GetoptError as er:
    usage()
    print_message('ERR', f'\n\n{60 * "="}\n{er}\n{60 * "="}\n\n', flush=True)
    sys.exit(2)

# Validate parameters.
if event_id is not None and not any([event_mag, event_date_time, event_lon, event_lat, event_depth]):
    usage()
    print_message('ERR', f'Cannot provide both event ID and event parameter(s)')
    sys.exit(2)
elif event_id is None and not all([event_mag, event_date_time, event_lon, event_lat, event_depth]):
    usage()
    print_message('ERR', f'Must provide either an event ID or all event parameters')
    sys.exit(2)

if event_id is None:
    try:
        eq_magnitude = float(event_mag)
    except Exception as ex:
        usage()
        print_message('ERR', f'Invalid magnitude {event_mag}\n{ex}')
        sys.exit(2)

    try:
        eq_lon = float(event_lon)
        if eq_lon > 180:
            eq_lon -= 360
        if 180 < eq_lon or  eq_lon < -180:
            raise Exception('Longitude must be between -180/180')
    except Exception as ex:
        usage()
        print_message('ERR', f'Invalid longitude {event_lon}\n{ex}')
        sys.exit(2)

    try:
        eq_lat = float(event_lat)
        if 90 < eq_lat or eq_lat < -90:
            raise Exception('Latitude must be between -90/90')
    except Exception as ex:
        usage()
        print_message('ERR', f'Invalid latitude {event_lat}\n{ex}')
        sys.exit(2)

    try:
        eq_depth = float(event_depth)
    except Exception as ex:
        usage()
        print_message('ERR', f'Invalid depth {event_depth}\n{ex}')
        sys.exit(2)

    try:
        eq_date_time = event_date_time
        eq_datetime = UTCDateTime(eq_date_time)
    except Exception as ex:
        usage()
        print_message('ERR', f'Invalid date-time {eq_date_time}\n{ex}')
        sys.exit(2)

# Start logging.
if log_to_screen:
    log_file = sys.stdout
else:
    log_dir = lib.mkdir(param.log_dir)
    log_file_name = os.path.join(log_dir, script.replace('.py', ''))
    log_file_name = f"{log_file_name}_{vn_name}_{datetime.now().strftime('%Y-%m-%d')}"
    log_file_name = '.'.join([log_file_name, 'log'])
    log_file = open(log_file_name, 'a')
    error_file_name = '.'.join([log_file_name, 'err'])
    sys.stderr = open(error_file_name, 'a')

print_message('INFO', f'Event: M{eq_magnitude} at ({eq_lat}, {eq_lon}, {eq_depth}) {eq_date_time}UTC', log=log_file)
if vn_name not in vn.keys():
    usage()
    print_message('ERR', f'Invalid virtual network {vn_name}\n'
                  f'Must be one of {list(vn.keys())}')
    sys.exit(2)

bp_filter = param.bp_filter[vn_name]
max_sta_count = param.virtual_networks[vn_name]['max_sta_count']
sta_too_close_deg = param.virtual_networks[vn_name]['sta_too_close_deg']

if not any([create_animation, create_summary, create_synthetics]):
    usage()
    print_message('ERR', f'All outputs are off, please turn at least one on '
                  f'(create_animation, create_summary, create_synthetics)')
    sys.exit(2)

if timing:
    t0 = time()
    t1 = lib.time_it(t0, 'START', log=log_file)

# Get the waveform request interval.
request_start_date_time, request_start_datetime, request_end_date_time, request_end_datetime = \
    lib.get_bp_time_window(eq_date_time, eq_magnitude)

# Get event's list of bulk requests using fedcatalog.
eq_content = lib.get_fedcatalog_stations(request_start_date_time, request_end_date_time,
                                         eq_lat, eq_lon, eq_min_radius, eq_max_radius,
                                         net=vn[vn_name]['network'], req='text', service='station', log=log_file)

if eq_content is None:
    print_message('ERR', f'get_fedcatalog_stations did not return any station list. Returned {eq_content} ')
    sys.exit(1)

eq_stations = lib.split_fedcatalog_stations(eq_content, log=log_file)

if eq_content is None:
    print_message('ERR', f'Station request for event failed, no stations to request!', log=log_file)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None, anim_tag='syn')
    sys.exit(4)

print_message('INFO', f'Requesting data from {len(eq_stations)} stations.', log=log_file)

# Find the azimuth and distance from the earthquake to the center of the continent the BP is made for.
center_dist, center_azim, center_back_azim = gps2dist_azimuth(eq_lat, eq_lon, vn[vn_name]['lat'], vn[vn_name]['lon'])
print_message('INFO', f'Center of the {vn_name} network ({vn[vn_name]["lat"]}, {vn[vn_name]["lon"]}) azimuth is '
                      f'{center_azim} degrees.', log=log_file)

# Distance is the great circle distance in m, convert to km.
center_dist /= 1000.0

# Before we start, check to see if this network is at a reasonable distance from the event.
center_dist_degrees = kilometer2degrees(center_dist)
if param.vn_check_dist and center_dist_degrees > param.qc_vn_max_distance and \
        vn_name not in param.vn_azimuth_exception_list:
    print_message('ERR', f'{vn_name} virtual network is {center_dist_degrees:0.2f} '
                         f'degrees from the earthquake location.'
                         f' Too far (> {param.qc_vn_max_distance}) for back Projection!', log=log_file)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None, anim_tag='syn')
    sys.exit(4)

bulk_request_info = dict()
station_coordinates = dict()
chunk = 0

# Go through each contributing FDSN data center (DC) and clean up the station list and create  actual data requests.
# eq_stations is a DC-based dictionary.
request_list = dict()

for dc in eq_stations:
    _service = None
    _list = list()

    for _line in eq_stations[dc].bulk:
        items = _line.split('|')
        _net, _sta, _loc, _chan, _lat, _lon = items[0:6]
        net_sta_key = f'{_net}.{_sta}'

        # Do not include, if already in the list:
        if net_sta_key in request_list:
            continue

        if _loc.strip() == '':
            _loc = '--'
        _chan = items[3]
        _lat = float(_lat)
        _lon = float(_lon)
        _dist, _azim, _back_azim = gps2dist_azimuth(eq_lat, eq_lon, _lat, _lon)
        _dist /= 1000.0
        _dist = kilometer2degrees(_dist)

        # Only accept stations that are within vn_azimuth degrees from the center line, except for selected networks.
        if vn_name not in param.vn_azimuth_exception_list:
            _angle = abs(center_azim - _azim)

            # See if the station azimuth is within the range.
            if _angle <= vn_azimuth or _angle >= 360.0 - vn_azimuth:
                if verbose:
                    print_message('INFO', f'Station {_net}.{_sta} azimuth {_azim:0.2f} '
                                          f'from center {center_azim:0.2f} is '
                                          f'{_angle:0.2f} and is within the '
                                          f'{vn_azimuth} bounds.', log=log_file)
            else:
                if verbose:
                    print_message('WARN', f'Station {_net}.{_sta} azimuth {_azim:0.2f} '
                                          f'from center {center_azim:0.2f} is '
                                          f'{_angle:0.2f} and is outside the '
                                          f'{vn_azimuth} bounds.', log=log_file)
                continue

        request_list[net_sta_key] = dict()
        request_list[net_sta_key]['lat'] = _lat
        request_list[net_sta_key]['lon'] = _lon

dense_patch_list = ''
sparsify = 0
sta_too_close_deg_inc = param.sta_too_close_deg_inc
search_radius_edge = lib.case(len(request_list), param.sta_too_close_deg_init)
print_message('INFO', f'There are {len(request_list)} stations in the virtual network, setting search radius to '
                      f'{search_radius_edge} degrees, increasing by {sta_too_close_deg_inc} degrees', log=log_file)

awhile = ' (this may take a while)'

# Do not set the initial search radius beyond the allowed sta_too_close_deg value for the virtual network.
if search_radius_edge > sta_too_close_deg:
    search_radius_edge = sta_too_close_deg

if len(request_list) <= max_sta_count:
    print_message('INFO', f'{len(request_list)} stations in the virtual network <= {max_sta_count}', log=log_file)
else:
    while search_radius_edge <= sta_too_close_deg:
        # Make the network sparse for this search radius.
        while dense_patch_list is not None:
            sparsify += 1
            print_message('INFO', f'Pre-request sparse iteration #{sparsify}, looking for dense patches '
                                  f'between {len(request_list)} stations within {search_radius_edge:0.2} degrees'
                                  f'{awhile}',
                          log=log_file)
            awhile = ''
            intra_station_dist, dense_patch_list = lib.find_dense_sta_patch(intra_station_dist, vn_name, request_list,
                                                                            search_radius_edge, log=log_file)
            if dense_patch_list is not None:
                for key in dense_patch_list:
                    if dense_patch_list[key]:
                        print_message('WARN', f'{key} removed to make network sparse.', log=log_file)
                        request_list.pop(key)
                if len(request_list) <= max_sta_count:
                    print_message('INFO', f'Pre-request sparse done!', log=log_file)
                    break
            else:
                print_message('INFO', f'Network {vn_name} with {len(request_list)} stations is sparse within '
                                      f'{search_radius_edge:0.2} degrees', log=log_file)

        # Is the number of remaining stations still too high?
        if len(request_list) > max_sta_count and (search_radius_edge < sta_too_close_deg or sta_too_close_deg <= 0):
            search_radius_edge += sta_too_close_deg_inc
            sparsify = 0

            # Do not set the  search radius beyond the allowed sta_too_close_deg value for the virtual network.
            if search_radius_edge > sta_too_close_deg > 0:
                if search_radius_edge > sta_too_close_deg:
                    search_radius_edge = sta_too_close_deg
            print_message('INFO', f'We still have {len(request_list)} > {max_sta_count} '
                                  f'stations in the virtual network '
                                  f'{vn_name}, increasing the search radius to {search_radius_edge:0.2} degrees.',
                          log=log_file)
            dense_patch_list = ''
        else:
            print_message('INFO', f'Pre-request sparse done!', log=log_file)
            break

print_message('INFO', f'The search area is st at {search_radius_edge} degrees!', log=log_file)
if timing:
    t1 = lib.time_it(t1, f'Pre-request sparse', log=log_file)

dc_service_url = dict()
for dc in eq_stations:
    _service = None
    _list = list()

    for _line in eq_stations[dc].bulk:
        items = _line.split('|')
        _net, _sta, _loc, _chan, _lat, _lon = items[0:6]
        if _loc.strip() == '':
            _loc = '--'
        _chan = items[3]
        _lat = float(_lat)
        _lon = float(_lon)
        net_sta_key = f'{_net}.{_sta}'
        if net_sta_key not in request_list:
            continue

        # Add the station to the request list.
        _list.append(' '.join([_net, _sta, _loc, _chan, request_start_date_time, request_end_date_time]))

        # Do we know the dataselect service address for this DC? If not, get it using this station.
        if _service is None:
            _service = lib.get_dc_dataselect_url(request_start_date_time, request_end_date_time,
                                                 _net, _sta, _loc, _chan, dc, log=log_file)
            dc_service_url[dc] = _service

    # We got the list for this DC. Clean it up and break it into chunks for easier request.
    _req = lib.get_request_items(_list)
    print_message('INFO', f'DC {dc}_{chunk}, {len(_req)} channels  out of possible {len(_list)}', log=log_file)
    if timing:
        t1 = lib.time_it(t1, f'Station list for {dc}', log=log_file)

    chunk_list = list()
    for _item in _req:
        chunk_list.append(_item)
        if len(chunk_list) >= chunk_count:

            bulk_request_info[f'{dc}_{chunk}'] = ObjDict({'url': eq_stations[dc].url,
                                                          'dataselect': _service,
                                                          'bulk': chunk_list.copy()})
            chunk += 1
            chunk_list = list()

    # Reached the end for this DC. Create request for the remaining, if any.
    if chunk_list:
        bulk_request_info[f'{dc}_{chunk}'] = ObjDict({'url': eq_stations[dc].url,
                                                      'dataselect': _service,
                                                      'bulk': chunk_list.copy()})
        chunk += 1
        chunk_list = list()

# Make data request to each DC.
trace_list = dict()
for dc in bulk_request_info:
    data_center = dc.split('_')[0]
    st = None
    if data_center in dc_to_exclude:
        print_message('WARN', f'skipped data from {data_center} because it is in the exclude list', log=log_file)
        continue
    else:
        if verbose:
            print_message('INFO', 'Sending requests for:', log=log_file)
            for line in bulk_request_info[dc].bulk:
                print(f'\t{line}\n', file=log_file)
            print('\n', file=log_file, flush=True)

        print('\n', file=log_file, flush=True)
        if len(bulk_request_info[dc].bulk) <= 0:
            print_message('WARN', f'Skipping DC {dc}, no stations to request!\n', log=log_file)
            continue

        print_message('INFO', f'Requesting data from {data_center} via '
                              f'{bulk_request_info[dc].dataselect}\n', log=log_file)

    # Set the client up for this data center.
    try:
        print_message('INFO', f'Sending {len(bulk_request_info[dc].bulk)} requests', log=log_file)
        client = Client(lib.get_service_url(bulk_request_info, dc))
        st = client.get_waveforms_bulk(bulk_request_info[dc].bulk, attach_response=True)
        print_message('INFO', f'Received {len(st)} traces', log=log_file)
    except Exception as er:
        if verbose:
            print_message('ERR', f'Request:\n\n{bulk_request_info[dc].bulk}\n\n from {data_center} '
                                 f'({bulk_request_info[dc].dataselect}) failed: {er}', log=log_file)
        else:
            print_message('ERR', f'({bulk_request_info[dc].dataselect}) failed: {er}', log=log_file)
        continue

    # Fill the missing values with zeros.
    if param.merge_gaps:
        st.merge(fill_value=0.0)

    # Work on individual traces in the stream.
    for index, tr in enumerate(st):

        # Reject traces with NaN and inf values.
        if True in np.isnan(tr.data) or True in np.isinf(tr.data):
            if verbose:
                print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta '
                                      f'{net_sta_key} from {data_center} '
                                      f'because of bad values', log=log_file)
            continue

        # Get the individual trace times.
        trace_times = list(tr.times())

        # Calculate trace length based on the time difference between the last and the first samples.
        trace_length = trace_times[-1] - trace_times[0]
        net_sta_key = f'{tr.stats.network}.{tr.stats.station}'

        # Ignore the short traces.
        min_trace_length = 0.9 * (request_end_datetime - request_start_datetime)
        if trace_length < min_trace_length:
            if verbose:
                print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta_key} '
                                      f'from {data_center} '
                                      f'because it is shorter than {min_trace_length} s', log=log_file)
            continue

        # Ignore stations with too high of STD before event time.
        if param.std_check:
            std_start = (param.request_time_before - param.std_offset - param.std_window) * tr.stats.sampling_rate
            std_end = (std_start + param.std_window) * tr.stats.sampling_rate
            _tr = tr.copy()
            _tr.detrend()
            _data = _tr.data
            trace_max = max(abs(_data))
            std_data = _data[int(std_start):int(std_end)] / trace_max
            std = np.std(np.array(std_data))
            if std > vn[vn_name]['std_max']:
                print_message(f'WARN',
                              f'Rejected {data_center}.{net_sta_key} due to high STD of '
                              f'{std:0.4f} > {vn[vn_name]["std_max"]}', log=log_file)
                continue

        # Get the station coordinates. Here we assume there is a possibility that station contains gaps, so
        # we may get more than one trace. We will get the station information from the first segment.
        # NOTE: We ignore the second trace
        if net_sta_key in station_coordinates:
            if verbose:
                print_message('WARN', f'Multiple trace, already have data from {net_sta_key} '
                                      f'for channel {station_coordinates[net_sta_key][2]}.', log=log_file)
            continue

        else:
            if verbose:
                print_message('INFO', f'Getting information for station {tr.stats.network}.{tr.stats.station}.'
                                      f'{tr.stats.location}.{tr.stats.channel}.', log=log_file)
        try:

            inventory = client.get_stations(network=tr.stats.network, station=tr.stats.station,
                                            location=tr.stats.location, channel=tr.stats.channel, level="station")
            _lon = inventory.networks[0].stations[0].longitude
            _lat = inventory.networks[0].stations[0].latitude
            station_coordinates[net_sta_key] = (_lon, _lat, tr.stats.channel)
        except Exception as er:
            print_message('ERR', f'Request error {data_center}.{net_sta_key} failed: {er}', log=log_file)
            continue

        # Make sure seismogram has the P-waves in it.
        _dist, _azim, _back_azim = gps2dist_azimuth(eq_lat, eq_lon, _lat, _lon)

        # Delay from predicted travel time.
        tt_cache, phase_delay = lib.get_phase_delay(tt_cache, _dist, eq_depth)

        if phase_delay is None:
            if verbose:
                print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta_key} '
                                      f'from {data_center} '
                                      f'because no P-wave arrival', log=log_file)
            continue

        if not lib.has_phase(tr, eq_datetime, phase_delay):
            if verbose:
                print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta_key} '
                                      f'from {data_center} '
                                      f'because it does not cover the {param.phase} arrival time', log=log_file)
            continue

        # Remove the response.
        try:
            if verbose:
                print_message('INFO', f'Removing response from {data_center}.{net_sta_key}', log=log_file)
            tr.remove_response(output=param.output_type, zero_mean=True, taper=True,
                               taper_fraction=param.taper_fraction, pre_filt=param.pre_filter)
        except ValueError as er:
            print_message('ERR', f'Removing response from {data_center}.{net_sta_key} failed: {er}', log=log_file)
            continue
        except Exception as er:
            print_message('ERR', f'Removing response from {data_center}.{net_sta_key} failed: {er}', log=log_file)
            continue

        # To make sure all is OK after response correction, reject traces with NaN and inf values.
        if True in np.isnan(tr.data) or True in np.isinf(tr.data):
            if verbose:
                print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta_key} '
                                      f'from {data_center} '
                                      f'because of bad values', log=log_file)
            continue

        # Demean the trace (#587).
        tr.detrend(type='demean')

        # Taper the trace.
        tr.taper(param.taper_fraction, type='hann')

        # Make sure pre-Phase to phase SNR is at least 2 on the unfiltered trace. (#525)
        snr = lib.p_wave_snr(tr, eq_datetime, phase_delay)
        if snr < param.min_snr:
            print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta_key} from {data_center} '
                                  f'because its P-wave SNR is {snr:0.2f} < {param.min_snr} ', log=log_file)
            continue
        else:
            print_message('INFO', f'Channel {tr.stats.channel} of Net.Sta {net_sta_key} from {data_center} '
                                  f'has acceptable unfiltered P-wave SNR of {snr:0.2f} > {param.min_snr} ',
                          log=log_file)

        # Make sure all traces have the same sampling (trace_sampling_frequency).
        try:
            if not math.isclose(tr.stats.delta,  trace_sampling_frequency):
                tr.resample(trace_sampling_frequency)
        except Exception as ex:
            print_message('ERR', f'Failed to re-sample {data_center}.{net_sta_key} to '
                                 f'{param.trace_sampling_frequency}: '
                                 f'{ex}', log=log_file)
            continue

        # Filter/normalize the trace (#715).
        tr_filter = lib.preprocess_trace(tr, bp_filter, eq_datetime, phase_delay)

        # No filter, just normalize and demean.
        tr = lib.preprocess_trace(tr, None, eq_datetime, phase_delay)

        # One last time, make sure pre-Phase to phase SNR is at least 1.5 on the filtered trace.(#added)
        if param.min_filtered_snr is not None:
            tr_ = tr_filter

            snr = lib.p_wave_snr(tr_, eq_datetime, phase_delay)
            if snr < param.min_filtered_snr:
                print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta_key} '
                                      f'from {data_center} '
                                      f'because its {bp_filter} filtered P-wave SNR is '
                                      f'{snr:0.2f} < {param.min_filtered_snr} ', log=log_file)
                continue
            else:
                print_message('INFO', f'Channel {tr.stats.channel} of Net.Sta {net_sta_key} from {data_center} '
                                      f'has {bp_filter} filtered P-wave SNR of '
                                      f'{snr:0.2f} > {param.min_filtered_snr} ', log=log_file)

            # Phase align traces (af2 is the tapered / filtered and normalized trace).
            trace_list[net_sta_key] = {'tr_final': tr.copy(), 'tr_filter': tr_filter.copy(),
                                       'weight': 0.0, 'dc': data_center,
                                       'lat': inventory.networks[0].stations[0].latitude,
                                       'lon': inventory.networks[0].stations[0].longitude, 'azim': _azim,
                                       'phase_delay': phase_delay}
    if timing:
        t1 = lib.time_it(t1, f'Pre-processed station data for {data_center}', log=log_file)

print_message('INFO', f'Pre-processing done', log=log_file)

# BP parameters.
stack_start, stack_end, bp_t_offset, bp_t_increment, bp_t_total, t_avg = lib.set_time_parameters(eq_magnitude)

# See if we have enough stations to continue.
if len(trace_list) < param.min_num_sta:
    if verbose:
        print_message('ERR', f'Only {len(trace_list)} stations remain, must be >= {param.min_num_sta},  '
                             f'will not continue', log=log_file)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, search_radius_edge)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, search_radius_edge, anim_tag='syn')
    sys.exit(4)

# remove dense patch of stations (#764).
sparse_title = ''
sparse_title = ' (sparse)'

dense_patch_list = ''
sparsify = 0

while dense_patch_list is not None:
    sparsify += 1
    print_message('INFO', f'Sparse iteration #{sparsify} looking for dense patches '
                          f'between {len(trace_list)} stations within {search_radius_edge} degrees', log=log_file)
    intra_station_dist, dense_patch_list = lib.find_dense_sta_patch(intra_station_dist, vn_name, trace_list,
                                                                    search_radius_edge, log=log_file)
    if dense_patch_list is not None:
        for key in dense_patch_list:
            if dense_patch_list[key]:
                print_message('WARN', f'{key} removed to make network sparse.', log=log_file)
                trace_list.pop(key)
    else:
        print_message('INFO', f'Network is sparse', log=log_file)

if timing:
    t1 = lib.time_it(t1, f'Sparse', log=log_file)

# Compute the cross-correlation window length (#734).
cc_window_length = lib.xcorr_window(trace_list, eq_magnitude, eq_datetime)
if timing:
    t1 = lib.time_it(t1, f'X-corr window set to {cc_window_length:0.2f}s', log=log_file)

# Optimal cross-correlation.
print_message('INFO', f'Computing optimal MCCC window for {len(trace_list)} stations', log=log_file)

trace_list, nsta, optimal_mccc = lib.find_optimal_mccc_window(eq_datetime, trace_list, vn[vn_name], log=log_file)

if optimal_mccc is None:
    print_message('ERR', f' Insufficient number of traces ({nsta} < {param.min_num_sta}) '
                         f'to set the optimal MCCC window!', log=log_file)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, search_radius_edge)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, search_radius_edge, anim_tag='syn')
    sys.exit(4)

if timing:
    t1 = lib.time_it(t1, f'Optimal MCC window set ', log=log_file)

# Double cross-correlation.
#trace_list, double_mccc = lib.mccc_double(eq_datetime, trace_list, cc_window_length)
#if timing:
#    t1 = lib.time_it(t1, f'Double MCCC done', log=log_file)

# Assign station weights #816.
print_message('INFO', f'Assign station weights', log=log_file)

# Compute a weight to each active station based on the distance and azimuth of the neighboring stations.
intra_station_dist, weight = lib.station_weight(trace_list, vn_name, intra_station_dist)
weight_count = 0
for net_sta_key in weight.keys():
    if weight[net_sta_key] > 0.0:
        weight_count += 1

# See if we have enough stations to continue.
if weight_count < param.min_num_sta:
    if verbose:
        print_message('WARN', f'Only {weight_count} stations remain, must be >= {param.min_num_sta},  will not continue'
                      , log=log_file)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, search_radius_edge)
    make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, search_radius_edge, anim_tag='syn')
    sys.exit(4)

for net_sta_key in weight.keys():
    if weight[net_sta_key] > 0.0:
        trace_list[net_sta_key]['weight'] = weight[net_sta_key]

        # Apply the station weight to its traces.
        trace_list[net_sta_key]['tr_final'].data = trace_list[net_sta_key]['tr_final'].data * weight[net_sta_key]
        trace_list[net_sta_key]['tr_filter'].data = trace_list[net_sta_key]['tr_filter'].data * weight[net_sta_key]
    else:
        trace_list.pop(net_sta_key)
if timing:
    t1 = lib.time_it(t1, f'Station weights assigned', log=log_file)

# Cut or zero pad the seismograms so that all of them start at origin time and have the MCCC delay added in.
# Now, any delay due to change in the source location can be applied to these traces.
print_message('INFO', f'Trim traces to start {bp_t_offset}s before P and continue for {bp_t_total} seconds.',
              log=log_file)

# Set the trace start and length.
_t0 = max(bp_t_offset, - param.min_summary_trace_time)
_t = max(bp_t_total, _t0 + param.max_summary_trace_time)
for tr_key in trace_list.keys():

    _tr = trace_list[tr_key]['tr_final'].copy()
    _tr_filter = trace_list[tr_key]['tr_filter'].copy()

    trace_list[tr_key]['tr_final'] = lib.trace_trimmer(_tr, eq_datetime, _t0, _t,
                                                       trace_list[tr_key]['phase_delay'],
                                                       trace_list[tr_key]['mccc_delay']
                                                       )
    trace_list[tr_key]['tr_filter'] = lib.trace_trimmer(_tr_filter, eq_datetime, _t0, _t,
                                                        trace_list[tr_key]['phase_delay'],
                                                        trace_list[tr_key]['mccc_delay']
                                                        )

if create_synthetics:
    synthetics = lib.gen_synthetics(trace_list, eq_datetime, vn_name, create_synthetics=create_synthetics)
    # This synthetics should have the same time history as the
    # regular trace, so should follow the same processing steps.
    for tr_key in trace_list.keys():
        trace_list[tr_key]['tr_syn'] = synthetics[tr_key].copy()

debug = False
if debug:
    for tr_key in trace_list.keys():
        trace_1 = trace_list[tr_key]['tr_filter']
        trace_2 = trace_list[tr_key]['tr_syn']
        stream = Stream(traces=[trace_1, trace_2])
        stream.plot()

# Stack traces.
if timing:
    t1 = lib.time_it(t1, f'Stack Start', log=log_file)
if True:
    stack, syn_stack = lib.stacker(trace_list, vn_name, bp_t_offset, t_avg, eq_datetime, eq_lat, eq_lon, eq_magnitude,
                                   eq_depth, bp_t_increment, bp_t_total, tt_cache, create_synthetics=create_synthetics,
                                   grid_factor=grid_factor, log=log_file, verbose=verbose)
    if timing:
        t1 = lib.time_it(t1, f'Stack END', log=log_file)

    # 1760
    # Raise stacks to the Nth power to do Nth-root stacking
    debug = False
    if debug:
        print_message('DEBUG', f'Stack before Nth power', log=log_file)
        plot_stack(stack)
        print_message('DEBUG', f'Syn stack before Nth power', log=log_file)
        plot_stack(syn_stack)
    stack, syn_stack = lib.stack_root(stack, syn_stack, param.stacking_root, create_synthetics=create_synthetics)

    if debug:
        print_message('DEBUG', f'Stack after Nth power', log=log_file)
        plot_stack(stack)
        print_message('DEBUG', f'Syn stack after Nth power', log=log_file)
        plot_stack(syn_stack)

    if timing:
        t1 = lib.time_it(t1, f'Stack_root done', log=log_file)

    # 1771
    # Average the beams by averaging the beams in a running beam_average_seconds window
    print_message('INFO', f'Averaging the beams', log=log_file)
    stack, syn_stack = lib.hann(stack, syn_stack, stack_start, stack_end, bp_t_increment,
                                param.beam_average_seconds,
                                create_synthetics=create_synthetics, log=log_file)

    if debug:
        print_message('DEBUG', f'Running average', log=log_file)
        plot_stack(stack)
        print_message('DEBUG', f'Syn Running average', log=log_file)
        plot_stack(syn_stack)

    if timing:
        t1 = lib.time_it(t1, f'Beam averaging done', log=log_file)

    # 1815
    # Time average the stacks.
    print_message('INFO', f'Time averaging the stacked traces and final resampling', log=log_file)

    # The output sample rate will only exactly match the selected decimation rate if the ratio of original to final
    # rate is factorable by 2,3,5 and 7. Otherwise, the closest factorable rate will be chosen.
    step = 1
    if decimate is not None:
        step = int(round(decimate / bp_t_increment))

    print_message('INFO', f'Stack decimation every {step} samples', log=log_file)
    # Average the beams and do final averaging using a running t_avg window

    stack_final, syn_stack_final = lib.hann(stack, syn_stack, stack_start, stack_end, bp_t_increment, t_avg,
                                            resamp=step, create_synthetics=create_synthetics, log=log_file)

    if debug:
        print_message('DEBUG', f'Time average', log=log_file)
        plot_stack(stack)
        print_message('DEBUG', f'Syn Time average', log=log_file)
        plot_stack(syn_stack)

    t1 = lib.time_it(t1, f'Time averaging done', log=log_file)

    # 1879
    # Compute power, square the stacks.
    stack_final, syn_stack_final = lib.stack_root(stack_final, syn_stack_final, 2,
                                                  create_synthetics=create_synthetics)

    if timing:
        t1 = lib.time_it(t1, f'Stack root done', log=log_file)

    stack_amp = dict()
    stack_amp_loc = dict()
    stack_max = dict()
    stack_median = dict()
    stack_max_loc = dict()
    global_max = None

    syn_stack_amp = dict()
    syn_stack_amp_loc = dict()
    syn_stack_max = dict()
    syn_stack_max_loc = dict()
    syn_global_max = None

    # Record amplitude and location of the stack at each grid point for the given time key.
    amp_tag = f'BP_PeakAmps_{param.trace_filter[param.bp_filter[vn_name]]["label"].replace(" ", "_")}'
    amp_file_tag = '_'.join([amp_tag, vn_name, lib.file_name_tag(eq_datetime)])
    amp_file_name = f'{amp_file_tag}.txt'

    sta_tag = f'BP'
    sta_file_tag = '_'.join([sta_tag, vn_name, lib.file_name_tag(eq_datetime)])
    sta_list_file_name = f'{sta_file_tag}_stations.txt'
    with open(os.path.join(data_dir, sta_list_file_name), 'w') as fp:
        for _key in trace_list:
            _dc = trace_list[_key]['dc']
            if _dc not in inventory_list:
                inventory_list[_dc] = list()
            _tr = trace_list[_key]['tr_final']
            _loc = _tr.stats.location.strip()
            if not _loc:
                _loc = '--'
            inventory_list[_dc].append(f'{_tr.stats.network} {_tr.stats.station} {_loc} '
                                       f'{_tr.stats.channel} {request_start_date_time} '
                                       f'{request_end_date_time}')
        for _dc in inventory_list:
            fp.write(f'DATACENTER={_dc}\nDATASELECTSERVICE={dc_service_url[_dc]}\n')
            for _sta in inventory_list[_dc]:
                fp.write(f'{_sta}\n')

    amp_file = os.path.join(data_dir, amp_file_name)
    fp = open(amp_file, 'w')
    if create_synthetics:
        fp.write(f'{"Time":>10}{"Lat":>10}{"Lon":>10}{"DataPeak":>15}{"DataMedian":>15}{"ARFPeak":>15}\n')
    else:
        fp.write(f'{"Time":>10}{"Lat":>10}{"Lon":>10}{"DataPeak":>15}{"DataMedian":>15}\n')

    # Create time slices.
    for grid_key in stack_final:
        for t_key in stack_final[grid_key]:
            if float(t_key) < stack_start or float(t_key) > stack_end:
                continue
            if t_key not in stack_amp:
                stack_amp[t_key] = list()
                stack_amp_loc[t_key] = list()

                # Synthetic.
                if create_synthetics:
                    if t_key not in syn_stack_amp:
                        syn_stack_amp[t_key] = list()
                        syn_stack_amp_loc[t_key] = list()

            _val = stack_final[grid_key][t_key]
            stack_amp[t_key].append(_val)
            _lat, _lon = grid_key.split('_')
            stack_amp_loc[t_key].append((float(_lat), float(_lon)))

            # Synthetic.
            if create_synthetics:
                _val = syn_stack_final[grid_key][t_key]
                syn_stack_amp[t_key].append(_val)
                _lat, _lon = grid_key.split('_')
                syn_stack_amp_loc[t_key].append((float(_lat), float(_lon)))

    # Here we scan all the values for a given time key and find the maximum save its value and location.
    for t_key in stack_amp:
        stack_max[t_key] = max(stack_amp[t_key])
        _index = stack_amp[t_key].index(stack_max[t_key])
        stack_max_loc[t_key] = stack_amp_loc[t_key][_index]
        stack_median[t_key] = np.median(stack_amp[t_key])

        # Synthetic.
        if create_synthetics:
            syn_stack_max[t_key] = max(syn_stack_amp[t_key])
            _index = syn_stack_amp[t_key].index(syn_stack_max[t_key])
            syn_stack_max_loc[t_key] = syn_stack_amp_loc[t_key][_index]
    # Find the global maximum.
    global_max = max(stack_max.values())

    global_max_location = None
    for _key, _value in stack_max.items():
        if _value == global_max:
            global_max_location = stack_max_loc[_key]
            break

    # The final QC to make sure BP for this VN is reasonable. The criteria is just a simple
    # check to see if the peak maxima is too far from the event location.
    if global_max_location is not None:
        print_message('INFO', f'Global max location {global_max_location}', log=log_file)
        _dist, _azim, _back_azim = gps2dist_azimuth(eq_lat, eq_lon,
                                                    global_max_location[0],
                                                    global_max_location[1])
        if _dist / 1000.0 > param.peak_offset_max_km:
            print_message('ERR', f'The global max is {_dist / 1000.0:0.2f} km from the event > '
                                 f'{param.peak_offset_max_km}', log=log_file)
            print_message('ERR', f'Station request for event failed, no stations to request!', log=log_file)
            make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None)
            make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None, anim_tag='syn')
            sys.exit(4)
        else:
            print_message('INFO', f'The global max is {_dist / 1000.0:0.2f} km from the event < '
                                  f'{param.peak_offset_max_km}', log=log_file)

    else:
        print_message('WARN', f'Could not locate the  global max  {global_max:0.2f} km', log=log_file)

    # Synthetic.
    if create_synthetics:
        syn_global_max = max(syn_stack_max.values())

    for t_key in stack_max:
        if create_synthetics:
            fp.write(f'{t_key:>10}{stack_max_loc[t_key][0]:>10.3f}{stack_max_loc[t_key][1]:>10.3f}'
                     f'{stack_max[t_key]:>15.5e}{stack_median[t_key]:>15.5e}{syn_stack_max[t_key]:>15.5e}\n')
        else:
            fp.write(f'{t_key:>10}{stack_max_loc[t_key][0]:>10.3f}{stack_max_loc[t_key][1]:>10.3f}'
                     f'{stack_max[t_key]:>15.5e}{stack_median[t_key]:>15.5e}\n')
    fp.close()

    if timing:
        t1 = lib.time_it(t1, f'Maximum amplitude found', log=log_file)

    # Compute a cumulative stack.
    cumulative_stack = dict()
    cumulative_stack_max = None
    for grid_key in stack_final:
        if grid_key not in cumulative_stack:
            cumulative_stack[grid_key] = 0.0
            for time_key in stack_final[grid_key]:
                cumulative_stack[grid_key] = np.add(cumulative_stack[grid_key], stack_final[grid_key][time_key])
            if cumulative_stack_max is None:
                cumulative_stack_max = abs(cumulative_stack[grid_key])
            else:
                cumulative_stack_max = max(cumulative_stack_max, abs(cumulative_stack[grid_key]))
    if timing:
        t1 = lib.time_it(t1, f'Cumulative stack computed', log=log_file)

# Create animation  for the virtual network, if requested.
if create_animation:
    print_message('INFO', f'Creating the animation', log=log_file)
    make_animation(stack_start, stack_end, stack_amp, stack_amp_loc, stack_max, stack_max_loc, global_max,
                   search_radius_edge, grid_factor=grid_factor)
    if create_synthetics:
        make_animation(stack_start, stack_end, syn_stack_amp, syn_stack_amp_loc, syn_stack_max, syn_stack_max_loc,
                       syn_global_max, search_radius_edge, grid_factor=grid_factor, anim_tag='syn')

# Create summary Plot for the virtual network.
# Initialise the summary plot.
if create_summary:
    # What type of media is generated?
    media = 'image'

    fig = plt.figure(figsize=param.figure_size, facecolor='white')
    fig.subplots_adjust(top=0.8)

    # Figure layout.
    subplot_columns = 2
    subplot_tall_rows = 2
    subplot_short_rows = 1
    tall_to_short_height = 3

    # 1. Top Left Subplot, stack sum.
    print_message('INFO', f'Stack sum', log=log_file)
    ax7 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                           (0, 0), rowspan=tall_to_short_height, colspan=1)

    # Extract grid and values.
    lat = list()
    lon = list()
    val = list()

    for grid_key in cumulative_stack:
            _lat, _lon = grid_key.split('_')
            lat.append(float(_lat))
            lon.append(float(_lon))
            value = cumulative_stack[grid_key] / cumulative_stack_max
            val.append(value)

    # Must be np arrays for grid
    lat_list = np.array(lat)
    lon_list = np.array(lon)
    value_list = np.array(val)

    # Find the min and max of coordinates.
    lon_min = lon_list.min()
    lat_min = lat_list.min()
    lon_max = lon_list.max()
    lat_max = lat_list.max()

    latitude, longitude = lib.set_grid(eq_lat, eq_lon, eq_magnitude, grid_factor=grid_factor)

    # Now let's grid the data. Find the number of grid points in each direction.
    lon_num = pcolormesh_grid_factor * int(((lon_max - lon_min) / longitude['inc']) + 1)
    lat_num = pcolormesh_grid_factor * int(((lat_max - lat_min) / latitude['inc']) + 1)

    width = lat_max - lat_min
    width = degrees2kilometers(width) * 1000.0
    lon_0 = eq_lon
    lat_0 = eq_lat
    bm = Basemap(width=width, height=width, projection=basemap_projection,
                 lat_0=lat_0, lon_0=lon_0, resolution=param.basemap_resolution)

    coast_alpha = 1
    if param.fill_continents:
        coast_alpha = param.fill_continents_alpha
        bm.fillcontinents(color=param.fill_continents_color, alpha=coast_alpha)

    # Network's name on the upper left
    plt.text(0.1, 0.95, vn_name, fontsize=font_size[media]['network'], horizontalalignment='center',
             verticalalignment='center', transform=ax7.transAxes, color='red',
             path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    scalebar_deg = kilometer2degrees(param.scalebar_km)
    _lat0 = lat_min + (lat_max - lat_min) / 4.0
    _lon0 = lon_max - (lon_max - lon_min) / 10.0
    _lat1, _lon1 = lib.get_location(_lat0, _lon0, 0, scalebar_deg)
    _x, _y = bm((_lon0, _lon0), (_lat0, _lat1))
    print_message('INFO', f'Map scale between: ({_lat0:0.2f}, {_lon0:0.2f}) and ({_lat1:0.2f}, {_lon1:0.2f}) / '
                          f'({_x[0]:0.2f}, {_y[0]:0.2f}) and ({_x[1]:0.2f}, {_y[1]:0.2f})', log=log_file)
    # Use the same _X to ensure scale is vertical.
    bm.plot([_x[0], _x[0]], _y, color='black', linewidth=2)
    plt.text(_x[0], _y[1], f'  {param.scalebar_km} km', fontsize=font_size[media]['legend'],
             horizontalalignment='center', rotation=90,
             verticalalignment='bottom',
             path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    # Create a uniform mesh for contouring. First transfer lon, lat to map units (x, y).
    x_old, y_old = bm(lon_list, lat_list)
    x_new = np.linspace(min(x_old), max(x_old), lon_num)
    y_new = np.linspace(min(y_old), max(y_old), lat_num)

    # Basic mesh in map's x, y.
    grid_x, grid_y = np.meshgrid(x_new, y_new)

    try:
        # Interpolate at the points in lon_new, lat_new.
        # Method : {'linear', 'nearest', 'cubic'}, optional.
        grid_v = griddata((x_old, y_old), value_list, (grid_x, grid_y), method='cubic')
    except Exception as _er:
        print_message('ERR', f'Griding failed: {_er}', log=log_file)
        sys.exit(2)

    # Make the custom color map.
    bp_cmap = lib.make_cmap(param.bp_colors, bit=False, log=log_file)

    # Mark the earthquake location.
    xpt, ypt = bm(eq_lon, eq_lat)
    bm.plot([xpt], [ypt], marker=earthquakes['marker'], markerfacecolor=earthquakes['color'], markeredgecolor='white',
            markersize=15, label='event')

    # Create a pseudocolor plot.
    bm.pcolormesh(grid_x, grid_y, grid_v, cmap=bp_cmap, alpha=0.7, shading='auto', linewidths=0,
                  vmin=0, vmax=1)

    # Avoid areas without coastlines.
    try:
        bm.drawcoastlines(color=param.fill_continents_color)
    except Exception as ex:
        if not coastline_skipped:
            print_message('WARN', f'Skipped drawcoastlines:\n{ex}', flush=True, log=log_file)
            coastline_skipped = True
        pass

    if basemap_countries:
        bm.drawcountries(color=param.fill_continents_color)
        if basemap_states:
            bm.drawstates(color=param.fill_continents_color)

    # labels = [left,right,top,bottom].
    bm.drawparallels(np.arange(int(lat_min), int(lat_max), 1), labels=[1, 0, 0, 0], fontsize=font_size[media]['label'],
                     linewidth=0.0,
                     labelstyle='+/-', fmt='%0.0f')
    bm.drawmeridians(np.arange(int(lon_min), int(lon_max), 2), labels=[0, 0, 0, 1], rotation=45,
                     fontsize=font_size[media]['label'],
                     linewidth=0.0,
                     labelstyle='+/-', fmt='%0.0f')

    trench_x, trench_y = lib.read_global_trenches(bmap=bm)
    bm.plot(trench_x, trench_y, color='black', linestyle=param.trench_linestyle, linewidth=0.5)

    # plt.ylabel('Latitude', labelpad=param.ylabel_pad, fontsize=font_size[media]['label'])

    title = f'Back Projection cumulative stack\n{param.trace_filter[param.bp_filter[vn_name]]["label"]}'
    plt.title(title, fontsize=font_size[media]['title'])

    # 2. Top Right Subplot, event and station locations.
    print_message('INFO', f'Summary plot /Station', log=log_file)
    ax2 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                           (0, 1), rowspan=tall_to_short_height, colspan=1)

    width = param.map_width
    sta_map = Basemap(width=width, height=width, projection=basemap_projection,
                      lat_0=lat_0, lon_0=lon_0)

    for dist_index, dist in enumerate(distance_circles):
        # Retrieve X and Y radius values using earthquake location as center point, draw dist degrees out.
        lon_1, lat_1 = lib.create_circle(lat_0, lon_0, dist)
        X, Y = sta_map(lon_1, lat_1)
        sta_map.plot(X, Y, marker=None, color='black', linestyle='--', linewidth=1)

        label_lat, label_lon = lib.get_location(lat_0, lon_0, 180.0, dist)
        x_l, y_l = sta_map(label_lon, label_lat)
        ax2.annotate(distance_circle_labels[dist_index], (x_l, y_l), xytext=(0, 0),
                     textcoords='offset points', c='black', backgroundcolor='white')

    sta_x = list()
    sta_y = list()
    print_message('INFO', f'Total {len(trace_list.keys())} stations selected', log=log_file)

    for sta_key in trace_list.keys():
        x, y = sta_map(trace_list[sta_key]['lon'], trace_list[sta_key]['lat'])
        sta_x.append(x)
        sta_y.append(y)

    if sta_x:
        # Plot stations.pt
        sta_map.plot(sta_x, sta_y, marker=vn[vn_name]['marker'],
                     markerfacecolor=vn[vn_name]['color'], linestyle='None',
                     markeredgewidth=0.3, markeredgecolor='None', alpha=1.0, markersize=6, zorder=3,
                     label=f'{vn_name}')

    # Avoid areas without coastlines.
    try:
        sta_map.drawcoastlines(linewidth=0.5, color=param.sta_map_color)
    except Exception as ex:
        if not coastline_skipped:
            print_message('WARN', f'Skipped drawcoastlines:\n{ex}', flush=True, log=log_file)
            coastline_skipped = True
        pass

    # Place a star at the center.
    xpt, ypt = sta_map(lon_0, lat_0)
    sta_map.plot([xpt], [ypt], marker=earthquakes['marker'], markerfacecolor=earthquakes['color'],
                 markeredgecolor='white', markersize=15,
                 label='event')
    x0, y0 = sta_map(lon_0, lat_0)
    lat_1, lon_1 = lib.get_2nd_point(lat_0, lon_0, eq_min_radius)
    if vn_name not in param.vn_azimuth_exception_list:
        lat_1, lon_1 = lib.get_2nd_point(lat_0, lon_0, eq_max_radius, bearing=center_azim - vn_azimuth)
        x1, y1 = sta_map(lon_1, lat_1)
        sta_map.plot((x0, x1), (y0, y1), '--', color=earthquakes['color'], lw=1)

        lat_2, lon_2 = lib.get_2nd_point(lat_0, lon_0, eq_max_radius, bearing=center_azim + vn_azimuth)
        x2, y2 = sta_map(lon_2, lat_2)
        sta_map.plot((x0, x2), (y0, y2), '--', color=earthquakes['color'], lw=1)

    # Network's name on the upper left
    plt.text(0.1, 0.95, vn_name, fontsize=font_size[media]['network'], horizontalalignment='center',
             verticalalignment='center', transform=ax2.transAxes, color='red',
             path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    title = f'{eq_date_time.replace("T", " ")}\nM{eq_magnitude} Z={eq_depth}km'
    plt.title(title, fontsize=font_size[media]['title'])

    # Save just the portion _inside_ the first axis's boundaries for the animation screen.
    #extent = full_extent(ax7, pad=0.25).transformed(fig.dpi_scale_trans.inverted())
    #anim_tag = 'BP'
    #file_tag = '_'.join([anim_tag, vn_name, lib.file_name_tag(eq_datetime)])
    #plot_file_name = f'{file_tag}_screen.png'
    #plot_file_name = os.path.join(param.video_dir, plot_file_name)
    #fig.savefig(plot_file_name, bbox_inches=extent)

    # 3. Middle left Subplot, local maxima.
    print_message('INFO', f'Summary plot / Local Maxima', log=log_file)
    ax3 = fig.add_subplot(323)
    ax3 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                           (3, 0), rowspan=tall_to_short_height, colspan=1)

    width = lat_max - lat_min
    width = degrees2kilometers(width) * 1000.0
    lon_0 = eq_lon
    lat_0 = eq_lat
    bm = Basemap(width=width, height=width, projection=basemap_projection,
                 lat_0=lat_0, lon_0=lon_0, resolution=param.basemap_resolution)

    coast_alpha = 1
    if param.fill_continents:
        coast_alpha = param.fill_continents_alpha
        bm.fillcontinents(color=param.fill_continents_color, alpha=coast_alpha)

    scalebar_deg = kilometer2degrees(param.scalebar_km)
    _lat0 = lat_min + (lat_max - lat_min) / 4.0
    _lon0 = lon_max - (lon_max - lon_min) / 10.0
    _lat1, _lon1 = lib.get_location(_lat0, _lon0, 0, scalebar_deg)
    _x, _y = bm((_lon0, _lon1), (_lat0, _lat1))
    print_message('INFO', f'Map scale between: ({_lat0:0.2f}, {_lon0:0.2f}) and ({_lat1:0.2f}, {_lon1:0.2f}) / '
                          f'({_x[0]:0.2f}, {_y[0]:0.2f}) and ({_x[1]:0.2f}, {_y[1]:0.2f})', log=log_file)
    # Use the same _X to ensure scale is vertical.
    bm.plot([_x[0], _x[0]], _y, color='black', linewidth=2)
    plt.text(_x[0], _y[1], f'  {param.scalebar_km} km', fontsize=font_size[media]['legend'],
             horizontalalignment='center', rotation=90,
             verticalalignment='bottom',
             path_effects=[pe.withStroke(linewidth=2, foreground='white')])

    # Make the custom color map.
    bp_cmap = lib.make_cmap(param.bp_colors, bit=False, log=log_file)

    # Mark the earthquake location.
    xpt, ypt = bm(eq_lon, eq_lat)
    zorder = 10000
    bm.plot([xpt], [ypt], marker=earthquakes['marker'], markerfacecolor=earthquakes['color'], markeredgecolor='white',
            markersize=15, label='USGS epicenter', linestyle='None', zorder=zorder)
    bm.plot(trench_x, trench_y, color='black', linestyle=param.trench_linestyle, linewidth=0.5, label='Trench')

    # Avoid areas without coastlines.
    try:
        bm.drawcoastlines(color=param.fill_continents_color)
    except Exception as ex:
        if not coastline_skipped:
            print_message('WARN', f'Skipped drawcoastlines:\n{ex}', flush=True, log=log_file)
            coastline_skipped = True
        pass

    if basemap_countries:
        bm.drawcountries(color=param.fill_continents_color)
        if basemap_states:
            bm.drawstates(color=param.fill_continents_color)

    # Post a note about the dot colors.
    plt.annotate(f'dot: local maxima\ncolor: time {down_arrow}', (1.0, .0), xycoords='axes fraction',
                 horizontalalignment='right', xytext=(-10, 10), fontsize=font_size[media]['legend'],
                 textcoords='offset points')

    # labels = [left,right,top,bottom].
    bm.drawparallels(np.arange(int(lat_min), int(lat_max), 1), labels=[1, 0, 0, 0], fontsize=font_size[media]['label'],
                     linewidth=0.0,
                     labelstyle='+/-', fmt='%0.0f')

    time_keys = list(stack_max.keys())
    times = list(np.array(time_keys, dtype=float))
    values = list(stack_max.values())

    # Find the local maxima and plot these within 50% of the maximum.
    # Note that the return value is a tuple even when data is 1-D.
    maxima_list = argrelextrema(np.array(values), np.greater)
    maxima_list = list(maxima_list[0])
    maxima_max = max(values)
    maxima_base = maxima_max * param.peak_marker_base_factor

    # Before we finish, check to see if this network has too many significant peaks.
    _count = 0
    for _max in maxima_list:
        if _max >= maxima_max * param.qc_max_peak_factor:
            _count += 1
    if _count > param.qc_max_peak_count:
        print_message('ERR', f'{vn_name} virtual network has {_count} '
                             f'peaks above {param. qc_max_peak_count} max. peak of {maxima_max}'
                             f' > {param.qc_max_peak_count}', log=log_file)
        make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None)
        make_insufficient_data_animation(eq_lat, eq_lon, eq_magnitude, None, anim_tag='syn')
        sys.exit(4)
    else:
        print_message('INFO', f'{vn_name} virtual network has {_count} '
                             f'peaks above {param. qc_max_peak_count} max. peak of {maxima_max}'
                             f' <= {param.qc_max_peak_count}', log=log_file)

    print_message('INFO', f'Peak maxima index list: {maxima_list}, maxima base: {maxima_base}', log=log_file)
    x = list()
    y = list()
    s = list()
    c = list()
    z = list()
    count = -1
    for max_index in maxima_list:
        if values[max_index] < maxima_base:
            continue
        t_index = time_keys[max_index]
        _x, _y = bm(stack_max_loc[t_index][1], stack_max_loc[t_index][0])
        x.append(_x)
        y.append(_y)
        val = values[max_index] / maxima_max
        count += 1
        c.append(count)
        # Scale marker area based on the amplitude size.
        marker_size = param.peak_marker_size_max * val * val
        s.append(marker_size)
        z.append(zorder - int(1000 * val))
    # Make sure the largest markers (lowest z) are plotted first.
    x = np.array(x)
    y = np.array(y)
    c = np.array(c)
    s = np.array(s)
    z = np.array(z)
    order = np.argsort(z)
    norm = colors.Normalize(vmin=min(c), vmax=max(c))
    if param.time_colors:
        c = np.array(param.time_colors[0:len(order)])
        bm.scatter(x[order], y[order], c=c[order], s=s[order], alpha=param.time_alpha, marker='o', linewidths=0,
                   edgecolor='black', linewidth=1, norm=norm, zorder=1000)
    else:
        bm.scatter(x[order], y[order], cmap=param.time_cmap, c=c[order], alpha=param.time_alpha,
                   s=s[order], marker='o', linewidths=0, edgecolor='black', linewidth=1,
                   norm=norm, zorder=1000)

    legend = plt.legend(loc='lower left', framealpha=1.0, frameon=False, facecolor='white', edgecolor=None,
                        fontsize=font_size[media]['legend'])

    # plt.ylabel('Latitude', labelpad=param.ylabel_pad, fontsize=font_size[media]['label'])

    # 4. Middle Right Subplot, traces.
    ax4 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                           (3, 1), rowspan=tall_to_short_height, colspan=1)

    # Plot normalized and  phase-aligned traces.
    # Unfiltered data.
    stacked = None
    stacked_filtered = None
    trace_count = len(trace_list.keys())
    for tr_key in trace_list.keys():

        _tr = trace_list[tr_key]['tr_final'].copy()
        _tr_filtered = trace_list[tr_key]['tr_filter'].copy()

        _tr_times = _tr.times(reftime=eq_datetime)

        # Stack with no weight applied.
        if stacked is None:
            stacked = _tr.data.copy()
            stacked_filtered = _tr_filtered.data.copy()
        else:
            # This is to address the issue that in rare cases, there is one sample missing.
            # This is a partial solution for now.
            if len(stacked) < len(_tr.data):
                _tr_copy = _tr.data.copy()
                _tr_filtered_copy = _tr_filtered.data.copy()
                stacked = np.add(stacked, _tr_copy[:len(stacked)])
                stacked_filtered = np.add(stacked_filtered, _tr_filtered_copy[:len(stacked_filtered)])
            else:
                stacked = np.add(stacked, _tr.data.copy())
                stacked_filtered = np.add(stacked_filtered, _tr_filtered.data.copy())
        # Unstacked raw traces.
        ax4.plot(_tr_times, _tr.data / trace_list[tr_key]['weight'] + 4.0, "k-", lw=0.4, clip_on=True, zorder=2000)

        # Unstacked filtered traces.
        ax4.plot(_tr_times, _tr_filtered.data / trace_list[tr_key]['weight'] + 0.0, "k-", lw=0.4,
                 clip_on=True, zorder=2000)

    trace_end = param.max_summary_trace_time

    print_message('INFO', f'Stack time range from {stack_start}s to {trace_end}s', log=log_file)
    ax4.set_xlim(param.min_summary_trace_time, param.max_summary_trace_time)
    ax4.set_ylim(-1, 7)
    plt.tick_params(left=False)
    ax4.set_xlabel('Time (second)')
    ax4.axes.yaxis.set_ticklabels([])

    label_y_position = param.label_y_position

    title = param.prefilter_label

    # Unstacked raw traces' title.
    t1 = plt.text(trace_end - 5, label_y_position + 4, title, fontsize=font_size[media]['label'],
                  horizontalalignment='right', verticalalignment='center', color='black',
                  path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=2100)

    # Stacked raw traces' title.
    title = f'stacked - {param.prefilter_label}'
    t2 = plt.text(trace_end - 5, label_y_position + 6, title, fontsize=font_size[media]['label'],
                  horizontalalignment='right', verticalalignment='center', color='red',
                  path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=1100)

    stacked_max = np.max(np.abs(stacked))
    stacked /= stacked_max

    # Stacked un filtered trace.
    ax4.plot(_tr_times, stacked + 6.0, "r-", lw=0.4, clip_on=True, zorder=1000)

    stacked_max = np.max(np.abs(stacked_filtered))
    stacked_filtered /= stacked_max

    # Stacked filtered trace.
    ax4.plot(_tr_times, stacked_filtered + 2.0, "r-", lw=0.4, clip_on=True, zorder=1000)

    # Unstacked traces title.
    title = f'filtered {param.trace_filter[param.bp_filter[vn_name]]["label"]}'
    t3 = plt.text(trace_end - 5, label_y_position, title, fontsize=font_size[media]['label'],
                  horizontalalignment='right', verticalalignment='center', color='black',
                  path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=2100)

    # Stacked trace title.
    title = f'stacked - {param.trace_filter[param.bp_filter[vn_name]]["label"]}'
    t4 = plt.text(trace_end - 5, label_y_position + 2, title, fontsize=font_size[media]['label'],
                  horizontalalignment='right', verticalalignment='center', color='red',
                  path_effects=[pe.withStroke(linewidth=2, foreground='white')], zorder=1100)

    title = 'Time-shifted and aligned'
    plt.title(title, fontsize=font_size[media]['title'])

    # Get info on this subplot so we can set its height the same as the map to the lef.
    # We always want to adopt the map's height since it is dynamic.
    map_bb = ax3.get_position()
    trace_bb = ax4.get_position()
    ax4.set_position([trace_bb.x0, trace_bb.y0, trace_bb.width, map_bb.height])

    # 5. Local maxima time.
    ax5 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                           (6, 0), rowspan=1, colspan=1)
    times = list(stack_max.keys())
    times = list(np.array(times, dtype=float))
    values = list(stack_max.values())
    ax5.fill_between(times, 0, values, facecolor=param.beam_power_fill_color)

    x = list()
    y = list()
    c = list()
    s = list()

    count = -1
    for max_index in maxima_list:
        if values[max_index] < maxima_base:
            continue
        _t = times[max_index]
        x.append(_t)
        y.append(values[max_index])
        count += 1
        c.append(count)
        s.append(param.peak_marker_size)

    # Set the color map range.
    if param.time_colors:
        c = np.array(param.time_colors[0:len(x)])
        plt.scatter(x, y, c=c, s=s, marker='o', linewidths=0, alpha=param.time_alpha, edgecolor='black', linewidth=1,
                    norm=norm, zorder=1000)
    else:
        plt.scatter(x, y, cmap=param.time_cmap, c=c, alpha=param.time_alpha, s=s, marker='o', linewidths=0,
                    edgecolor='black', linewidth=1, norm=norm, zorder=1000)

    ax5.axes.yaxis.set_ticklabels([])

    # Title for the subplot above and this subplot.
    title = f'{up_arrow} location map / {down_arrow} time plot\nPeak BP stack amplitudes'
    plt.title(title, fontsize=font_size[media]['label'])
    ax5.text(0.05, 0.2, f'dot: local maxima; color: time', horizontalalignment='left',
             fontsize=font_size[media]['legend'],
             verticalalignment='top', transform=ax5.transAxes)

    ax5.set_xlim(stack_start, stack_end)
    ax5.set_ylim(bottom=0.0)
    ax5.yaxis.set_ticks_position('none')
    ax5.set_xlabel('time relative to origin (sec.)')

    # 6. Lover right, logo, etc..
    ax6 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                           (6, 1), rowspan=1, colspan=1)
    plt.axis('off')

    # Plot the logo.
    delimiter = ' '
    if os.path.isfile(logo_image):
        logo_image = np.array(Image.open(logo_image))
        im = OffsetImage(logo_image, zoom=logo_zoom, alpha=logo_alpha, zorder=1000)
        b_box = AnnotationBbox(im, xy=logo_location, xycoords=logo_coords, box_alignment=(0.0, 0.0),
                               boxcoords='offset pixels', frameon=False)

        # Add the AnnotationBbox artist.
        ax6.add_artist(b_box)
        delimiter = '\n'

    # Production date time stamp.
    production_date = lib.version_timestamp(script_version, search_radius_edge, delimiter=delimiter)
    ax6.annotate(production_date, xy=logo_location, xycoords='axes pixels',
                 xytext=(logo_location[0] + logo_width + logo_padding, logo_location[1] + logo_height / 2.0),
                 textcoords='offset pixels', horizontalalignment='left', fontsize=production_label_font_size,
                 verticalalignment='center')

    # Place the tight_layout after most of the elements are in, so the layout can be configured properly.
    # Info - pad: Padding between the figure edge and the edges of subplots, as a fraction of the font size.
    # Info - h_pad, w_pad : Padding (height/width) between edges of adjacent subplots, as a fraction of the font size.
    #        Defaults to pad.
    plt.text(x=0.5, y=0.6, s=f'{vn_name} Virtual Network ({trace_count} stations)\nBack Projection Summary Plot',
             fontsize=font_size[media]['title'], ha='center')

    plt.tight_layout(pad=2, h_pad=0, w_pad=0, rect=None)

    # Save the figure to a .png file
    _tag = 'BP_summary'
    file_tag = '_'.join([_tag, vn_name, lib.file_name_tag(eq_datetime)])
    plot_file_name = f'{file_tag}.png'
    plot_file_name = os.path.join(param.image_dir, plot_file_name)
    plt.savefig(plot_file_name, bbox_inches='tight', dpi=param.dpi, facecolor='white')
    print_message('INFO', f'Summary plot for the {vn_name} virtual network saved as {plot_file_name}', log=log_file)
    plt.close()

    # Generate a screen image for animations.
    if create_animation:
        for bp_plot_type in ('BP_screen', 'BP_syn_screen'):
            # ax7 Plot screen plots for the video.
            # Video frame layout.
            subplot_columns = 1
            subplot_tall_rows = 1
            subplot_short_rows = 1
            tall_to_short_height = 3

            fig = plt.figure(figsize=param.video_size, facecolor='white')
            ax7 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                                   (0, 0), rowspan=tall_to_short_height, colspan=1)

            # Extract grid and values.
            lat = list()
            lon = list()
            val = list()

            if bp_plot_type == 'BP_syn_screen':
                # Compute a cumulative stack.
                cumulative_stack = dict()
                cumulative_stack_max = None
                for grid_key in syn_stack_final:
                    if grid_key not in cumulative_stack:
                        cumulative_stack[grid_key] = 0.0
                        for time_key in syn_stack_final[grid_key]:
                            cumulative_stack[grid_key] = np.add(cumulative_stack[grid_key],
                                                                syn_stack_final[grid_key][time_key])
                        if cumulative_stack_max is None:
                            cumulative_stack_max = abs(cumulative_stack[grid_key])
                        else:
                            cumulative_stack_max = max(cumulative_stack_max, abs(cumulative_stack[grid_key]))

            for grid_key in cumulative_stack:
                    _lat, _lon = grid_key.split('_')
                    lat.append(float(_lat))
                    lon.append(float(_lon))
                    value = cumulative_stack[grid_key] / cumulative_stack_max
                    val.append(value)

            # Must be np arrays for grid
            lat_list = np.array(lat)
            lon_list = np.array(lon)
            value_list = np.array(val)

            # Find the min and max of coordinates.
            lon_min = lon_list.min()
            lat_min = lat_list.min()
            lon_max = lon_list.max()
            lat_max = lat_list.max()

            latitude, longitude = lib.set_grid(eq_lat, eq_lon, eq_magnitude, grid_factor=grid_factor)

            # Now let's grid the data. Find the number of grid points in each direction.
            lon_num = pcolormesh_grid_factor * int(((lon_max - lon_min) / longitude['inc']) + 1)
            lat_num = pcolormesh_grid_factor * int(((lat_max - lat_min) / latitude['inc']) + 1)

            width = lat_max - lat_min
            width = degrees2kilometers(width) * 1000.0
            lon_0 = eq_lon
            lat_0 = eq_lat
            bm = Basemap(width=width, height=width, projection=basemap_projection,
                         lat_0=lat_0, lon_0=lon_0, resolution=param.basemap_resolution)

            coast_alpha = 1
            if param.fill_continents:
                coast_alpha = param.fill_continents_alpha
                bm.fillcontinents(color=param.fill_continents_color, alpha=coast_alpha)

            # Network's name on the upper left
            plt.text(0.1, 0.95, vn_name, fontsize=font_size[media]['network'], horizontalalignment='center',
                     verticalalignment='center', transform=ax7.transAxes, color='red',
                     path_effects=[pe.withStroke(linewidth=2, foreground='white')])

            scalebar_deg = kilometer2degrees(param.scalebar_km)
            _lat0 = lat_min + (lat_max - lat_min) / 4.0
            _lon0 = lon_max - (lon_max - lon_min) / 10.0
            _lat1, _lon1 = lib.get_location(_lat0, _lon0, 0, scalebar_deg)
            _x, _y = bm((_lon0, _lon0), (_lat0, _lat1))
            print_message('INFO', f'Map scale between: ({_lat0:0.2f}, {_lon0:0.2f}) and ({_lat1:0.2f}, {_lon1:0.2f}) / '
                                  f'({_x[0]:0.2f}, {_y[0]:0.2f}) and ({_x[1]:0.2f}, {_y[1]:0.2f})', log=log_file)
            # Use the same _X to ensure scale is vertical.
            bm.plot([_x[0], _x[0]], _y, color='black', linewidth=2)
            plt.text(_x[0], _y[1], f'  {param.scalebar_km} km', fontsize=font_size[media]['legend'],
                     horizontalalignment='center', rotation=90,
                     verticalalignment='bottom',
                     path_effects=[pe.withStroke(linewidth=2, foreground='white')])

            # Create a uniform mesh for contouring. First transfer lon, lat to map units (x, y).
            x_old, y_old = bm(lon_list, lat_list)
            x_new = np.linspace(min(x_old), max(x_old), lon_num)
            y_new = np.linspace(min(y_old), max(y_old), lat_num)

            # Basic mesh in map's x, y.
            grid_x, grid_y = np.meshgrid(x_new, y_new)

            try:
                # Interpolate at the points in lon_new, lat_new.
                # Method : {'linear', 'nearest', 'cubic'}, optional.
                grid_v = griddata((x_old, y_old), value_list, (grid_x, grid_y), method='cubic')
            except Exception as _er:
                print_message('ERR', f'Griding failed: {_er}', log=log_file)
                sys.exit(2)

            # Make the custom color map.
            bp_cmap = lib.make_cmap(param.bp_colors, bit=False, log=log_file)

            # Mark the earthquake location.
            xpt, ypt = bm(eq_lon, eq_lat)
            bm.plot([xpt], [ypt], marker=earthquakes['marker'], markerfacecolor=earthquakes['color'],
                    markeredgecolor='white',
                    markersize=15, label='event')

            # Create a pseudocolor plot.
            bm.pcolormesh(grid_x, grid_y, grid_v, cmap=bp_cmap, alpha=0.7, shading='auto', linewidths=0,
                          vmin=0, vmax=1)

            # Avoid areas without coastlines.
            try:
                bm.drawcoastlines(color=param.fill_continents_color)
            except Exception as ex:
                if not coastline_skipped:
                    print_message('WARN', f'Skipped drawcoastlines:\n{ex}', flush=True, log=log_file)
                    coastline_skipped = True
                pass

            if basemap_countries:
                bm.drawcountries(color=param.fill_continents_color)
                if basemap_states:
                    bm.drawstates(color=param.fill_continents_color)

            # labels = [left,right,top,bottom].
            bm.drawparallels(np.arange(int(lat_min), int(lat_max), 1), labels=[1, 0, 0, 0],
                             fontsize=font_size[media]['label'],
                             linewidth=0.0,
                             labelstyle='+/-', fmt='%0.0f')
            bm.drawmeridians(np.arange(int(lon_min), int(lon_max), 2), labels=[0, 0, 0, 1], rotation=0,
                             fontsize=font_size[media]['label'],
                             linewidth=0.0,
                             labelstyle='+/-', fmt='%0.0f')

            trench_x, trench_y = lib.read_global_trenches(bmap=bm)
            bm.plot(trench_x, trench_y, color='black', linestyle=param.trench_linestyle, linewidth=0.5)

            # plt.ylabel('Latitude', labelpad=param.ylabel_pad, fontsize=font_size[media]['label'])

            if bp_plot_type == 'BP_screen':
                title = (f"{eq_datetime.strftime('%Y-%m-%d %H:%M:%S')} M{eq_magnitude} Z={eq_depth}km\n"
                         f"{param.trace_filter[param.bp_filter[vn_name]]['label']}")
            else:
                title = (f"{eq_datetime.strftime('%Y-%m-%d %H:%M:%S')} M{eq_magnitude} Z={eq_depth}km\n"
                         f'Array Response Function (synthetics) {param.trace_filter[param.bp_filter[vn_name]]["label"]}')
            plt.title(title, fontsize=font_size[media]['title'])

            # Get info on this subplot so we can align the one below it.
            map_bbox = ax7.get_position()

            # Plot the beam power.
            ax8 = plt.subplot2grid((subplot_tall_rows * tall_to_short_height + subplot_short_rows, subplot_columns),
                                   (3, 0), rowspan=1, colspan=1)

            if bp_plot_type == 'BP_syn_screen':
                times = list(syn_stack_max.keys())
                times = list(np.array(times, dtype=float))
                values = list(syn_stack_max.values())
            else:
                times = list(stack_max.keys())
                times = list(np.array(times, dtype=float))
                values = list(stack_max.values())
            max_value = max(values)
            values = np.array(values) / max_value
            values = list(values)
            ax8.fill_between(times, 0, values, facecolor=param.beam_power_fill_color)
            ax8.set_xlim(stack_start, stack_end)
            ax8.set_ylim(bottom=0.0)
            ax8.set_xlabel('time relative to origin (sec.)')
            ax8.set_ylabel('beam power')
            ax8.yaxis.set_ticklabels([])

            # Get info on this subplot so we can align the one above it.
            # We always want to adopt the map's width since it is dynamic.
            power_bbox = ax8.get_position()
            ax8.set_position([map_bbox.x0, power_bbox.y0, map_bbox.width, power_bbox.height])

            # Save the figure to a .png file
            _tag = bp_plot_type
            file_tag = '_'.join([_tag, vn_name, lib.file_name_tag(eq_datetime)])
            plot_file_name = f'{file_tag}.png'
            plot_file_name = os.path.join(param.video_dir, plot_file_name)
            plt.savefig(plot_file_name, bbox_inches='tight', dpi=param.dpi, facecolor='white')
            print_message('INFO', f'Screen plot for the {vn_name} virtual network saved as {plot_file_name}',
                          log=log_file)
            plt.close()


