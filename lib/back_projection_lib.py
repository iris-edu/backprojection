import sys
import os

import math

from obspy.geodetics import degrees2kilometers
# NOTE: gps2dist_azimuth will check if you have installed the Python module geographiclib. It will be better to
#       have it installed.
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy import UTCDateTime
from obspy import read
from obspy.taup import TauPyModel
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.signal.filter import envelope

from obspy.core.stream import Stream

import matplotlib as mpl

from datetime import datetime

import numpy as np
import time

from urllib.request import urlopen

# Import the back projection parameters and libraries.
_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
param_dir = os.path.join(_dir, 'param')

sys.path.append(param_dir)

import back_projection_param as param

"""
    Description:

    A Python utility library used by the main script.

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

# Parameters.
vn = param.virtual_networks
vn_list = list(vn.keys())
vn_min_radius = param.vn_min_radius
vn_max_radius = param.vn_max_radius
vn_azimuth = param.vn_azimuth

model = TauPyModel(model=param.travel_time_model)

channel_order = param.channel_order
channel_list = channel_order.keys()

eq_min_radius = param.eq_min_radius
eq_max_radius = param.eq_max_radius

dc_to_exclude = param.dc_to_exclude

chunk_count = param.chunk_count

earthquakes = param.earthquakes

fedcatalog_service_url = param.fedcatalog_service_url

log_file = sys.stdout
verbose = param.verbose

timing = param.timing


class ObjDict(dict):
    """Accessing dictionary items as object attributes:
            https://goodcode.io/articles/python-dict-object/
    """

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: {}".format(name))

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: {}".format(name))


def version_timestamp(version, search_radius, delimiter=' '):
    current_time = datetime.utcnow()
    if search_radius is None:
        second_line = f'Chan: {param.request_channel}'
    else:
        second_line = f'Chan: {param.request_channel}, sparse search: {search_radius}°'

    timestamp = f'{param.production_label}{delimiter}({version}/{current_time.strftime("%Y-%m-%d %H:%M")} UTC)'\
                f'\n{second_line}'
    return timestamp


def file_name_tag(event_date, date=None):
    """Create a file name tag from the event's date and time."""
    if date is not None:
        event_date = datetime.strptime(event_date, '%Y-%m-%dT%H:%M:%S')
    tag = f"{event_date.strftime('%Y.%m.%d.%H.%M')}"
    return tag


def print_message(flag, text, flush=True, log=sys.stdout, end='\n'):
    """Print out a message. Force the flush and write to the file handle"""

    if flag == 'ERR':
        print(f'\n\n{60 * "="}\n', file=log, flush=flush, end=end)
        print(f'[{flag}] {text}', sep='\n', file=log, flush=flush, end=end)
        print(f'\n\n{60 * "="}\n', file=log, flush=flush, end=end)
    elif log is not None:
        print(f'[{flag}] {text}', file=log, flush=flush, end=end)
    else:
        print(f'[{flag}] {text}', flush=flush, end=end)


def read_global_trenches(bmap=None, log=sys.stdout):
    """Read location of the global trenches from the data file."""
    global_trench_file = os.path.join(param.assets_dir, param.global_trenches_file)
    try:
        fp = open(global_trench_file, 'r')
    except Exception as ex:
        print_message('ERR', f'problem reading the trench file {global_trench_file}, will not plot trenches\n{ex}', log=log)
        return None, None

    data = fp.read()
    fp.close()
    lines = data.split('\n')
    trench_x = list()
    trench_y = list()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        values = line.split()
        if values[0] == 'NaN':
            trench_x.append(float('NaN'))
            trench_y.append(float('NaN'))
        elif map is None:
            trench_x.append(float(values[0]))
            trench_y.append(float(values[1]))
        else:
            _x, _y = bmap(float(values[0]), float(values[1]))
            trench_x.append(_x)
            trench_y.append(_y)

    return trench_x, trench_y


def get_location(lat1, lon1, brng, distance_degrees):
    """find location of a point at a given distance from lat0, lon0
    source: https://stochasticcoder.com/2016/04/06/python-custom-distance-radius-with-basemap/"""
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    distance_radians = math.radians(distance_degrees)
    brng = math.radians(brng)

    lat2 = math.asin(math.sin(lat1) * math.cos(distance_radians)
                     + math.cos(lat1) * math.sin(distance_radians) * math.cos(brng))

    lon2 = lon1 + math.atan2(math.sin(brng) * math.sin(distance_radians)
                             * math.cos(lat1), math.cos(distance_radians) - math.sin(lat1) * math.sin(lat2))

    lon2 = math.degrees(lon2)
    lat2 = math.degrees(lat2)

    return lat2, lon2


def create_circle(lat0, lon0, radius_degrees):
    """Create a circle around lat0, lon0 for a given radius
    source: https://stochasticcoder.com/2016/04/06/python-custom-distance-radius-with-basemap/"""
    lat_list = list()
    lon_list = list()

    for brng in range(0, 360):
        lat2, lon2 = get_location(lat0, lon0, brng, radius_degrees)
        lat_list.append(lat2)
        lon_list.append(lon2)

    return lon_list, lat_list


def make_cmap(colors, position=None, bit=False, log=sys.stdout):
    """
    Source: Chris Slocum (with some modifications)  http://schubert.atmos.colostate.edu/~cslocum/custom_cmap.html

    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    """

    bit_rgb = np.linspace(0, 1, 256)
    if position is None:
        position = np.linspace(0, 1, len(colors))
    else:
        if len(position) != len(colors):
            print_message('ERR', f'position length ({len(position)}) '
                                 f'must be the same as colors ({len(colors)})', log=log)
            sys.exit(2)
        elif position[0] != 0 or position[-1] != 1:
            print_message('ERR', f'position must start with 0 and end with 1 ({position[0]},'
                                 f' {position[-1]})', log=log)
            sys.exit(2)
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red': list(), 'green': list(), 'blue': list()}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('bp_colormap', cdict, 256)
    return cmap


def set_float_digits(value, places=4):
    """Set a float number precision"""
    form = f'0.{places}f'
    new_value = float(f'{value:{form}}')

    return new_value


def case(this_val, case_dict):
    """Simulates a case statement to check a value against given ranges in a dictionary."""

    default = case_dict['default']
    condition = case_dict['condition']
    ranges = case_dict['ranges']

    reverse = False
    if condition in ('>', '>='):
        reverse = True
    value_items = sorted(ranges.items(), reverse=reverse)

    for limit, value in value_items:
        if condition == '<' and this_val < limit:
            return value
        elif condition == '<=' and this_val <= limit:
            return value
        elif condition == '==' and this_val == limit:
            return value
        elif condition == '>' and this_val > limit:
            return value
        if condition == '>=' and this_val >= limit:
            return value
        elif condition == '!=' and this_val != limit:
            return value
    return default


def time_it(t, stage='', end='\n', log=sys.stdout):
    """Compute elapsed time since the last call."""
    if stage is None:
        return t
    t1 = time.time()
    dt = t1 - t
    if dt < param.timing_threshold:
        return t
    if dt < 1.0:
        print_message('TIME', f'{stage} in {dt:0.3f} second', end=end, log=log)
    elif dt < 5.0:
        print_message('TIME', f'{stage} in {dt:0.2f} seconds', end=end, log=log)
    elif dt < 100.0:
        print_message('TIME', f'{stage} in {dt:0.1f} seconds', end=end, log=log)
    else:
        print_message('TIME', f'{stage} in {dt:0.0f} seconds', end=end, log=log)
    t = t1
    return t


def mkdir(target_directory, log=sys.stderr):
    """ Make a directory if it does not exist."""
    directory = None
    try:
        directory = target_directory
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
    except Exception as _er:
        print_message('ERR', f'failed to create directory {directory}\n{_er}', log=log)
        return None


def sign(number):
    """Find sign of a number"""

    num_sign = math.copysign(1, number)
    return num_sign


def stfs_from_mccc_traces(eq_datetime, eq_magnitude, trace_list, optimal_mccc, double_mccc):
    """Create STFs based on the MCCCC traces."""

    # Window length before the phase.
    pre_phase_seconds = 20.0

    # STF length in seconds.
    stf_length = case(eq_magnitude, param.stf_length_dict)

    for net_sta in trace_list.keys():
        trace = trace_list[net_sta]

        p_time = trace['phase_delay']
        # xcorr_win is the xcorr window starting pre_phase_seconds before the phase.
        segment_1 = trace['tr_final'].slice(starttime=eq_datetime + p_time - pre_phase_seconds,
                                            endtime=eq_datetime + p_time - pre_phase_seconds + stf_length)
        segment_1.normalize()

        segment_2 = trace['tr_final'].slice(starttime=eq_datetime + p_time - pre_phase_seconds,
                                            endtime=eq_datetime + p_time - pre_phase_seconds +
                                            optimal_mccc[net_sta]['delay'])
        segment_2.normalize()

        segment_3 = trace['tr_final'].slice(starttime=eq_datetime + p_time - pre_phase_seconds,
                                            endtime=eq_datetime + p_time - pre_phase_seconds +
                                            double_mccc[net_sta]['delay'])
        segment_3.normalize()


def station_weight(trace_list, vn_name, intra_station_dist):
    """Compute wight for each active station based on the distance and azimuth of the neighboring stations."""

    sta_weight_dist = param.virtual_networks[vn_name]['sta_weight_dist']
    sta_weight_azim = param.virtual_networks[vn_name]['sta_weight_azim']
    weight = dict()
    # Check one station at a time.
    for center_key in trace_list.keys():
        center = trace_list[center_key]

        neighbor_dist_count = 0
        neighbor_azim_count = 0

        center_lat = center['lat']
        center_lon = center['lon']

        center_azim = center['azim']

        weight[center_key] = 1.0

        # Find how far are other stations from this station.
        for sta_net_key in trace_list.keys():
            if center_key == sta_net_key:
                continue

            neighbor = trace_list[sta_net_key]
            center_neighbor_key = f'{center_key}_{sta_net_key}'

            if center_neighbor_key in intra_station_dist:
                _dist, _distk,  _azim, _back_azim = intra_station_dist[center_neighbor_key]
            else:
                _dist, _azim, _back_azim = gps2dist_azimuth(center_lat, center_lon,
                                                            neighbor['lat'],
                                                            neighbor['lon'])
                _distk = _dist / 1000.0
                _dist = kilometer2degrees(_distk)
                intra_station_dist[f'{center_key}_{neighbor}'] = (_dist, _distk, _azim, _back_azim)
                intra_station_dist[f'{neighbor}_{center_key}'] = (_dist, _distk, _azim, _back_azim)

            if _distk < sta_weight_dist:
                neighbor_dist_count += 1

            azim_diff_1 = abs(center_azim - _azim)
            azim_diff_2 = abs((360 - center_azim) - _azim)
            azim_diff_3 = abs(center_azim - (_azim - 360))
            azim_diff = min(azim_diff_1, azim_diff_2, azim_diff_3)

            if azim_diff < sta_weight_azim:
                neighbor_azim_count += 1
        _weight = max(neighbor_dist_count, neighbor_azim_count, neighbor_dist_count * neighbor_azim_count)
        if _weight == 0:
            weight[center_key] = 1.0
        else:
            weight[center_key] = 1.0 / _weight

    return intra_station_dist, weight


def get_slice(eq_datetime, trace, pre_sec, total_sec, delay=0.0):
    """Get a slice of a trace."""

    # Calculate the Phase arrival time with correction_sec representing additional correction.
    p_time = eq_datetime + trace['phase_delay'] + delay
    _start = p_time - pre_sec
    _end = _start + total_sec

    # Should use filtered or raw traces.
    if param.filtered_trace_align:
        trace_slice = trace['tr_filter'].copy()
    else:
        trace_slice = trace['tr_final'].copy()

    # Now slice the trace.
    trace_slice.trim(starttime=_start, endtime=_end, pad=True, nearest_sample=True, fill_value=0.0)
    return trace_slice


def mccc(eq_datetime, trace_list, xcorr_win, shift_t, correction_sec=0.0):
    """Multi-channel cross-correlation to estimate time delays using multiple shift."""
    mccc_cc_max = dict()
    mccc_cc_shift_sec = dict()
    pre_phase_seconds = param.pre_phase_seconds

    # shift represents the number of samples to shift for cross correlation.
    shift = int(shift_t * param.trace_sampling_frequency)

    bad_traces = list()
    for net_sta_a in trace_list.keys():

        # xcorr_win is the xcorr window starting pre_phase_seconds before the phase.
        segment_a = get_slice(eq_datetime, trace_list[net_sta_a], pre_phase_seconds,
                              xcorr_win, correction_sec)

        cc_max_list = list()
        cc_max_index_list = list()

        for net_sta_b in trace_list.keys():
            if net_sta_a == net_sta_b:
                continue
            segment_b = get_slice(eq_datetime, trace_list[net_sta_b], pre_phase_seconds,
                                  xcorr_win, correction_sec)

            cc = correlate(segment_a, segment_b, shift, demean=False)

            # Use the absolute value to find the max.
            shift_index, shift_value = xcorr_max(cc, abs_max=True)

            # Keep index for time shift calculation and value for comarison.
            if not math.isnan(shift_value):
                cc_max_index_list.append(shift_index)
                cc_max_list.append(abs(shift_value))

        # MCCC for net_sta_a is done, so save the necessary values.
        if cc_max_list:
            # Calculate mean of the CC parameters.
            mean_cc_max = np.mean(cc_max_list)
            mean_cc_max_index = int(np.mean(cc_max_index_list))
            mccc_cc_shift_sec[net_sta_a] = mean_cc_max_index / param.trace_sampling_frequency
            mccc_cc_max[net_sta_a] = mean_cc_max

        else:
            # Station with bad mccc.
            bad_traces.append(net_sta_a)

    return bad_traces, mccc_cc_max, mccc_cc_shift_sec


def find_optimal_mccc_window(eq_datetime, trace_list, vnet, log=sys.stdout):
    """Performs MCCC for multiple windows."""
    # 1103.
    optimal_mccc_window = dict()
    stf_pre_phase_seconds = param.stf_pre_phase_seconds
    stf_post_phase_seconds = param.stf_post_phase_seconds
    stf_search_seconds = param.stf_search_seconds
    shift = param.xcorr_shift_default

    # Go through all xcorr windows and see which one produces the maximum stacking amplitude.
    for xc_window in param.xcorr_window:
        # Create empty stack stream.
        stack_stream = read()
        stack_stream.clear()

        # Compute CC and time shift for this window.
        print_message('INFO', f'Evaluating mccc of {len(trace_list)} traces for xc_window of {xc_window}s.', log=log)
        bad_traces, cc_max, cc_delay = mccc(eq_datetime, trace_list, xc_window, shift)

        rejected_list = bad_traces.copy()
        # Inspect CC for each station.
        for tr_key in cc_max.keys():
            if tr_key in bad_traces:
                continue

            # Only apply the correction if improvement is more than minimum. If this window is
            # selected, then the low CC stations will be removed. Segment is a little longer than
            # needed to avoid end effects.
            if cc_max[tr_key] >= vnet['xcorr_min']:
                _segment = get_slice(eq_datetime, trace_list[tr_key], param.pre_phase_seconds,
                                     xc_window + param.pre_phase_seconds, cc_delay[tr_key])

                stack_stream.append(_segment)
            else:
                rejected_list.append(tr_key)
                print_message('WARN', f'Skipped {tr_key}  for MCCC {xc_window} '
                                      f'because CC {abs(cc_max[tr_key]):0.2f} is < {vnet["xcorr_min"]}', log=log)
        # Did anything go into stack_stream?
        if not stack_stream:
            print_message('WARN', f'For xcorr_window of {xc_window} s. The stack is empty!', log=log)
            continue

        stack = np.sum([tr.data for tr in stack_stream], axis=0) / len(stack_stream)
        max_search_index = int(round(stf_search_seconds * param.trace_sampling_frequency))
        stf = stack[0:max_search_index].copy()
        stack_max = np.max(abs(stf))
        print_message('INFO', f'Evaluated xcorr_window of {xc_window} s. The stack maximum is {stack_max:0.2f} for '
                              f'{len(stack_stream)} traces.', log=log)

        if not optimal_mccc_window:
            optimal_mccc_window['window_length'] = xc_window
            optimal_mccc_window['cc_max'] = cc_max
            optimal_mccc_window['delay'] = cc_delay.copy()
            optimal_mccc_window['stack_max'] = stack_max
            optimal_mccc_window['rejected_list'] = rejected_list
            optimal_mccc_window['stf'] = stf
        elif stack_max > optimal_mccc_window['stack_max']:
            optimal_mccc_window['window_length'] = xc_window
            optimal_mccc_window['cc_max'] = cc_max
            optimal_mccc_window['delay'] = cc_delay.copy()
            optimal_mccc_window['stack_max'] = stack_max
            optimal_mccc_window['rejected_list'] = rejected_list
            optimal_mccc_window['stf'] = stf

    # Did we find optimal_mccc_window?
    if not optimal_mccc_window:
        return trace_list, len(stack_stream), None

    print_message('INFO', f'Optimum X-Corr shift is {optimal_mccc_window["window_length"]:0.1f} seconds', log=log)

    # Only select traces that contributed to the optimum x-corr, if their amplitude above the
    # minimum designated value. Remove traces with bad MCCCC.
    for key in optimal_mccc_window['rejected_list']:
        if key in optimal_mccc_window['cc_max'].keys():
            print_message('WARN', f'Removing {key} from trace list in "find_optimal_mccc_window" '
                                  f'due to bad MCCC', log=log)
            optimal_mccc_window['cc_max'].pop(key)
            optimal_mccc_window['delay'].pop(key)
        if key in trace_list:
            trace_list.pop(key)

    # Reset the stack stream.
    stack_stream = read()
    stack_stream.clear()

    xcorr = int(min(25.0, optimal_mccc_window['window_length']) * param.trace_sampling_frequency) + 1
    for tr_key in optimal_mccc_window['cc_max']:
        cc_delay = optimal_mccc_window['delay'][tr_key]
        trace_list[tr_key]['mccc_delay'] = cc_delay
        print_message('INFO', f'{tr_key} MCCC {optimal_mccc_window["window_length"]} '
                              f'CC {abs(cc_max[tr_key]):0.2f} is >= {vnet["xcorr_min"]}',
                      log=log)

        # Do one last pass at aligning by CCing w the best STF.
        #tr_slice = get_slice(eq_datetime, trace_list[tr_key], param.stf_pre_phase_seconds,
        #                     stf_post_phase_seconds + param.pre_phase_seconds, cc_delay)
        #cc = correlate(tr_slice, optimal_mccc_window['stf'], shift)
        #index_shift, value = xcorr_max(cc, abs_max=True)
        #correction = index_shift / param.trace_sampling_frequency
        #print_message('INFO', f'{tr_key} additional {correction} seconds correction from STF', log=log)
        #optimal_mccc_window['mccc_delay'][tr_key] += correction

        # Apply the correction and stack.
        tr_slice = get_slice(eq_datetime, trace_list[tr_key], stf_pre_phase_seconds, stf_post_phase_seconds,
                             cc_delay)
        stack_stream.append(tr_slice.copy())

    # Sufficient number of traces left?
    if len(stack_stream) < param.min_num_sta:
        return trace_list, len(stack_stream), None

    stack = np.sum([tr.data for tr in stack_stream], axis=0) / len(stack_stream)
    times = (stack_stream[0]).times()

    optimal_mccc_window['stack'] = stack
    optimal_mccc_window['time'] = times

    return trace_list, len(stack_stream), optimal_mccc_window


def find_dense_sta_patch(intra_station_dist, vnet, trace_list, search_radius_edge, log=sys.stdout):
    """Find patch of dense stations. Currently, we take the first station with this many neighbors to match Alex's
    code. We may have to revisit this later."""

    # Number of neighboring stations.
    neighbor_count = dict()

    # vnet parameters
    vnet_info = param.virtual_networks[vnet]

    # Stations that are within this distance are considered close neighbors.
    if search_radius_edge <= 0:
        return intra_station_dist, None

    # Check one station at a time.
    for center_key in trace_list.keys():
        center = trace_list[center_key]

        neighbor_count[center_key] = 0
        center_lat = center['lat']
        center_lon = center['lon']

        # Find how far are other stations from this station.
        for sta_net_key in trace_list.keys():
            if center_key == sta_net_key:
                continue
            neighbor = trace_list[sta_net_key]
            center_neighbor_key = f'{center_key}_{sta_net_key}'
            neighbor_center_key = f'{sta_net_key}_{center_key}'
            if center_neighbor_key in intra_station_dist:
                _dist, _distk, _azim, _back_azim = intra_station_dist[center_neighbor_key]
            else:
                _dist, _azim, _back_azim = gps2dist_azimuth(center_lat, center_lon,
                                                            neighbor['lat'],
                                                            neighbor['lon'])
                _distk = _dist / 1000.0
                _dist = kilometer2degrees(_distk)
                intra_station_dist[center_neighbor_key] = (_dist, _distk, _azim, _back_azim)
                intra_station_dist[neighbor_center_key] = (_dist, _distk, _azim, _back_azim)

            # If any station is within the search radius, it is considered a close neighbor.
            if _dist < search_radius_edge:
                neighbor_count[center_key] += 1
    if neighbor_count:
        sparse_patch_count = max(neighbor_count.values())
    else:
        sparse_patch_count = 0

    # Do we have a dense patch?
    dense_patch_list = dict()
    if sparse_patch_count <= vnet_info['sparse_patch_count']:
        print_message('INFO', f'No dense patches with more than {vnet_info["sparse_patch_count"]} neighbors found.',
                      log=log)
        return intra_station_dist, None
    else:
        # We take the first station with this many neighbors to match Alex's code. We may have to revisit this.
        s = ''
        if sparse_patch_count > 1:
            s = 's'
        print_message('INFO', f'The densest patch has {sparse_patch_count} station{s}.', log=log)
        for center_key in neighbor_count.keys():
            if sparse_patch_count == neighbor_count[center_key]:
                center_lat = trace_list[center_key]['lat']
                center_lon = trace_list[center_key]['lon']

                # Find how far are other stations from this station.
                for sta_net_key in trace_list.keys():
                    if center_key == sta_net_key:
                        dense_patch_list[sta_net_key] = False
                        continue
                    neighbor = trace_list[sta_net_key]
                    center_neighbor_key = f'{center_key}_{sta_net_key}'
                    neighbor_center_key = f'{sta_net_key}_{center_key}'

                    if f'{center_key}_{sta_net_key}' in intra_station_dist:
                        _dist, _distk, _azim, _back_azim = intra_station_dist[f'{center_key}_{sta_net_key}']
                    else:
                        _dist, _azim, _back_azim = gps2dist_azimuth(center_lat, center_lon,
                                                                    neighbor['lat'],
                                                                    neighbor['lon'])
                        _distk = _dist / 1000.0
                        _dist = kilometer2degrees(_distk)
                        intra_station_dist[center_neighbor_key] = (_dist, _distk, _azim, _back_azim)
                        intra_station_dist[neighbor_center_key] = (_dist, _distk, _azim, _back_azim)

                    # If any station is within the search radius, it is considered a close neighbor.
                    if _dist < search_radius_edge:
                        dense_patch_list[sta_net_key] = True
                    else:
                        dense_patch_list[sta_net_key] = False
                break
        return intra_station_dist, dense_patch_list


def p_wave_snr(trace, eq_datetime, p_time):
    """Computes SNR for the p-wave"""
    snr_window = param.snr_window
    noise_segment = trace.slice(starttime=eq_datetime + p_time + snr_window['pre-p'][0],
                                endtime=eq_datetime + p_time + snr_window['pre-p'][1])
    noise_max = np.max(np.abs(noise_segment.data))
    signal_segment = trace.slice(starttime=eq_datetime + p_time + snr_window['p'][0],
                                 endtime=eq_datetime + p_time + + snr_window['p'][1])

    signal_max = np.max(np.abs(signal_segment.data))
    snr = signal_max / noise_max

    return snr


def p_wave_normalize(trace, eq_datetime, p_time):
    """Normalize on a window around the phase."""
    signal_segment = trace.slice(starttime=eq_datetime + p_time - param.pre_phase_seconds,
                                 endtime=eq_datetime + p_time + param.post_phase_seconds)
    signal_max = np.max(np.abs(signal_segment.data))

    return signal_max


def xcorr_window(trace, eq_mag, eq_datetime):
    """Calculate length of the xcorr window starting pre_phase_seconds before the phase."""
    # 732.
    window_sum = list()
    for key in trace.keys():

        # Select a window around the phase.
        p_time = trace[key]['phase_delay']
        signal_segment = trace[key]['tr_final'].slice(starttime=eq_datetime + p_time - param.pre_phase_seconds,
                                                      endtime=eq_datetime + p_time + param.post_phase_seconds)
        signal_segment.normalize()

        # Get the trace envelope.
        data_envelope = envelope(signal_segment.data)

        # Select the window from start of the segment to just before the amplitude drops by 75%.
        # Note that the trace is already normalized in this window, so the max amplitude should be 1.
        times = (signal_segment.times()).copy()
        peak_flag = False
        for i, data in enumerate(data_envelope):
            if i > 0 and data >= 0.95:
                peak_flag = True

            # Want to check segment of the envelope that has a downward slope and is past the peak location.
            if peak_flag and data < 0.25:
                window_sum.append(times[i - 1] - times[0])
                break
    average_win = np.mean(window_sum)

    # We also have a window length based on the event's magnitude.
    mag_win = case(eq_mag, param.mag_window_length_dict)

    # The largest of the two window lengths is selected.
    xcorr_win_length = max(average_win, mag_win)

    return xcorr_win_length


def preprocess_trace(tr, filter_type, eq_datetime, phase_delay, normalize=True):
    """Preprocess the trace"""
    # Filter/normalize traces.
    tr_filter = tr.copy()
    tr_filter.detrend("linear")
    tr_filter.taper(max_percentage=0.05, type="hann")
    if filter_type is not None:
        freq = param.trace_filter[filter_type]['corners']
        tr_filter.filter('bandpass', freqmin=freq[0], freqmax=freq[1], corners=4, zerophase=True)
    tr_filter.detrend(type='demean')
    if normalize:
        tr_filter.data /= p_wave_normalize(tr_filter, eq_datetime, phase_delay)

    return tr_filter


def gen_synthetics(trace_list, eq_datetime, vn_name, create_synthetics=False):
    """Create synthetic seismograms."""
    # 1646
    synthetics = dict()
    if not create_synthetics:
        return synthetics

    for net_sta in trace_list.keys():
        trace = trace_list[net_sta]

        phase_delay = 0.0
        # Copy the structure of the filtered trace since this is the one we will be using for stacking.
        tr = trace[f'tr_filter'].copy()

        times = tr.times(reftime=eq_datetime)
        triangle_width = param.triangle_width
        triangle_half_width = triangle_width / 2.0
        syn_data = list()

        for t_index, t in enumerate(times):
            # Left half.

            if t < phase_delay - triangle_half_width:
                value = 0.0
            elif phase_delay - triangle_half_width <= t <= phase_delay:
                value = 1.0 + (t - phase_delay) * 1. / triangle_half_width

            # Right half
            elif phase_delay < t <= triangle_half_width + phase_delay:
                value = 1.0 - (t - phase_delay) * 1. / triangle_half_width
            else:
                value = 0.0

            syn_data.append(value)

        # Filter/normalize the trace.
        tr.data = np.array(syn_data.copy()) * trace['weight']
        tr_filter = preprocess_trace(tr, param.bp_filter[vn_name], eq_datetime, phase_delay, normalize=False)

        # For debugging only.
        debug = False
        if debug:
            trace_1 = tr
            trace_2 = tr_filter
            stream = Stream(traces=[trace_1, trace_2])
            stream.plot()
            sys.exit()

        synthetics[net_sta] = tr_filter.copy()

    return synthetics


def has_phase(trace, eq_datetime, phase_delay):
    """Check if seismogram covers a particular phase based on the phase_delay."""
    if (trace.stats.starttime <= eq_datetime + phase_delay - param.seconds_before_p) and \
            (trace.stats.endtime >= eq_datetime + phase_delay + param.seconds_after_p):
        return True
    return False


def set_time_parameters(eq_magnitude):
    """Set computation time parameters based on the earthquake magnitude
    """
    bp_t_offset = param.bp_time_offset

    bp_t_total = case(eq_magnitude, param.bp_time_total_dict)

    bp_t_increment = case(eq_magnitude, param.bp_time_inc_dict)

    # Total seconds of the Hanning averaging taper.
    t_avg = case(eq_magnitude, param.hann_t_avg_dict)

    stack_start = - bp_t_offset
    stack_end = stack_start + bp_t_total

    return stack_start, stack_end, bp_t_offset, bp_t_increment, bp_t_total, t_avg


def set_grid(eq_latitude, eq_longitude, eq_magnitude, grid_factor=1):
    """Set computation latitude grid increment."""

    decimal_places = param.grid_decimal_places

    # Set the latitude range and increment.
    latitude_inc = round(set_float_digits(case(eq_magnitude, param.grid_latitude_inc_dict)), decimal_places)
    latitude_half = round(set_float_digits(case(eq_magnitude, param.grid_latitude_half_dict)), decimal_places)
    # We assume latitude_half_count represents number of latitude_incs.
    latitude_half_count = int(round(latitude_half / latitude_inc))

    # Just make sure the half latitude is set as a multiple of latitude inc.
    latitude_half = round(latitude_inc * latitude_half_count, decimal_places)

    latitude_start = eq_latitude - latitude_half
    latitude_end = eq_latitude + latitude_half

    # Set the longitude range and increment base on the latitude.
    longitude_half = round(float(latitude_half / abs(math.cos(math.radians(eq_latitude)))), decimal_places)
    longitude_inc = round(float(longitude_half / latitude_half_count), decimal_places)
    # We assume longitude_half_count represents number of longitude_incs.
    longitude_half_count = int(round(longitude_half / longitude_inc))

    # Just make sure the half longitude is set as a multiple of longitude inc.
    longitude_half = round(longitude_inc * longitude_half_count, decimal_places)

    longitude_start = eq_longitude - longitude_half
    longitude_end = eq_longitude + longitude_half

    latitude = {'start': min(latitude_start, latitude_end), 'end': max(latitude_start, latitude_end),
                'inc': latitude_inc * grid_factor}
    longitude = {'start': min(longitude_start, longitude_end), 'end': max(longitude_start, longitude_end),
                 'inc': longitude_inc * grid_factor}

    return latitude, longitude


def phase_time(source_depth_km, source_distance):
    """Returns time of a phase in seconds."""

    arrivals = model.get_travel_times(source_depth_in_km=source_depth_km, distance_in_degree=source_distance,
                                      phase_list=param.phase)
    if arrivals:
        return arrivals[0].time
    else:
        return None


def get_bp_time_window(eq_date_time, eq_magnitude):
    """Set the request window based on the event magnitude and time."""

    # Base on "Get_date_for_FetchData_BP.f".
    # Start seconds before the event.
    t_before = param.request_time_before

    # Duration seconds after the event.
    t_after = case(eq_magnitude, param.request_time_after)

    event_datetime = UTCDateTime(eq_date_time)
    request_start_datetime = event_datetime - t_before
    request_start = request_start_datetime.strftime('%Y-%m-%d %H:%M:%S')
    request_start_date_time = request_start.replace(' ', 'T')
    request_end_datetime = event_datetime + t_after
    request_end = request_end_datetime.strftime('%Y-%m-%d %H:%M:%S')
    request_end_date_time = request_end.replace(' ', 'T')

    return request_start_date_time, request_start_datetime, request_end_date_time, request_end_datetime


def is_number(n):
    """Check if the input string input is a number.
    """
    try:
        float(n)
    except ValueError:
        return False
    return True


def is_loc_higher_priority(loc1, loc2):
    if loc1 in ['--', '00', '', '  ']:
        return False
    elif loc2 in ['--', '00', '', '  ']:
        return True
    elif is_number(loc1) and is_number(loc2):
        if int(loc2) < int(loc1):
            return True
        else:
            return False
    elif loc2 < loc1:
        return True
    return False


def get_2nd_point(lat1, lon1, radius, bearing=0):
    """Find latitude and longitude of a second point distance radius degrees from the first one."""

    # Radius of the Earth
    earth_radius = param.earth_radius

    # Bearing is converted to radians.
    bearing = math.radians(bearing)

    # Distance in km.
    d = degrees2kilometers(radius)

    # Current lat,lon point converted to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)

    lat2 = math.asin(math.sin(lat1) * math.cos(d / earth_radius) +
                     math.cos(lat1) * math.sin(d / earth_radius) * math.cos(bearing))

    lon2 = lon1 + math.atan2(math.sin(bearing) * math.sin(d / earth_radius) * math.cos(lat1),
                             math.cos(d / earth_radius) - math.sin(lat1) * math.sin(lat2))

    # The new point.
    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)

    return lat2, lon2


def get_service_url(ws_catalog, ws_dc):
    """Extract the service URL from dataselect service URL."""
    ws_service_url = ws_catalog[ws_dc].dataselect.split('/fdsn')[0]
    return ws_service_url


def get_request_items(req_lines):
    """Split a request line to its components."""
    net_sta_rec = dict()
    bulk_rec = list()
    for _line in req_lines:
        _net, _sta, _loc, _chan, _start, _end = _line.strip().split()
        net_sta_key = f'{_net}-{_sta}'
        if net_sta_key not in net_sta_rec.keys():
            net_sta_rec[net_sta_key] = [_net, _sta, _loc, _chan, _start, _end]
        else:
            if channel_order[_chan] < channel_order[net_sta_rec[net_sta_key][3]]:
                net_sta_rec[net_sta_key] = [_net, _sta, _loc, _chan, _start, _end]
            elif channel_order[_chan] == channel_order[net_sta_rec[net_sta_key][3]]:
                if is_loc_higher_priority(net_sta_rec[net_sta_key][2], _loc):
                    net_sta_rec[net_sta_key] = [_net, _sta, _loc, _chan, _start, _end]
    for _key in net_sta_rec.keys():
        bulk_rec.append(net_sta_rec[_key])

    return bulk_rec


def is_net_temporary(net):
    """Exclude temporary networks."""
    if len(net) <= 2:
        if net[0].isdigit():
            return True
        if net[0].lower() in ['x', 'y', 'z']:
            return True
    return False


def read_url(target_url, log=sys.stdout, verbose=False):
    """Read content of a URL."""
    if verbose:
        print_message('INFO', f'Opening URL: {target_url}', log=log)

    with urlopen(target_url) as url:
        content = url.read().decode()
    return content


def get_fedcatalog_stations(start, end, lat, lon, rmin, rmax, net='*', dc='*', req='request', service='dataselect', log=sys.stdout):
    url = f'{fedcatalog_service_url}net={net}&cha={param.request_channel}&starttime={start}&endtime={end}' \
          f'&targetservice={service}&level=channel&datacenter={dc}&format={req}&includeoverlaps=false' \
          f'&endafter={start}&lat={lat}&lon={lon}&minradius={rmin}&maxradius={rmax}&includerestricted=false&&nodata=404'

    print_message('INFO', f'Requesting: {url}', log=log)
    try:
        content = read_url(url, log=log)
    except Exception as _er:
        print_message('ERR', f'Request  {url}: {_er}', log=log)
        return None

    return content


def get_dc_dataselect_url(start, end, net, sta, loc, chan, dc='*', req='request', service='dataselect', log=sys.stdout):
    url = f'{fedcatalog_service_url}starttime={start}&endtime={end}' \
          f'&targetservice={service}&level=channel&datacenter={dc}&format={req}&includeoverlaps=false' \
          f'&net={net}&sta={sta}&cha={chan}&loc={loc}&nodata=404'

    print_message('INFO', f'Requesting: {url}', log=log)
    try:
        content = read_url(url, log=log)
    except Exception as _er:
        print_message('ERR', f'Request  {url}: {_er}', log=log)
        return None
    for _line in content.split('\n'):
        if 'DATASELECTSERVICE' in _line:
            _url = _line.replace('DATASELECTSERVICE=', '')
            return _url.strip()


def split_fedcatalog_stations(station_data, log=sys.stdout):
    """Get station list from fedcatalog service."""

    # This dictionary provides a template for fetdatalog creation.
    catalog_info = dict()

    bulk_list = dict()

    # Go through the station lines and split them to data centers.
    _lines = station_data.split('\n')

    _line_index = -1
    previous_dc = None
    dc_name = None
    for _line in _lines:
        _line_index += 1

        # Skip the blank and the comment lines.
        if not _line.strip():
            continue
        if _line.startswith('#') and '=' not in _line:
            continue

        # From the parameter=value lines, we are interested in the DATACENTER and DATASELECTSERVICE lines.
        elif _line.startswith('#') and '=' in _line:
            _par, _value = _line.split('=')

            # Found the data center name.
            if _par == '#DATACENTER':
                if dc_name is not None:
                    previous_dc = dc_name
                print_message('INFO', f'from the {_value} data center', log=log)
                dc_name, dc_url = _value.strip().split(',')

                # Initialize the data center information.
                if dc_name not in catalog_info.keys():
                    print_message('INFO', f'Initiating fedcatalog request for {dc_name}', log=log)
                    catalog_info[dc_name] = ObjDict({'url': dc_url, 'dataselect': '', 'bulk': []})

                # if this is not the first data center, save the previous data center's bulk list
                if bulk_list:
                    catalog_info[previous_dc].bulk = list()
                    for _key in bulk_list:
                        catalog_info[previous_dc].bulk.append(bulk_list[_key]['line'])

                    # The list is saved. Now, reset the bulk_list.
                    bulk_list = dict()

                continue
            # Found the dataselect service address.
            elif _par == 'DATASELECTSERVICE':
                # Save the dataselect service address in the catalog for this DC.
                catalog_info[dc_name].dataselect = _value.strip()
                print_message('INFO', f'dataselect service is {_value.strip()}', log=log)
                continue
            elif _par == 'STATIONSERVICE':
                # Save the dataselect service address in the catalog for this DC.
                catalog_info[dc_name].dataselect = _value.strip()
                print_message('INFO', f'station service is {_value.strip()}', log=log)
                continue
            else:
                # Ignore the other definitions.
                continue

        # The rest are the station lines.
        else:
            # Skip the blank lines.
            _line = _line.strip()
            if not _line:
                continue

            # Insert channels based on the channel priority list.
            try:
                items = _line.split('|')
                _net, _sta, _loc, _chan = items[0:4]
            except Exception as ex:
                print_message('ERR', f'Failed to split station line {_line} to (_net, _sta, _loc, _chan)\n{ex}',
                              log=log)
                continue
            _key = '_'.join([_net, _sta])
            if _key not in bulk_list:
                bulk_list[_key] = {'chan': _chan, 'line': _line}
            else:
                # If the station exists, just save the lower ordered channel.
                if param.channel_order[_chan] < param.channel_order[bulk_list[_key]['chan']]:
                    bulk_list[_key] = {'chan': _chan, 'line': _line}

    # Save the last data center's bulk list.
    if bulk_list:
        catalog_info[dc_name].bulk = list()
        for _key in bulk_list:
            for _chan in bulk_list[_key]['chan']:
                catalog_info[dc_name].bulk.append(bulk_list[_key]['line'])

    return ObjDict(catalog_info)


def trace_amp(tr_times, tr_data, t1, sampling=None):
    """Get trace amplitude at a give time."""
    if sampling is None:
        t_index = math.floor((t1 - tr_times[0]) * sampling)

    return tr_data[t_index]


def hann(stack, syn_stack, start, end, bp_t_increment, average_seconds, resamp=1, create_synthetics=False,
         log=sys.stdout):
    """Integrate the beams in a running beam_average_seconds  window"""

    print_message('INFO', f'In Hann: Integrate {len(stack)} beams over a {average_seconds} s window', log=log)
    # The Hanning function of length num_points is used to perform Hanning smoothing
    num_points = int(round(average_seconds / bp_t_increment)) + 1

    # Half points. We always want odd numbers, so Hanning amplitude of 1 falls on the current sample.
    if num_points % 2 == 0:
        num_half_points = int(num_points / 2)
        num_points = num_points + 1
    else:
        num_half_points = int((num_points - 1) / 2)

    # The Hanning window.
    _hann = np.hanning(num_points)

    # Move the start and end by half of the averaging window length.
    half_window = num_half_points * bp_t_increment
    stack_start = start - half_window
    stack_end = end + half_window
    print_message('INFO', f'In Hann: half_window {half_window}, stack_start {stack_start}, '
                          f'stack_end {stack_end}', log=log)

    smooth_stack = dict()
    syn_smooth_stack = dict()

    # Loop over grid points.
    for grid_key in stack:
        smooth_stack[grid_key] = dict()
        if create_synthetics:
            syn_smooth_stack[grid_key] = dict()

        # Apply smoothing. First get a list of sample times at each grid point.
        time_key_list = list(stack[grid_key])
        step_count = 0

        # Step through the time samples.
        for time_index, time_key in enumerate(time_key_list):

            # Skip samples (resample) if necessary. But always start from the first sample.
            step_count += 1
            if time_index == 0:
                step_count = 0
            elif step_count != resamp:
                continue
            else:
                step_count = 0

            # Initialize.
            smooth_stack[grid_key][time_key] = 0.0
            if create_synthetics:
                syn_smooth_stack[grid_key][time_key] = 0.0

            # Loop over all samples within the trace window.
            if stack_start <= float(time_key) <= stack_end:
                # Then Hanning index.
                hann_index = 0
                # Index of data point with respect to the current point index.
                this_point = - num_half_points
                # Counter.
                count = 0

                # Going from -num_half_points to +num_half_points.
                while this_point <= num_half_points:
                    debug = False
                    if debug:
                        print("this_point:", this_point)
                        print("time_index:", time_index)
                        print("time_key:", time_key)
                        print("time_key_list:", time_key_list)

                    # Make sure not to go beyond the list length.
                    if 0 <= time_index + this_point <= len(time_key_list) - 1:
                        t_key = time_key_list[time_index + this_point]
                        # cos(pi*x/L)**2    abs(x) <= L/2
                        try:
                            smooth_stack[grid_key][time_key] = smooth_stack[grid_key][time_key] + \
                                                               stack[grid_key][t_key] * _hann[hann_index]
                        except Exception as ex:
                            print_message('ERR', f'KeyError: {t_key}, {time_key}', log=log)

                        if create_synthetics:
                            try:
                                syn_smooth_stack[grid_key][time_key] = syn_smooth_stack[grid_key][time_key] + \
                                                                       syn_stack[grid_key][t_key] * _hann[hann_index]
                            except Exception as ex:
                                print_message('ERR', f'Synthetics KeyError: {t_key}, {time_key}', log=log)
                    this_point += 1
                    count += 1
                    hann_index += 1

                smooth_stack[grid_key][time_key] /= count
                if create_synthetics:
                    syn_smooth_stack[grid_key][time_key] /= count

    return smooth_stack, syn_smooth_stack


def smooth(x, window_length, log=sys.stdout):
    """smooth the data using a Hanning window of the requested size."""
    if x.size < window_length:
        print(f'Input vector length {x.size} is smaller than the '
              f'requested window length of {window_length}.')
        raise ValueError

    print_message('INFO', f'In Smooth: window_length {window_length}', log=log)

    # np.r_ stacks the comma-separated  arraysalong their first axis.
    signal = np.r_[x[window_length - 1:0:-1], x, x[-2:-window_length - 1:-1]]
    weights = np.hanning(window_length)
    smooth_signal = np.convolve(weights / weights.sum(), signal, mode='valid')
    return smooth_signal


def stack_root(stack, syn_stack, stacking_root, create_synthetics=False):
    """If necessary, raise stacks to the Nth power to do Nth-root stacking."""
    if stacking_root != 1.0:
        for grid_key in stack:
            for t_key in stack[grid_key]:
                stack[grid_key][t_key] = np.power(stack[grid_key][t_key], stacking_root)
                if create_synthetics:
                    syn_stack[grid_key][t_key] = np.power(syn_stack[grid_key][t_key],
                                                          stacking_root)

    return stack, syn_stack


def trace_trimmer(trace, origin_time, bp_t_offset, bp_t_total, phase_delay, mccc_delay):
    """Trim a trace and pad with zeros and include delays"""
    trimmed_trace = trace.copy()

    # Trim the trace based on the origin time and delays and BP offset.
    _start = origin_time + phase_delay + mccc_delay - bp_t_offset
    _end = _start + bp_t_total
    trimmed_trace.trim(starttime=_start, endtime=_end, pad=True, nearest_sample=False, fill_value=0.0)

    # Reset the start time without including the delays so the corrections will always be included.
    trimmed_trace.stats.starttime = origin_time - bp_t_offset
    return trimmed_trace


def shift_trace(trace, origin_time, bp_t_offset, bp_t_total, shift_time):
    """Shift a trace by the time given shift_time"""
    shifted_trace = trace.copy()

    # Trim the trace based on the origin time and delays.
    _start = trace.stats.starttime
    _start += shift_time
    _end = _start + bp_t_total
    shifted_trace.trim(starttime=_start, endtime=_end, pad=True, nearest_sample=False, fill_value=0.0)

    # Reset the start time without including the delays so the corrections will always be included.
    shifted_trace.stats.starttime = origin_time - bp_t_offset
    return shifted_trace


def get_phase_delay(tt_cache, dist, depth):
    """Obtauin the phase delay for a give distance either directly or from a cache."""
    # Testing shows that 2 decimal places (0.01 deg resolution) improves the match and speeds up the code while it has
    # little affect on the quality.
    dist /= 1000.0
    dist = kilometer2degrees(dist)
    dist_key = f'{dist:0.2f}'
    if dist_key in tt_cache:
        phase_delay = tt_cache[dist_key]
    else:
        # Delay from predicted travel time.
        phase_delay = phase_time(depth, dist)
        tt_cache[dist_key] = phase_delay
    return tt_cache, phase_delay


def stacker(trace_list, vn_name, bp_t_offset, t_avg, eq_datetime, eq_lat, eq_lon, eq_mag, eq_depth, bp_t_increment,
            bp_t_total, tt_cache, create_synthetics=False, grid_factor=1, log=sys.stdout, verbose=False):
    """Stack traces using the BP window and step."""
    # inum=nint((ittotal+tavg)/fitincrement)+1
    # ittotal = bp_t_total
    # fitincrement = bp_time_inc_dict
    # tavg = t_avg  Total seconds of the Hanning averaging taper
    # i1len = len(lon)
    # i2len = len(lat)
    # ntotal = total number of stations
    # npts = number of samples
    # dt() trace sampling interval
    # 1705
    # BP start one increment back.
    # bp_initial_seconds, seconds before the event time when BP starts.

    # Delta time (seconds) before the event time when BP starts.
    seconds_before = math.ceil(param.bp_initial_seconds - 1.0 * bp_t_offset - t_avg / 2.0 - bp_t_increment)

    # Delta time (seconds) after the event time when BP ends.
    seconds_after = math.ceil(seconds_before + bp_t_total + t_avg)

    print_message('INFO', f'start {seconds_before}s and end {seconds_after}s relative to the '
                          f'event time of {eq_datetime}, using bp_t_offset of '
                          f'{bp_t_offset}s,  '
                          f't_avg {t_avg}s, and bp_t_total {bp_t_total}s',
                  log=log)
    stacking_root = param.stacking_root

    # Set the grid points around the earthquake location.
    latitude, longitude = set_grid(eq_lat, eq_lon, eq_mag, grid_factor=grid_factor)
    print_message('INFO', f'Grid latitude {latitude}, longitude {longitude}', log=log)
    if verbose:
        print_message('INFO', f'Earthquake at {eq_lat}, {eq_lon}')
        lat_ = latitude['start'] - latitude['inc']
        lat_str = ''
        while lat_ < latitude['end']:
            lat_ += latitude['inc']
            lat_str = f'{lat_str}, {lat_:0.3f}'

        lon_ = longitude['start'] - longitude['inc']
        lon_str = ''
        while lon_ < longitude['end']:
            lon_ += longitude['inc']
            lon_str = f'{lon_str}, {lon_:0.3f}'
        print_message('INFO', f'Grid latitudes: {lat_str}')
        print_message('INFO', f'Grid longitudes: {lon_str}')

    global_max = None
    t0 = time.time()

    warn_list = list()
    stack_list = list()

    stack = dict()
    syn_stack = dict()

    # Individual stations.
    sta_counter = 0

    # Get ready to loop through the grid.
    lat_start = latitude['start']
    lat_inc = latitude['inc']
    lat_end = latitude['end']

    lon_start = longitude['start']
    lon_inc = longitude['inc']
    lon_end = longitude['end']

    net_sta_list = list()

    # Loop through station list and extract information on active traces. Here we want to avoid
    # repeating this over every grid point.
    for net_sta in trace_list:

        # Check the trace to process to see if it is active. (#1727)
        trace = trace_list[net_sta]

        # The station location.
        lat0, lon0 = trace['lat'], trace['lon']

        # Keep track of stations in the stack.
        if net_sta not in stack_list:
            stack_list.append(net_sta)
            sta_counter += 1
            print_message('INFO', f'Adding {net_sta} to stack list ({sta_counter}), weight = {trace["weight"]:0.2f}',
                          log=log)

        # Trace information.
        tr = trace[f'tr_filter']
        tr_starttime = tr.stats.starttime
        tr_endtime = tr.stats.endtime
        net_sta_list.append((net_sta, tr.copy(), float(lat0), float(lon0), tr_starttime,
                             tr_endtime, trace['phase_delay']))

        # Synthetic trace information.
        if create_synthetics:
            syn_tr = trace['tr_syn'].copy()

    grid_key_list = list()

    # Loop over the grid latitudes.
    grid_index = 0
    grid_lat = lat_start - lat_inc
    while grid_lat < lat_end:
        grid_lat += lat_inc

        grid_lon = lon_start - lon_inc
        # 1721
        while grid_lon < lon_end:
            grid_index += 1
            grid_lon += lon_inc

            # Grid points.
            grid_key = f'{grid_lat:0.3f}_{grid_lon:0.3f}'
            grid_key_list.append((grid_index, grid_key, grid_lat, grid_lon))

    n_grids = len(grid_key_list)
    grid_count = 0
    show_message = True
    for (grid_index, grid_key, grid_lat, grid_lon) in grid_key_list:
        grid_count += 1

        # Loop through grid longitudes.
        shifted_trace = None
        if grid_key not in stack:
            _stack = None
        stack[grid_key] = dict()

        if create_synthetics:
            if grid_key not in syn_stack:
                _syn_stack = None
            syn_stack[grid_key] = dict()

        # For each new source location (grid_lat, grid_lon), loop through stations and stack.
        for (net_sta, tr, lat0, lon0, tr_starttime, tr_endtime, phase_delay0) in net_sta_list:

            # Assume event is at the grid point. Get the distance from the station.
            _dist, _azim, _back_azim = gps2dist_azimuth(grid_lat, grid_lon, lat0, lon0)

            tt_cache, phase_delay = get_phase_delay(tt_cache, _dist, eq_depth)

            if phase_delay is None:
                if net_sta not in warn_list:
                    warn_list.append(net_sta)
                    print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta} '
                                          f'because no P-wave arrival', log=log)
                continue

            # The time_shift variable represents the number of seconds difference in P- travel time due to
            # event relocation.
            time_shift = phase_delay - phase_delay0

            # Skip if trace window is not covered by BP.
            if eq_datetime + seconds_before >= tr_endtime + time_shift or \
                    eq_datetime + seconds_after <= tr_starttime + time_shift:
                print_message('WARN', f'Skipped, channel {tr.stats.channel} of Net.Sta {net_sta} '
                                      f'because it does not fall in BP window', log=log)
                continue

            # Here we shift the trace in time based on the time_shift.
            shifted_trace = tr.copy()

            # time_shift < 0 trace must be shifted to the left (arrives sooner, trim window slides right).
            # time_shift > 0 trace must be shifted to the right (arrives later, trim window slides left).
            # So we use a negative sign to shift the trim window accordingly.
            _start = tr_starttime + time_shift
            _end = tr_endtime + time_shift
            # Slide the trace based on the delta time shift due to new event location.
            shifted_trace.trim(starttime=_start, endtime=_end, pad=True, nearest_sample=False, fill_value=0.0)
            if verbose and show_message:
                print_message('INFO', f'Slide the trace based on the delta time shift due to new event '
                                      f'location start:{_start} and end:{_end}', log=log)

            # Now, reset  the start time to that of the original trace so the time shift will be implicit.
            shifted_trace.stats.starttime = tr_starttime

            # Set the trace window to the BP time.
            shifted_trace.trim(starttime=eq_datetime + seconds_before, endtime=eq_datetime + seconds_after,
                               pad=True, nearest_sample=False, fill_value=0.0)
            if verbose and show_message:
                show_message = False
                print_message('INFO', f'Set the trace window to the BP time '
                                      f'start:eq_datetime -{abs(seconds_before)}s '
                                      f'and end:eq_datetime + {seconds_after}s',
                              log=log)
            # Resample to bp increment.
            shifted_trace.resample(1.0 / bp_t_increment)

            _data = shifted_trace.data
            _data_sign = np.sign(_data)
            _data_n = np.power(np.abs(_data), 1.0 / stacking_root)
            _data_n *= _data_sign

            # Perform stacking.
            if _stack is None:
                _stack = _data_n.copy()
            else:
                _stack = np.add(_stack, _data_n)

            # For debugging only.
            debug = False
            if debug:
                import matplotlib.pyplot as plt
                fig = plt.figure(figsize=(9.5, 11))
                ax = fig.add_subplot(1, 1, 1)
                ax.plot(shifted_trace.times(reftime=eq_datetime),
                        _stack, "k-", lw=0.4, label=f'{grid_lat}, {grid_lon}/ {eq_lat}, {eq_lon}')
                plt.legend()
                plt.show()
                plt.close()

            # Do the same as above with synthetics.
            if create_synthetics:
                shifted_syn = syn_tr.copy()
                shifted_syn.trim(starttime=_start, endtime=_end, pad=True, nearest_sample=False, fill_value=0.0)
                shifted_syn.stats.starttime = tr_starttime

                # Set the trace window to the BP time.
                shifted_syn.trim(starttime=eq_datetime + seconds_before, endtime=eq_datetime + seconds_after,
                                 pad=True, nearest_sample=False, fill_value=0.0)

                shifted_syn.resample(1.0 / bp_t_increment)

                _data = shifted_syn.data
                _data_sign = np.sign(_data)
                _data_n = np.power(np.abs(_data), 1.0 / stacking_root)
                _data_n *= _data_sign

                if _syn_stack is None:
                    _syn_stack = _data_n.copy()
                else:
                    _syn_stack = np.add(_syn_stack, _data_n)

                # For debugging only.
                debug = False
                if debug:
                    import matplotlib.pyplot as plt
                    fig = plt.figure(figsize=(9.5, 11))
                    ax = fig.add_subplot(1, 1, 1)
                    ax.plot(shifted_syn.times(reftime=eq_datetime),
                            _syn_stack, "k-", lw=0.4, label=f'{grid_lat}, {grid_lon}/ {eq_lat}, {eq_lon}')
                    plt.legend()
                    plt.show()
                    plt.close()

        # Store stack for this grid point.
        # 1726
        if shifted_trace is not None:
            stack[grid_key] = {f'{t :0.2f}': _stack[t_index] for t_index, t in
                               enumerate(shifted_trace.times(reftime=eq_datetime))}
        else:
            print_message('WARN', f'stack trace empty!', log=log)

        if create_synthetics:
            if shifted_syn is not None:
                syn_stack[grid_key] = {f'{t :0.2f}': _syn_stack[t_index] for t_index, t in
                                       enumerate(shifted_syn.times(reftime=eq_datetime))}
            else:
                print_message('WARN', f'synthetic stack trace empty!', log=log)

        if grid_count == 100:
            t0 = time_it(t0, stage=f'All stations for grid points {grid_index}/{n_grids}', log=log)
            grid_count = 0

    return stack, syn_stack

