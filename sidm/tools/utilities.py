"""Module to define miscellaneous helper methods"""

import os
import yaml
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import hist.intervals


def print_list(l):
    """Print one list element per line"""
    print('\n'.join(l))

def print_debug(name, val, print_mode=True):
    """Print variable name and value"""
    if print_mode:
        print(f"{name}: {val}")

def partition_list(l, condition):
    """Given a single list, return separate lists of elements that pass or fail a condition"""
    passes = []
    fails = []
    for x in l:
        if condition(x):
            passes.append(x)
        else:
            fails.append(x)
    return passes, fails

def flatten(x):
    """Flatten arbitrarily nested list or dict"""
    # https://stackoverflow.com/questions/2158395/
    flattened_list = []
    def loop(sublist):
        if isinstance(sublist, dict):
            sublist = sublist.values()
        for item in sublist:
            if isinstance(item, (dict, list)):
                loop(item)
            else:
                flattened_list.append(item)
    loop(x)
    return flattened_list

def add_unique_and_flatten(flattened_list, x):
    """Flatten arbitrarily nested list or dict, keeping only unique items"""
    # https://stackoverflow.com/questions/2158395/
    def loop(sublist):
        if isinstance(sublist, dict):
            sublist = sublist.values()
        for item in sublist:
            if isinstance(item, (dict, list)):
                loop(item)
            elif item not in flattened_list:
                flattened_list.append(item)
    loop(x)
    return flattened_list

def as_int(array):
    """Return array with values converted to ints"""
    return ak.values_astype(array, "int64")

def dR(obj1, obj2):
    """Return dR between obj1 and the nearest obj2; returns None if no obj2 is found"""
    return obj1.nearest(obj2, return_metric=True)[1]

def dR_outer(obj1, obj2):
    """Return dR between outer tracks of obj1 and obj2"""
    return np.sqrt((obj1.outerEta - obj2.outerEta)**2 + (obj1.outerPhi - obj2.outerPhi)**2)

def drop_none(obj):
    """Remove None entries from an array (not available in Awkward 1)"""
    return obj[~ak.is_none(obj, axis=1)] # fixme: not clear why axis=1 works and axis=-1 doesn't

def matched(obj1, obj2, r):
    """Return set of obj1 that have >=1 obj2 within r; remove None entries before returning"""
    return drop_none(obj1[dR(obj1, obj2) < r])

#def lxy(obj):
#    """Return transverse distance between production and decay vertices"""
#    return (obj - ak.firsts(obj.children)).r

def lxy(obj):
    """Return transverse distance between production and decay vertices"""
    return (obj - ak.firsts(obj.children, axis=2)).r

def rho(obj, ref=None):
    """Return transverse distance between object and reference (default reference is 0,0)"""
    ref_x = ref.x if ref is not None else 0.0
    ref_y = ref.y if ref is not None else 0.0
    return np.sqrt((obj.x - ref_x)**2 + (obj.y - ref_y)**2)

def set_plot_style(style='cms', dpi=50):
    """Set plotting style using mplhep"""
    if style == 'cms':
        plt.style.use(hep.style.CMS)
    else:
        raise NotImplementedError
    plt.rcParams['figure.dpi'] = dpi

def plot(hists, skip_label=False, **kwargs):
    """Plot using hep.hist(2d)plot and add cms labels"""
    dim = len(hists[0].axes) if isinstance(hists, list) else len(hists.axes)
    if dim == 1:
        h = hep.histplot(hists, **kwargs)
    elif dim == 2:
        h = hep.hist2dplot(hists, **kwargs)
    else:
        raise NotImplementedError(f"Cannot plot {dim}-dimensional hist")
    if not skip_label:
        hep.cms.label()
    return h

def get_eff_hist(num_hist, denom_hist):
    """Returns the histogram of num_hist/denom_hist and a 2D numpy array of the up/down errors on the efficiency. Plot the errors using yerr=errors when plotting. """
    denom_vals  = denom_hist.values()
    num_vals   = num_hist.values()

    errors = hist.intervals.ratio_uncertainty(num_vals,denom_vals,'efficiency')
    small_num= 1e-9
    eff_values = num_vals/(denom_vals + small_num)
   

    eff_hist = hist.Hist(*num_hist.axes)
    eff_hist.values()[:] = eff_values
    return eff_hist, errors

def load_yaml(cfg):
    """Load yaml files and return corresponding dict"""
    cwd = os.path.dirname(os.path.abspath(__file__))
    with open(f"{cwd}/{cfg}", encoding="utf8") as yaml_cfg:
        return yaml.safe_load(yaml_cfg)

def make_fileset(samples, ntuple_version, max_files=-1, location_cfg="../configs/ntuple_locations.yaml"):
    """Make fileset to pass to processor.runner"""
    ntuple_versions = ["ffntuple_v2", "ffntuple_v4", "llpNanoAOD_v1", "llpNanoAOD_v2", "ffntuple_official", "ffntuple_private", "llpNanoAOD_v2_merged"]
    if ntuple_version not in ntuple_versions:
        raise NotImplementedError(f"Only {ntuple_versions} ntuples have been implemented")
    locations = load_yaml(location_cfg)[ntuple_version]
    fileset = {}
    for sample in samples:
        base_path = locations["path"] + locations["samples"][sample]["path"]
        file_list = [base_path + f for f in locations["samples"][sample]["files"]]
        if max_files != -1:
            file_list = file_list[:max_files]
        fileset[sample] = file_list
    return fileset

def check_bit(array, bit_num):
    """Return boolean stored in the bit_numth bit of array"""
    return (array & pow(2, bit_num)) > 0

def get_hist_mean(hist):
    """Return mean of 1D histogram"""
    return np.atleast_1d(hist.profile(axis=0).view())[0].value
