# A couple of Python classes for reading in MST data from the MDSplus
# servers and doing some simple analysis tasks.
#
# Written by Mark Nornberg 3/6/2012 
# (with some code copied from J.D. Lee's mdsplus module)

# We load all of the relevant modules from standard libraries...
from os import path, mkdir, remove
import shelve
import datetime

# and from third-party libraries...
from numpy import *
#from scipy.io import idl
from scipy.signal import resample, cspline1d, cspline1d_eval, medfilt
import pylab as p
from MDSplus import Connection

from analysis_util import *

# Here is where you should pick the location of your local data
# storage directory
store_data_locally = False
local_path = path.expanduser('~/p/mst/data')
if not path.isdir(local_path): mkdir(local_path)

# Here's where I shamelessly steal JD Lee's code for Global caches. 
# It turns out that opening and closing connections to the server can 
# cause the server to become unresponsive. This may be because the mdsip 
# processes on the server get locked up some how, I'm not sure. 
_SVR_CONNS = {} # Connection objects indexed by server address.
_SVR_SHOTS = {} # Shot numbers indexed by server address.
_SVR_TREES = {} # Tree names indexed by server address.

def _get_svr_cached(svr):
    """For the given server address, return the cached connection object,
    tree name, and loaded shot number.
    """
    try:
        return _SVR_CONNS[svr], _SVR_TREES[svr], _SVR_SHOTS[svr] 
    except:
        return None, None, None
    
    
def _update_svr_cache(svr, conn, tree, shot):
    _SVR_CONNS[svr] = conn
    _SVR_TREES[svr] = tree
    _SVR_SHOTS[svr] = shot


def get_server_for_shot(shot):
    """Get the server address for the given shot number. Shots for the 
    current day are in aurora, while previous days are on dave.
    """
    first = min_shot_for_date(datetime.date.today())
    if shot >= first:
        return 'aurora.physics.wisc.edu'
    else:
        return 'dave.physics.wisc.edu'

def min_shot_for_date(date):
    """Return the minimum shot number for the given date."""
    y = date.year
    m = date.month
    d = date.day
    return (y % 1000 + 100) * 10000000 + m * 100000 + d * 1000 + 1


def max_shot_for_date(date):
    """Return the maximum shot for the given date."""
    return min_shot_for_date(date) + 998


def shot_to_date(shot):
    """Return the date corresponding to the given shot number."""
    day = int(shot / 1000) % 100
    month = int(shot / 100000) % 100
    year = int(shot / 10000000) + 1900
    return datetime.date(year, month, day)


def shot_to_date_num(shot):
    """Given a shot, return a number corresponding to the full date. 
    For example, given 1100502001, return 20100502. 
    """
    return date_to_date_num(shot_to_date(shot))


def get_connection(shot, tree='mst'):

    """Get an MDSplus connection object connected to the appropriate 
    server and having the given tree and shot opened. 
    """
    svr = get_server_for_shot(shot)
    conn, tree_, shot_ = _get_svr_cached(svr)
        
    if conn is None:
        conn = Connection(svr)
        
    if tree != tree_ or shot != shot_:
        try:
            # This throws an exception if there are no open trees. 
            conn.closeAllTrees()
        except:
            pass
        conn.openTree(tree, shot)
        
    _update_svr_cache(svr, conn, tree, shot)
    
    return conn

        

class Signal(dict):

    """A signal class to provide an organiational scheme for storing
    timeseries data."""

    def __init__(self, shot, nodepath, dims = [],
                 tree='mst'):

        # Use either the defaults or specified values for the tree,
        # server, shot and nodepath to the data on the MDSPlus server
        self['shot'] = shot 
        self['tree'] = tree
        self['nodepath'] = nodepath
        self['server'] = get_server_for_shot(shot)

        # The name of this signal will be whatever is at the end of
        # the nodepath
        self['name'] = (nodepath.split(':'))[-1]

        # Retreive the data from the MDSPlus tree on dave.physics.wisc.edu

        try:
            conn = get_connection(shot, tree)
            sig_data = conn.get(nodepath)
            self['signal'] = sig_data.data()

        # If there is no data available in the node then there is no
        # point in doing anything more. Close the MDSPlus connection
        # and return.
        except: 
            print self['nodepath'] + " is not availble for shot " + \
            str(self['shot']) + " on server " + self['server'] + '.'
            return

        # Is the data we obtained a string? If so, we're done
        if type(self['signal']) == type(''): return

        # What if the data is more than 1D? Let's check and hope for
        # clarity in the dims parameter provided.
        ndims = self['signal'].shape

        if len(ndims) == 1: 
            try: 
                time_data = conn.get('dim_of({0})'.format(nodepath))
                self['time'] = time_data.data()
            except: pass
        elif len(ndims) == len(dims):
            for d in dims:
                self[d] = conn.get('dim_of({0},{1})'.format(nodepath, 
                                                            dims.index(d)))
        else:
            print "Number of dimensions of " + nodepath 
            + " is not the same as \'dims\' argument."

        # Let's see if there's any more meta data like units and error bars
        try: 
            self['units'] = conn.get('units_of({0})'.format(nodepath))
            # If the units are just whitespace then delete the object
            if self['units'].split() == []: del self['units']
        except: pass

        try:
            t_units_str = 'units_of(dim_of({0}))'.format(nodepath)
            self['time_units'] = conn.get(t_units_str)
            # If the units are just whitespace then delete the object
            if self['time_units'].split() == []: del self['time_units']
        except: pass

        try: 
            err_str = 'error_of({0})'.format(nodepath)
            self['error'] = conn.get(err_str)
            self['min'] = self['signal'] - self['error']
            self['max'] = self['signal'] + self['error']
        except: pass


    def plot(self, *args, **keywords):
        
        """Plot the signal. Arguments and keywords are passed to the
        matplotlib plot method.

        This method clears the current figure by default. 
        Set hold=True to overplot."""

        # Create a new signal plot by default
        if not 'hold' in keywords: keywords['hold'] = False

        # Define the label based on the signal name including units if
        # available
        if not 'label' in keywords: keywords['label'] = self['name']

        # Plot the signal
        p.plot(self['time'], self['signal'], *args, **keywords)
        
        if ('max' in self) and ('min' in self):
            band_keywords = {}
            if 'color' in keywords: 
                band_keywords['color'] = keywords['color']
            plot_band(self['time'], self['signal'], 
                      ymin=self['min'], ymax=self['max'], **band_keywords)

        # Add what annotation we can
        try: p.xlabel(self['time_units'])
        except: pass
        try: p.ylabel(self['units'])
        except: pass

        p.legend()


    def make_empty(name):
        
        """A method to create the signal object without the MDSPlus
        calls to fill the entries."""

        # Let's set some defaults for the members of our object
        self['tree'] = ''
        self['server'] = ''
        self['shot'] = None
        self['name'] = name
        self['signal'] = None
        self['time'] = None


    def resample(self, num, axis=0, window=None):

        """Resample timeseries data for a signal onto a different time
        base. This routine uses scipy.signal's resample function which
        creates an interpolation with the option of filtering by
        specifying the window keyword. See help(scipy.signal.resample)
        for details."""

        if self['time'] is not None:
            (self['signal'], self['time']) = resample(self['signal'], num, 
                                                      t=self['time'], axis=0, 
                                                      window=None)
        else: 
            self['signal'] = resample(self['signal'], num, axis=0, 
                                      window=None)


class MSTdata:

    """A class that allows one to analyize MST data from the MDSPlus
    tree in a python-friendly environment.

    Usage: Specify an MDSPlus tree and shot number to create the
    object to contain all of the data you want to analyze as well as
    the methods for analysis.
    """

    # Define how the MSTdata object is created. This will make an
    # object which will contain the signals requested in its 'data'
    # member. The signals will also be stored to a local database.

    def __init__(self, shot, tree='mst'):

        self.shot = shot  # Remember the shot number we're looking at
        self.tree = tree
        self.data = {}    # Create a dictionary to store the signals
        self.shelf_path() # Note that this method replaces itself with
                          # a path string
        self.server = get_server_for_shot(shot)
        
        # If there is already locally stored data then read it in.
        if self.shelf_exists(): self.read_shelf()
            
    def __getitem__(self, nodepath=''):

        """The idea here is to provide an easy method for grabbing signal
        data based on a dictionary reference.
        
        signal = MSTdata['subtree::nodepath']
        """        

        # Check if the item is already in our object
        subtree, node = nodepath.lstrip('\\').split('::')
        if self.data.has_key(subtree) and self.data[subtree].has_key(node):
            return self.data[subtree][node]
        else:
            self.get([node], subtree)
            return self.data[subtree][node]

    def __delitem__(self, nodes=None, subtree=None):
        
        """Here we provide a method of removing a signal both from
        local memory and local storage."""

        if subtree is not None:
            self.open_shelf()
            
            if nodes is None:
                del self.local_store[subtree]
                del self.data[subtree]
            else:
                for node in nodes:
                    del self.local_store[subtree][node]
                    del self.data[subtree][node]
            self.close_shelf()

    def __len__(self):
        
        """The length of the MSTdata structure will be the number of
        signals in the data object."""

        return len(self.data)

    def get(self, nodes=[''], subtree='', dims=[],
            local_store=store_data_locally):

        """For each node in the node list grab the signal from either
        the local database (if it's there) or from the MDSPlus
        database."""

        for node in nodes:            
            try: 
                if subtree != '': nodepath = '\\' + subtree + '::' + node
                else: nodepath = node
                if not self.data.has_key(subtree): self.data[subtree] = {}
                name = node.split(':')[-1]
                self.data[subtree][name] = Signal(self.shot, nodepath, dims)
            except:
                print 'Error reading ' + nodepath
                continue

        # Now update the local database
        if local_store: self.store()

    def get_node_list_recursive(self, subtree='top'):

        """Get a list of ALL the available nodes in a subtree."""
        node_list = []

        try: 
            conn = get_connection(shot)
            mdsplus_expr = 'getnci("\\{0}***","NODE_NAME")'.format(subtree)
            node_list = conn.get(mdsplus_expr)

            # Strip out the extra whitespace in the node names
            node_list = [ node.strip() for node in node_list ]

        except:
            print 'MDSPlus database currently unavailable.'

        return node_list

    def get_node_list(self, subtree='\\top'):

        """Get a list of the available members in the first level of a
        subtree. Defaults to \top."""
        
        node_list = []

        try:
            conn = get_connection(shot, tree)
            mdsplus_expr = 'getnci("\\{0}...","NODE_NAME")'.format(subtree)
            node_list = conn.get(mdsplus_expr)

            # Strip out the extra whitespace in the node names
            node_list = [ node.strip() for node in node_list ]

        except:
            print 'MDSPlus database currently unavailable.'

        return node_list

    def common_time_base(self, subtree=None, nodes=None, n_time=inf, 
                         window=None):

        """Put all of the signals in the list of tree nodes on a
        common time base. By default, the timebase of the lowest
        resolution signal is used."""

        if nodes is None: nodes = self.data[subtree].keys()

        # Determine the signal with the lowest resolution
        if dt is inf:
            for sig in nodes:
                n_time = min( [n_time, len(self.data[subtree][sig].time)] )

        for sig in nodes:
            self.data[subtree][sig].resample(n_time, window=window)


    def plot(self, subtree=None, nodes = None):
        
        """Plot some signals from the subtree in the dataset."""

        if subtree is None:
            subtree = self.data.iterkeys().next()
        if nodes is None:
            nodes = self.data[subtree].keys()
        n_nodes = len(nodes)
        p.figure()
        p.subplots_adjust(hspace=0.1)
        ax = {}
        
        for node in nodes:
            i = nodes.index(node)+1
            if i>1: ax[i] = p.subplot(n_nodes,1,i,sharex=ax[1])
            else: ax[i] = p.subplot(n_nodes,1,i)
            self.data[subtree][node].plot()
            if i<n_nodes: 
                p.setp(ax[i].get_xticklabels(),visible=False)
            else:
                p.xlabel('Time [s]')
            p.legend()

#    def 2Dplot(self, signal_list=[], x=None, y=None):
#        """Create a surface plot of a list of signals."""
#        pass

# Next we start a section that allows us to store the data locally
# using the shelve module. We will create a local path to the data

    def shelf_path(self):

        """Establish the location of the locally stored data."""

        self.shelf_path = path.join(local_path,str(self.shot))

    def shelf_exists(self):

        """Determine whether there is a locally stored version of the data."""

        return path.isfile(self.shelf_path)

    def open_shelf(self):
        
        """Open a locally stored version of the data."""

        self.local_store = shelve.open(self.shelf_path)
        
    def read_shelf(self):

        """Read in the locally stored data."""

        if self.shelf_exists():
            self.open_shelf()
            for subtree in self.local_store:
                self.data[subtree] = self.local_store[subtree]
            self.close_shelf()

    def close_shelf(self):

        """Close the locally stored version of the data."""

        self.local_store.close()

    def store(self, subtrees=None):
        
        """Store the data in a local file."""

        self.open_shelf()
        if subtrees is None:
            for subtree in self.data:
                self.local_store[subtree] = self.data[subtree]
        elif type(subtrees) is type([]):
                for subtree in subtrees:
                    self.local_store[subtree] = self.data[subtree]
        elif type(subtrees) is type(''):
                self.local_store[subtrees] = self.data[subtrees]
        self.close_shelf()

    def clear_shelf(self):
        
        """Remove locally stored data."""

        if self.shelf_exists():
            remove(self.shelf_path)

# Now we provide a series of methods to obtain common groups of signals.

    def get_proc_ops(self):

        """Get all the signals in \mst_ops:proc_ops."""

        subtree = 'mst_ops.proc_ops'
        node_list = self.get_node_list_recursive(subtree)
        self.get(node_list, 'mst_ops')

    def get_magnetic_modes(self, 
                           directions = [ 'bp', 'bt' ], 
                           modes = range(1,16), 
                           names = [ 'amp', 'vel', 'phs' ] ):

        """Read in the magnetic mode data from the pickup coil array.
        By default, the routine reads in amplitude, velocity, and
        phase for all of the modes for both poloidal and toroidal
        fields, but you can specify which direction, mode, and data to
        read by using the keywords:
        
        directions = [ 'bp', 'bt' ] (default)
        modes = range(1,16) (default)
        names = [ 'amp', 'vel', 'phs' ] (default)

        where 'amp' is mode amplitude in gauss, 'vel' is mode velocity
        in km/s, and 'phs' is the mode phase in radians.
        """

        # The nodes at which particular magnetic toroidal harmonics
        # are stored have the following naming convention:
        #   \(direction)_n(mode)_(name) where 
        #     direction is either bt or bp,
        #     mode is the mode number,
        #     and name is amp, vel, or phase
        #
        # e.g. \bp_n06_vel is the velocity of the n=6 mode of the
        # poloidal field.

        # Here we construct a list of the node names of the various
        # magnetic mode data.

        subtree = 'mst_mag'
        nodes = []
        
        for direction in directions:
            for mode in modes:
                for name in names:
                    nodes.append('%(direction)s_n%(#)02d_%(name)s' % \
                                     {'direction': direction,
                                      '#': mode,
                                      'name': name})

        # And here we retreive the data from the server.
        self.get(nodes, subtree)


    def get_ids(self):
        
        """Read in the processed IDS data. Output includes:
        x_lam: The assumed (or calibrated) central line wavelength.
        species: The impurity species used for the analysis.
        ids_amp: A timeseries of the impurity emission amplitude.
        ids_vel: The fitted impurity velocity timeseries in km/s.
        ids_temp: The fitted impurity temperature in eV.
        ids_chi2: A timeseries of the squared residuals.
        """

        subtree = 'mst_ids'
        nodes = ['x_lam', 'species', 'amp_17to32', 'vel_17to32', 
                 'temp_17to_32', 'chi2_17to32' ]

        self.get(nodes, subtree)

    def get_ts(self):

        """Read in the Thompson Scattering temperature data."""

        subtree = 'mst_tsmp'
        node = ['top.proc_tsmp.t_e']
        self.get(node, subtree, dims=['time', 'r'])
    
class Ensemble(dict):

    """A class that creates a grouping of common signals from a shot
    list to perform operations like calculating distribution
    functions, timeseries averaging, etc. The ensemble does not
    contain the raw data -- that is the job of the MSTdata class --
    only information on the ensemble."""

    def __init__(self, name='', shotlist=[], nodelist=[], condition=None):

        self.name = name
        self.shotlist = shotlist
        self.nodelist = nodelist
        self.data = {}
        self.condition = slice(0,Ellipsis)   # Default to using it all
        
        # Now we will implement the specified condition. In the
        # following set of cases, the condition will specify what
        # elements of the timeseries signals are of interest for
        # analysis.
        if condition is None:
            shot = self.shotlist[0]
            current = MSTdata(shot)['\\mst_ops::ip']['signal']
            self.condition = slice(0,len(current))
        if condition is 'flattop':
            for shot in shotlist:
                # The current flattop will be defined as the period
                # when the measured plasma current is within 2.5% of its
                # maximum
                current = MSTdata(shot)['\\mst_ops::ip']
                curr = current['signal'].max()
                curr_array = where((abs((current['signal'] -
                                         curr))/curr) < 0.05)
                self.condition = slice(max(curr_array[0], 
                                           self.condition.start), 
                                       min(curr_array[-1], 
                                           self.condition.stop))
                self.start_time = current['time'][self.condition.start]
                self.stop_time = current['time'][self.condition.stop]

    def __getitem__(self, node):

        # If the ensemble average has already been calculated then it
        # will be in the data dictionary
        if self.data.has_key(node): return self.data[node]

        # Grab the first shot off the shotlist to determine the length
        # of the arrays and to initialize the average, minimum, and
        # maximum arrays.
        shotlist = list(self.shotlist)
        shot = shotlist.pop()
        shot_data = MSTdata(shot)
        sig = shot_data[node]
        time_window = (sig['time'] >= self.start_time) & \
            (sig['time'] <= self.stop_time) 
        # Provide a clean exit incase there is no time overlap with
        # the ensemble condition
        if time_window.any() is False:
            print 'Data are unavailable for ' + node + \
                ' during condition used in ensemble.'
            return

        sig['time'] = sig['time'][time_window]
        sig['signal'] = sig['signal'][time_window]
        sig_sum = sig['signal']
        min_sig = sig['signal']
        max_sig = sig['signal']

        for shot in shotlist:
            shot_data = MSTdata(shot)
            sig = shot_data[node]
            time_window = (sig['time'] >= self.start_time) & \
                (sig['time'] <= self.stop_time) 
            sig['time'] = sig['time'][time_window]
            sig['signal'] = sig['signal'][time_window]
            min_sig = array([min_sig, sig['signal']]).min(axis=0)
            max_sig = array([max_sig, sig['signal']]).max(axis=0)
            sig_sum += sig['signal']

        # What we store and return is an enemble-averaged signal along
        # with the range of possible values
        sig['signal'] = sig_sum / len(self.shotlist)
        sig['max'] = max_sig
        sig['min'] = min_sig
        self.data[node] = sig

        return sig


    def __len__(self):

        """We define the length of an ensemble based on the length of
        the condition slice."""

        return self.condition.stop - self.condition.start


    def calc_distribution(self, **keywords):

        """Estimate the distribution function for a signal by creating
        a normalized histogram. Additional keywords are passed to
        scipy's 'histogram' function."""
        
        keywords['density'] = True
        if not keywords.has_key('bins'): 
            n_bins = 51
        else:
            if type(keywords['bins']) is int: n_bins = keywords['bins']
            else: n_bins = len(keywords['bins']) - 1

        # Make a copy of the shotlist so that we can pop off a single
        # shot to initialize the histogram arrays
        shotlist = list(self.shotlist)

        # Initialize the arrays of distributions and the bins used for
        # each signal by looking at just one shot in the shotlist
        shot = shotlist.pop()
        shot_data = MSTdata(shot)
        for node in self.nodelist:
            time_window = (shot_data[node]['time'] >= self.start_time) & \
                (shot_data[node]['time'] <= self.stop_time)
            hist, bin = histogram(shot_data[node]['signal'][time_window], 
                                  **keywords)
            self[node]['distribution'] = array([hist, bin])

        # Now proceed with the rest of the shots in the shotlist
        for shot in shotlist:
            shot_data = MSTdata(shot)
            for node in self.nodelist:
                time_window = (shot_data[node]['time'] >= self.start_time) & \
                    (shot_data[node]['time'] <= self.stop_time)
                keywords['bins'] = self[node]['distribution'][1]
                hist, bin = histogram(shot_data[node]['signal'][time_window],
                                      **keywords)
                self[node]['distribution'][0] += hist

        # Now normalize the distributions
        for node in self.nodelist:
            hist = self[node]['distribution'][0]
            hist /= hist.sum()

    def plot_distribution(self, node, **keywords):
        
        """Once the distributions have been calculated they can be
        plotted as a probability distribution curve."""

        if node in self.data and 'distribution' in self[node]:

            dist, bins = self[node]['distribution']
            height = 0.75*(bins[1]-bins[0])

            p.figure(figsize=(12,5))
            ax1 = p.subplot(121)
            self[node].plot()
            ylim = ax1.get_ylim()
            ax2 = p.subplot(122)
            p.barh(bins[0:-1], dist, height=height, label=self[node]['name'])
            p.xlabel('PDF')
            ax2.set_ylim(ylim)
 
# class MSTfit:

#     """A class that allows one to analyize MST data from an MSTfit
#     IDL save file."""

#     def __init__(filename = None):

#         try:
#             self.fit = idl.readsav(filename)
#         except:
#             print 'Can not open ' + filename

