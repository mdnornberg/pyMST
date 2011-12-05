import pylab as p
import scipy.signal

def where(condition):

    """Return the indices of an array that satisfies a condition.

    Input:
       condition (bool): Some comparison which generates a boolean array.

    Output:
       (array of ints): An array of indices where the condition is true."""

    return condition.nonzero()[0]

def close_to(array, value):

    """Return the indices of the array element closest to a given value.

    Input:
       array: An array of numbers.
       value: A number.

    Output: An array of indices where the array value is closest to the
       input value."""

    return where(abs(array - value) == abs(array - value).min())

def plot_band(x, y, yerr=None, ymin=None, ymax=None, alpha=0.25,
              **keywords):

    """Plots a shaded region to represent the uncertainty in time
    series data."""

    if 'alpha' not in keywords: keywords['alpha']=0.25
    if ymin is not None and ymax is not None:
            x0,y0 = p.poly_between(x,ymin,ymax)
    else:
        if yerr==None: yerr = 0.05*y
        x0,y0 = p.poly_between(x,y-y_err,y+y_err)
    p.fill(x0,y0,**keywords)

    return

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=scipy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=scipy.ones(window_len,'d')
    else:
        w=eval('scipy.'+window+'(window_len)')

    y=scipy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]
