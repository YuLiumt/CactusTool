from scipy import interpolate
import numpy as np
from .array import is_monotonically_increasing

class TimeSeries:
    """
    This class represents real or complex valued time series.

    :param t: Sampling times, need to be strictly increasing.
    :param y: Data samples, can be real or complex valued.
    """
    
    def __init__(self, t, y):
        assert is_monotonically_increasing(t), 'times need to be strictly increasing'
        assert len(t) == len(y), 'Times and Values length mismatch'
        self.t = np.array(t).copy()
        self.y = np.array(y).copy()
        self._complex = issubclass(self.y.dtype.type, complex)
        
    def __len__(self):
        return len(self.t)
    
    def tmin(self):
        return self.t[0]
    
    def tmax(self):
        return self.t[-1]
    
    def remove_mean(self):
        """
        Remove the mean value from the data.
        """
        self.y -= self.y.mean()
        
    def shift(self, tshift):
        self.t += tshift
        
    def clip(self, tmin=None, tmax=None):
        """
        Throws away data outside the time intarval [tmin, tmax]. If tmin or tmax are not specified or None, it does not remove anything from this side.

        :param tmin: Left boundary cut interval or None.
        :param tmax: Right boundary cut interval or None.
        """
        if tmin is not None:
            m = self.t >= tmin
            self.t = self.t[m]
            self.y = self.y[m]
        if tmax is not None:
            m = self.t <= tmax
            self.t = self.t[m]
            self.y = self.y[m]
    
    def conjugate(self):
        assert self._complex, 'Not complex-valued'
        self.y = np.conj(self.y)
    
    def absolute(self):
        assert self._complex, 'Not complex-valued'
        self.y = abs(self.y)
        self._complex = False
    
    def phase(self):
        """
        Compute the complex phase of a complex-valued signal such that
        no phase wrap-arounds occur.

        :returns:   Continuous complex phase.
        """
        phase = np.angle(self.y)/(2*np.pi)
        wind = np.zeros_like(phase)
        wind[1:] = np.rint(phase[1:] - phase[:-1])
        wind = np.cumsum(wind)
        self.y = 2*np.pi*(phase - wind)
    
    def real(self):
        assert self._complex, 'Not complex-valued'
        self.y = self.y.real
        self._complex = False

    def imag(self):
        assert self._complex, 'Not complex-valued'
        self.y = self.y.imag
        self._complex = False
    
    def resample(self, times, ext=0):
        """
        Resamples the timeseries to new times.

        :param times: New sample times.
        :param ext: How to handle points outside the time interval.
        :type ext: 0 for extrapolation, 1 for returning zero, 2 for ValueError.
        :returns: Resampled time series.  
        """
        assert is_monotonically_increasing(times), 'times need to be strictly increasing'
        if self._complex:
            spl_real = interpolate.splrep(self.t, self.y.real)
            spl_imag = interpolate.splrep(self.t, self.y.imag)
            real = interpolate.splev(times, spl_real, ext=ext)
            imag = interpolate.splev(times, spl_imag, ext=ext)
            self.y = (real + 1j*imag)
        else:
            spl = interpolate.splrep(self.t, self.y)
            self.y = interpolate.splev(times, spl, ext=ext)
        self.t = np.array(times).copy()
    
    def deriv(self, order):
        """
        Compute derivative of order<=5 using splines.
        
        :param int order: Order of differentiation.
        """
        assert order <= 5, 'Cannot compute differential of order {}'.format(order)
        if order > 3:
            ks = 5
        else:
            ks = 3
            
        if self._complex:
            spl_real = interpolate.splrep(self.t, self.y.real, k=ks)
            spl_imag = interpolate.splrep(self.t, self.y.imag, k=ks)
            td = self.t[2*ks:-2*ks]
            real = interpolate.splev(td, spl_real, der=order)
            imag = interpolate.splev(td, spl_imag, der=order)
            self.y = (real + 1j*imag)
        else:
            spl = interpolate.splrep(self.t, self.y, k=ks)
            td = self.t[2*ks:-2*ks]
            self.y = interpolate.splev(td, spl, der=order)
        self.t = td
    
    def integrate(self, a=None, b=None):
        """
        Compute the definite integral over the interval [a,b] using spline representations. If lower and/or upper bound is not specified, use boundary of the timeseries.

        :param a: Lower integration bound or None.
        :param b: Upper integration bound or None.
        """
        if a is None:
            a = self.tmin() 
        if b is None:
            b = self.tmax() 
            
        if a > b:
            a, b = b, a
            sf = -1
        else:
            sf = 1
        assert ((a >= self.tmin()) and (b <= self.tmax())), 'Integration bounds out of range.'
        
        if self._complex:
            spl_real = interpolate.splrep(self.t, self.y.real)
            spl_imag = interpolate.splrep(self.t, self.y.imag)
            real = interpolate.splint(a, b, spl_real)
            imag = interpolate.splint(a, b, spl_imag)
            return sf*(real + 1j*imag)
        else:
            spl = interpolate.splrep(self.t, self.y)
            return sf*interpolate.splint(a, b, spl)
    
    def copy(self):
        return TimeSeries(self.t, self.y)
    
    def save(self, fname):
        """
        Saves into simple ascii format with 2 collumns (t,y) for real valued
        data and 3 collumns (t, Re(y), Im(y)) for complex valued data.

        :param str fname: File name.
        """
        if self._complex:
            np.savetxt(fname, np.transpose((self.t, self.y.real, self.y.imag)))
        else:
            np.savetxt(fname, np.transpose((self.t, self.y)))