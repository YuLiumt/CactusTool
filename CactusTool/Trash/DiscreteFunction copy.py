from ..utils.units import UnitConversion
from ..Lib.kuibit import BaseNumerical
from scipy import integrate, interpolate
import numpy as np
import pandas as pd


class BaseSeries(BaseNumerical):
    """
    Base class for generic series data.
    """
    def __init__(self, t, y, CU=True):
        """Initialise a generic series.

        :param array_like t: Sampling times, need to be strictly increasing.
        :param array_like y: values
        :param bool CU: Cactus Units
        """
        assert len(t) >= 1, "initial_array must contain at least one sample."
        assert len(t) == len(y), "Times and Values length mismatch"
        assert np.diff(t).min() >= 0, "Timeseries is not monotonically increasing"

        self.t = np.array(t)
        self.y = np.array(y)

        # Cactus Unit
        self._CU = CU

        self.is_complex = issubclass(self.y.dtype.type, complex)

    def _apply_binary(self, other, function):
        """
        This is an abstract function that is used to implement mathematical operations with other series (if they have the same t).

        :param other: Other TimeSeries.
        :param callable function: Function to apply to the series.
        """
        # If the other object is of the same type
        if isinstance(other, type(self)):
            if (len(self.t) != len(other.t)) or (not np.allclose(other.t, self.t, atol=1e-14)):
                raise ValueError("The objects do not have the same t!")
            assert self._CU == other._CU, "Unit is different"
            return self.__class__(self.t, function(self.y, other.y), self._CU)

        # If it is a number
        if isinstance(other, (int, float, complex)):
            return self.__class__(self.t, function(self.y, other), self._CU)

        # If we are here, it is because we cannot add the two objects
        raise TypeError("I don't know how to combine these objects")

    def __eq__(self, other):
        """Check for equality up to numerical precision."""
        if isinstance(other, type(self)):
            return np.allclose(self.t, other.t, atol=1e-14) and np.allclose(self.y, other.y, atol=1e-14)
        return False

    def _apply_unary(self, function):
        """Apply a unary function to the data.

        :param callable function: Function to apply to the series.
        """
        return self.__class__(self.t, function(self.y), self._CU)

    def _apply_reduction(self, reduction):
        """Apply a reduction to the data.

        :param callable function: Function to apply to the series.
        """
        return reduction(self.y)

    def __len__(self):
        """The number of data points."""
        return len(self.t)

    def __iter__(self):
        for x, y in zip(self.t, self.y):
            yield x, y

    @cached_property
    def tstart(self):
        return self.t[0]

    @cached_property
    def tend(self):
        return self.t[-1]

    @cached_property
    def duration(self):
        """
        Return the length of the covered time interval.
        """
        return self.tend - self.tstart

    @cached_property
    def dt(self):
        dt = np.diff(self.t)
        assert np.allclose(dt, dt[0], atol=1e-14), "Timeseries is not is regularly sampled"
        return dt[0]

    def NearestTimeIndex(self, v):
        """Find nearest neighboring point (index) to t

        :param float t: nearest time
        """
        return (np.abs(self.t-v)).argmin()

    def shift(self, tshift):
        return self.__class__(self.t + tshift, self.y, self._CU)  

    def redshift(self, z):
        """Apply redshift to the data by rescaling the time so that the frequencies are redshifted by ``1 + z``.
        """
        factor = 1 / (1 + z)
        return self.__class__(self.t*factor, self.y, False)

    def crop(self, tmin, tmax):
        """
        Remove data outside the the interval ``[init, end]``. If ``init`` or ``end`` are not specified or None, it does not remove anything from this side.

        :param float tmin: Left boundary cut interval
        :param float tmax: Right boundary cut interval
        """
        t = self.t
        index = np.where((t >= tmin) & (t <= tmax))
        return self.__class__(t[index], self.y[index], self._CU)
    
    # Define aliases
    clip = crop

    def to(self, unitt, unity):
        assert self._CU, "Only from Cactus Units to other"
        t = self.t * UnitConversion(unitt)
        y = self.y * UnitConversion(unity)
        return self.__class__(t, y, False)     

class TimeSeries(BaseNumerical):
    """
    Class describing discrete functions
    """

    def __init__(self, t, y, CU=True):
        """Initialise a 1D transformation function

        :param t: time series
        :type t: array_like
        :param y: values
        :type y: array_like
        """
        if len(y) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        assert len(t) == len(y), "Times and Values length mismatch"

        self.t = np.array(t)
        self.y = np.array(y)

        # Cactus Unit
        self._CU = CU

    @property
    def tstart(self):
        return self.t[0]

    @property
    def tend(self):
        return self.t[-1]

    @property
    def dt(self):
        return self.t[1] - self.t[0]

    def __len__(self):
        return len(self.t)

    def __getitem__(self, idx):
        return self.y[idx]

    def to(self, unitt, unity):
        assert self._CU, "Only from CU to other"
        t = self.t * UnitConversion(unitt)
        y = self.y * UnitConversion(unity)
        return self.__class__(t, y, False)     

    def shift(self, tshift):
        return self.__class__(self.t + tshift, self.y)     

    def remove_mean(self):
        """
        Remove the mean value from the data.
        """
        f = self.y - self.y.mean()
        return self.__class__(self.t, f)

    def clip(self, tmin, tmax):
        """Throws away data outside the time intarval [tmin, tmax].

        :param tmin: Left boundary cut interval
        :type tmin: float
        :param tmax: Right boundary cut interval
        :type tmax: float
        """
        t = self.t
        index = np.where((t >= tmin) & (t <= tmax))
        return self.__class__(t[index], self.y[index])

    def smooth(self, window_len=11, window='hanning'):
        """Smooth the data by convoluting with a window function.

        :param window_len: Smoothing length, defaults to 11
        :type window_len: int, optional
        :param window: The window function, defaults to 'hanning'
        :type window: str, optional
        """
        assert window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman'], "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        if window_len < 3:
            return self.__class__(self.t, self.y)
        assert window_len % 2 == 1, "windowLen should be an odd integer"

        x = self.y
        assert x.size >= window_len, "Input data needs to be bigger than window size."
        s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]

        if window == 'flat':
            w = np.ones(window_len, 'd')
        else:
            w = eval('np.'+window+'(window_len)')
        
        y = np.convolve(w / w.sum(), s, mode='valid')
        return self.__class__(self.t, y[(window_len // 2):-(window_len // 2)])

    def NearestTimeIndex(self, v):
        """Find nearest neighboring point (index) to t

        :param t: nearest time
        :type t: float
        """
        return (np.abs(self.t-v)).argmin()

    def NearestValueIndex(self, v):
        """Find nearest neighboring point (index) to y

        :param t: nearest time
        :type t: float
        """
        return (np.abs(self.y-v)).argmin()

    def resample(self, t, **kwargs):
        """Resample field to new timeseries using :py:class:`scipy.interpolate.interp1d`. Default is 'cubic'.

        :param t: timeseries
        :type t: array_like
        """
        if 'kind' not in kwargs:
            kwargs['kind'] = 'cubic'
        f = interpolate.interp1d(self.t, self.y, **kwargs)
        return self.__class__(t, f(t))

    def cumtrapz(self, **kwargs):
        """
        :return: Cumulatively integrate is computed using :py:class:`scipy.integrate.cumtrapz`.
        :rtype: :py:class:`DiscreteFunction`
        """
        if 'initial' not in kwargs:
            kwargs['initial'] = 0
        f = integrate.cumtrapz(self.y, self.t, **kwargs)
        return self.__class__(self.t, f)

    def gradient(self, **kwargs):
        """
        :return: The gradient is computed using :py:class:`numpy.gradient` (second order accurate central differences in the interior points and either first or second order accurate one-sides differences at the boundaries).
        :rtype: :py:class:`DiscreteFunction`
        """
        f = np.gradient(self.y, self.t, **kwargs)
        return self.__class__(self.t, f)
        
    def polyfit(self, deg, **kwargs):
        return np.polyfit(self.t, self.y, deg, **kwargs)

    def psd(self):
        """ Calculate the power spectral density of this time series.

        Use the `pycbc.psd.welch` method to estimate the psd of this time segment.
        For more complete options, please see that function.

        Parameters
        ----------
        segment_duration: float
            Duration in seconds to use for each sample of the spectrum.
        kwds : keywords
            Additional keyword arguments are passed on to the `pycbc.psd.welch` method.

        Returns
        -------
        psd : FrequencySeries
            Frequency series containing the estimated PSD.
        """
        pass

    def preview(self, mode='lines+markers'):
        """
        Plot by :py:class:`plotly.graph_objects.go`
        """
        try:
            import plotly.graph_objects as go
        except:
            print("Package plotly not found. Install plotly to get nice plot.")

        fig = go.Figure(data=go.Scatter(
            x=self.t, 
            y=self.y,
            mode=mode,
        ))
        fig.show()

    def save(self, fname):
        """
        Saves into simple ascii format with 2 collumns (t,y) for real valued data.

        :param str fname: File name.
        """
        np.savetxt(fname, np.transpose((self.t, self.y)))


class ComplexSeries(TimeSeries):
    """
    Class describing discrete complex functions
    """

    def __init__(self, t, y, CU=True):
        """
        Initialise a 1D transformation function

        :param t: time series
        :type t: array_like
        :param y: complex number
        :type y: array_like
        """
        assert np.iscomplex(y).any(), "Discrete functions is not a complex number."
        super().__init__(t, y, CU)

    def __len__(self):
        return len(self.t)

    def __getitem__(self, idx):
        return self.y[idx]

    def to(self, unitt, unity):
        assert self._CU, "Only from CU to other"
        t = self.t * UnitConversion(unitt)
        y = self.y * UnitConversion(unity)
        return self.__class__(t, y, False) 

    @property
    def real(self):
        return TimeSeries(self.t, self.y.real)

    @property
    def imag(self):
        return TimeSeries(self.t, self.y.imag)

    @property
    def conjugate(self):
        return self.__class__(self.t, self.y.conjugate())

    @property
    def absolute(self):      
        return TimeSeries(self.t, np.absolute(self.y))

    @property
    def phase(self):
        return TimeSeries(self.t, np.unwrap(np.angle(self.y)))

    def preview(self, mode='lines+markers'):
        """
        Plot by :py:class:`plotly.graph_objects.go`
        """
        try:
            import plotly.graph_objects as go
        except:
            print("Package plotly not found. Install plotly to get nice plot.")

        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=self.t, 
            y=self.y.real,
            mode=mode,
            name ='real',
        ))
        fig.add_trace(go.Scatter(
            x=self.t, 
            y=self.y.imag,
            mode=mode,
            name='imag',
        ))

        fig.show()

    def save(self, fname):
        """
        Saves into simple ascii format with 3 collumns (t, Re(y), Im(y)) for complex valued data.

        :param str fname: File name.
        """
        np.savetxt(fname, np.transpose((self.t, self.y.real, self.y.imag)))


class VectorSeries:
    """
    Class describing discrete vector functions
    """

    def __init__(self, t ,y, CU=True):
        self._vnum = y.shape[1]
        self.t = t
        self.y = y

        # Cactus Unit
        self._CU = CU

    @property
    def tstart(self):
        return self.t[0]

    @property
    def tend(self):
        return self.t[-1]

    @property
    def dt(self):
        return self.t[1] - self.t[0]

    def __len__(self):
        return len(self.t)

    def __getitem__(self, idx):
        return self.y[idx]

    def to(self, unitt, unity):
        assert self._CU, "Only from CU to other"
        t = self.t * UnitConversion(unitt)
        y = self.y * UnitConversion(unity)
        return self.__class__(t, y, False) 

    @property
    def dim(self):
        return self.y.shape[1]

    def get_component(self, n):
        return TimeSeries(self.t, self.y[:, n])

    def shift(self, tshift):
        return self.__class__(self.t + tshift, self.y)   

    def remove_mean(self):
        """
        Remove the mean value from the data.
        """
        f = self.y - self.y.mean(axis=0)
        return self.__class__(self.t, f)

    def clip(self, tmin, tmax):
        """Throws away data outside the time intarval [tmin, tmax].

        :param tmin: Left boundary cut interval
        :type tmin: float
        :param tmax: Right boundary cut interval
        :type tmax: float
        """
        t = self.t
        index = np.where((t >= tmin) & (t <= tmax))
        return self.__class__(t[index], self.y[index])

    def NearestTimeIndex(self, v):
        """Find nearest neighboring point (index) to t

        :param t: nearest time
        :type t: float
        """
        return (np.abs(self.t-v)).argmin()

    def resample(self, t, **kwargs):
        """Resample field to new timeseries using :py:class:`scipy.interpolate.interp1d`. Default is 'cubic'.

        :param t: timeseries
        :type t: array_like
        """
        if 'kind' not in kwargs:
            kwargs['kind'] = 'cubic'

        foft = []
        for i in range(self._vnum):
            f = interpolate.interp1d(self.t, self.y[:, i], **kwargs)
            foft.append(f(t))
        return self.__class__(t, np.stack(foft, axis=1))

    def cumtrapz(self, **kwargs):
        """
        :return: Cumulatively integrate is computed using :py:class:`scipy.integrate.cumtrapz`.
        :rtype: :py:class:`DiscreteFunction`
        """
        if 'initial' not in kwargs:
            kwargs['initial'] = 0

        foft = []
        for i in range(self._vnum):
            foft.append(integrate.cumtrapz(self.y[:, i], self.t, **kwargs))
        return self.__class__(self.t, np.stack(foft, axis=1))

    def gradient(self, **kwargs):
        """
        :return: The gradient is computed using :py:class:`numpy.gradient` (second order accurate central differences in the interior points and either first or second order accurate one-sides differences at the boundaries).
        :rtype: :py:class:`DiscreteFunction`
        """
        foft = []
        for i in range(self._vnum):
            foft.append(np.gradient(self.y[:, i], self.t, **kwargs))
        return self.__class__(self.t, np.stack(foft, axis=1))

    @property
    def norm(self):
        f = np.linalg.norm(self.y, axis=1)
        return TimeSeries(self.t, f)

    def preview(self, mode='lines+markers'):
        """
        Plot by :py:class:`plotly.graph_objects.go`
        """
        try:
            import plotly.graph_objects as go
        except:
            print("Package plotly not found. Install plotly to get nice plot.")
        
        y = self.y

        fig = go.Figure()

        for i in range(self._vnum):
            fig.add_trace(go.Scatter(
                x=self.t, 
                y=self.y[:, i],
                mode=mode,
                name ='{} component'.format(i),
            ))

        fig.show()

    def save(self, fname):
        """
        Saves into simple ascii format with 2 collumns (t,y) for real valued data.

        :param str fname: File name.
        """
        np.savetxt(fname, np.vstack((self.t, self.y.T)).T)

class DataSeries:
    """
    Class describing discrete DataFrame
    """

    def __init__(self, t, y, columns, CU=True):
        self.t = t
        self.y = y
        self.columns = columns

        # Cactus Unit
        self._CU = CU

    @property
    def tstart(self):
        return self.t[0]

    @property
    def tend(self):
        return self.t[-1]

    @property
    def dt(self):
        return self.t[1] - self.t[0]

    def __len__(self):
        return len(self.t)

    def __getitem__(self, idx):
        return self.y[idx]

    def to(self, unitt, unity):
        assert self._CU, "Only from CU to other"
        t = self.t * UnitConversion(unitt)
        y = self.y * UnitConversion(unity)
        return self.__class__(t, y, self.columns, False) 

    @property
    def DataFrame(self):
        return pd.DataFrame(self.y.T, index=self.t, columns=self.columns)

class FrequencySeries:
    """
    Models a frequency series consisting of uniformly sampled scalar values.
    """
    def __init__(self, f, y):
        """AI is creating summary for __init__

        :param f: [description]
        :type f: [type]
        :param y: [description]
        :type y: [type]
        :raises ValueError: [description]
        """
        if len(y) < 1:
            raise ValueError('initial_array must contain at least one sample.')
        assert len(t) == len(y), "Frequency and Values length mismatch"
        self.f = f
        self.y = y

    # @property
    # def tstart(self):
    #     return self.t[0]

    # @property
    # def tend(self):
    #     return self.t[-1]

    @property
    def df(self):
        return self.f[1] - self.f[0]

    def to_timeseries(self, delta_t=None):
        """ Return the Fourier transform of this time series.

        Note that this assumes even length time series!

        Parameters
        ----------
        delta_t : {None, float}, optional
            The time resolution of the returned series. By default the
        resolution is determined by length and delta_f of this frequency
        series.

        Returns
        -------
        TimeSeries:
            The inverse fourier transform of this frequency series.
        """
        pass