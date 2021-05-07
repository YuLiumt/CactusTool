from ..utils.units import UnitConversion
from ..Lib.kuibit import BaseSeries
from functools import cached_property
from scipy import integrate, interpolate, signal
import numpy as np
import pandas as pd


class TimeSeries(BaseSeries):
    """
    This class represents real or complex valued time series.
    """

    def __init__(self, t, y, CU=True):
        """Initialise a 1D transformation function

        :param array_like t: Sampling times, need to be strictly increasing.
        :param array_like y: values
        :param bool CU: Cactus Units
        """
        super().__init__(t, y)
        # Cactus Unit
        self._CU = CU

    @cached_property
    def t(self):
        """
        Return the time.
        """
        return self.x

    @cached_property
    def tmin(self):
        """
        Return the starting time.
        """
        return self.xmin

    @cached_property
    def tmax(self):
        """
        Return the final time.
        """
        return self.xmax

    @cached_property
    def dt(self):
        """
        Return the timestep if the series is regularly sampled,
        otherwise raise error.
        """
        assert self.is_regularly_sampled(), "Timeseries is not regularly sampled"
        return self.x[1] - self.x[0]

    @property
    def time_length(self):
        """
        Return the length of the covered time interval (tmax - tmin).
        """
        return self.xmax - self.xmin

    duration = time_length

    def NearestTimeIndex(self, v):
        """Find nearest neighboring point (index) to t

        :param float t: nearest time
        """
        return (np.abs(self.x-v)).argmin()

    def NearestValueIndex(self, v):
        """Find nearest neighboring point (index) to y

        :param float v: nearest y
        """
        return (np.abs(self.y-v)).argmin()

    def redshift(self, z):
        """Apply redshift to the data by rescaling the time so that the frequencies are redshifted by ``1 + z``.
        """
        factor = 1 / (1 + z)
        return self.__class__(self.x*factor, self.y, False)

    def align_at_minimum(self):
        """
        Return a new timeseries with absolute minimum at t=0.
        """
        t = self.x[np.argmin(np.abs(self.y))]
        return self.shift(-t)

    def align_at_maximum(self):
        """
        Return a new timeseries with absolute maximum is at t=0.
        """
        t = self.x[np.argmax(np.abs(self.y))]
        return self.shift(-t)

    def phase_shift(self, pshift):
        """
        Shift the complex phase timeseries by ``pshift``. If the signal is real, it is turned complex with phase of ``pshift``.
        """
        return self.__class__(self.x, self.y * np.exp(1j * pshift), self._CU)

    def phase_angular_velocity(self, window_size, order=3):
        """
        Compute the phase angular velocity, i.e. the time derivative of the complex phase.
        """
        ret_value = self.phase.gradient()
        return ret_value.savgol_smooth(window_size, order)

    def phase_frequency(self, window_size, order=3):
        """
        Compute the phase frequency, i.e. the time derivative
        of the complex phase divided by 2 pi.
        """
        return self.phase_angular_velocity(window_size, order=3) / (2 * np.pi)

    # @property
    # def remove_nan(self):
    #     """Filter out nans/infinite values. Return a new series with finite values only.
    #     """
    #     msk = np.isfinite(self.y)
    #     return self.__class__(self.t[msk], self.y[msk], self._CU)

    def remove_mean(self):
        """
        Remove the mean value from the data.
        """
        f = self.y - self.y.mean()
        return self.__class__(self.x, f, self._CU)

    def resample(self, t, **kwargs):
        """
        Return a new series resampled from this to new t using :py:class:`scipy.interpolate.interp1d`. Default is 'cubic'.

        :param array_like t: New timeseries
        """
        if 'kind' not in kwargs:
            kwargs['kind'] = 'cubic'
        f = interpolate.interp1d(self.x, self.y, **kwargs)
        return self.__class__(t, f(t), self._CU)

    def regular_resample(self):
        """
        Resample the timeseries to regularly spaced times,
        with the same number of points.
        """
        t = np.linspace(self.xmin, self.xmax, len(self))
        return self.resample(t)

    def integrate(self, **kwargs):
        """
        Return a series that is the integral computed with method of :py:class:`scipy.integrate.cumtrapz`.
        """
        if 'initial' not in kwargs:
            kwargs['initial'] = 0
        f = integrate.cumtrapz(self.y, self.x, **kwargs)
        return self.__class__(self.x, f, self._CU)

    def gradient(self, **kwargs):
        """
        Return a series that is the derivative of the current one using :py:class:`numpy.gradient` (second order accurate central differences in the interior points and either first or second order accurate one-sides differences at the boundaries).
        """
        f = np.gradient(self.y, self.x, **kwargs)
        return self.__class__(self.x, f, self._CU)

    def smooth(self, window_len=11, window='hanning'):
        """Smooth the data by convoluting with a window function.

        :param int window_len: Smoothing length, defaults to 11
        :param str window: The window function, defaults to 'hanning'
        """
        assert window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman'], "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        if window_len < 3:
            return self.__class__(self.x, self.y)
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
        
    def savgol_smooth(self, window_size, order=3):
        """Return a smoothed series with a Savitzky-Golay filter with window of size ``window_size`` and order ``order``.

        :param int window_size: Number of points of the smoothing window (needs to be odd).
        :param int order: Order of the filter.
        """
        if self.is_complex:
            f = signal.savgol_filter(self.y.real, window_size, order) + 1j * signal.savgol_filter(self.y.imag, window_size, order)
        else:
            f = signal.savgol_filter(self.y, window_size, order)
        return self.__class__(self.x, f, self._CU)

    def window(self, window_function, *args, **kwargs):
        """
        Apply window_function to the data.

        ``window_function`` can be a string with the name of the window
        function that takes as first argument the number of points of the signal.
        """
        assert window_function in ['tukey', 'hamming', 'blackman'], "Window is one of 'tukey', 'hamming', 'blackman'"
        window_array = eval('signal.'+window_function+'(len(self), *args, **kwargs)')
        return self.__class__(self.t, self.y * window_array)

    def polyfit(self, deg, **kwargs):
        return np.polyfit(self.x, self.y, deg, **kwargs)

    def to_FrequencySeries(self):
        """
        Fourier transform of the timeseries. If the signal is not complex, only positive frequencies are kept. It will be resampled before transforming.

        To have meaningful results, you should consider removing the mean and windowing the signal before calling this method!
        """
        regular_ts = self.regular_resample()
        dt = regular_ts.dt

        if self.is_complex:
            frequencies = np.fft.fftfreq(len(regular_ts), d=dt)
            fft = np.fft.fft(regular_ts.y)

            f = np.fft.fftshift(frequencies)
            fft = np.fft.fftshift(fft)
        else:
            # Note the "r"
            f = np.fft.rfftfreq(len(regular_ts), d=dt)
            fft = np.fft.rfft(regular_ts.y)

        return FrequencySeries(f, fft * dt)

    def plot(self, fig, label='y', **kwargs):
        """
        Plot by :py:class:`matplotlib` or :py:class:`plotly.graph_objects.go`.
        """
        if hasattr(fig, 'plot'):
            if self.is_complex:
                fig.plot(self.x, self.y.real, label=f'Re[{label}]', **kwargs)
                fig.plot(self.x, self.y.imag, label=f'Im[{label}]', **kwargs)
            else:
                fig.plot(self.x, self.y, label=label, **kwargs)
        elif hasattr(fig, 'add_trace'):
            import plotly.graph_objects as go
            if self.is_complex:
                fig.add_trace(go.Scatter(
                    x=self.x, 
                    y=self.y.real,
                    name =f'Re[{label}]',
                    **kwargs,
                ))
                fig.add_trace(go.Scatter(
                    x=self.x, 
                    y=self.y.imag,
                    name=f'Im[{label}]',
                    **kwargs,
                ))
            else:
                fig.add_trace(go.Scatter(
                    x=self.x, 
                    y=self.y,
                    name=label, 
                    **kwargs,
                ))
        else:
            raise RuntimeError("fig must be 'plotly.graph_objects.go' or 'matplotlib.axes._subplots.AxesSubplot'")

    def preview(self, mode='lines+markers'):
        """
        Plot by :py:class:`plotly.graph_objects.go`
        """
        try:
            import plotly.graph_objects as go
        except:
            print("Package plotly not found. Install plotly to get nice plot.")

        fig = go.Figure()
        self.plot(fig, mode=mode)
        fig.show()

    def to(self, unitt, unity):
        assert self._CU, "Only from Cactus Units to other"
        t = self.x * UnitConversion(unitt)
        y = self.y * UnitConversion(unity)
        return self.__class__(t, y, False)     

    def copy(self):
        """
        Return a deep copy.
        """
        return self.__class__(self.x, self.y, self._CU)

    def save(self, fname, *args, **kwargs):
        """
        Saves into simple ascii format with 2 collumns (t,y) for real valued data.

        :param str fname: File name.
        """
        if self.is_complex:
            np.savetxt(fname, np.transpose((self.x, self.y.real, self.y.imag)), *args, **kwargs)
        else:
            np.savetxt(fname, np.transpose((self.x, self.y)), *args, **kwargs)


class VectorSeries(BaseSeries):
    """
    Class describing discrete vector functions
    """

    def __init__(self, t ,y, CU=True):
        super().__init__(t, y)
        # Cactus Unit
        self._CU = CU
        self.dim = self.y.shape[1]

    @cached_property
    def t(self):
        """
        Return the time.
        """
        return self.x

    @cached_property
    def tmin(self):
        """
        Return the starting time.
        """
        return self.xmin

    @cached_property
    def tmax(self):
        """
        Return the final time.
        """
        return self.xmax

    @cached_property
    def dt(self):
        """
        Return the timestep if the series is regularly sampled,
        otherwise raise error.
        """
        assert self.is_regularly_sampled(), "Timeseries is not regularly sampled"
        return self.x[1] - self.x[0]

    @property
    def time_length(self):
        """
        Return the length of the covered time interval (tmax - tmin).
        """
        return self.xmax - self.xmin

    duration = time_length

    def NearestTimeIndex(self, v):
        """Find nearest neighboring point (index) to t

        :param float t: nearest time
        """
        return (np.abs(self.x-v)).argmin()

    def remove_mean(self):
        """
        Remove the mean value from the data.
        """
        f = self.y - self.y.mean(axis=0)
        return self.__class__(self.x, f, self._CU)

    def resample(self, t, **kwargs):
        """Resample field to new timeseries using :py:class:`scipy.interpolate.interp1d`. Default is 'cubic'.

        :param t: timeseries
        :type t: array_like
        """
        if 'kind' not in kwargs:
            kwargs['kind'] = 'cubic'

        foft = []
        for i in range(self.dim):
            f = interpolate.interp1d(self.x, self.y[:, i], **kwargs)
            foft.append(f(t))
        return self.__class__(t, np.stack(foft, axis=1), self._CU)

    def regular_resample(self):
        """
        Resample the timeseries to regularly spaced times,
        with the same number of points.
        """
        t = np.linspace(self.xmin, self.xmax, len(self))
        return self.resample(t)

    def integrate(self, **kwargs):
        """
        :return: Cumulatively integrate is computed using :py:class:`scipy.integrate.cumtrapz`.
        :rtype: :py:class:`DiscreteFunction`
        """
        if 'initial' not in kwargs:
            kwargs['initial'] = 0

        foft = []
        for i in range(self.dim):
            foft.append(integrate.cumtrapz(self.y[:, i], self.x, **kwargs))
        return self.__class__(self.x, np.stack(foft, axis=1), self._CU)

    def gradient(self, **kwargs):
        """
        :return: The gradient is computed using :py:class:`numpy.gradient` (second order accurate central differences in the interior points and either first or second order accurate one-sides differences at the boundaries).
        :rtype: :py:class:`DiscreteFunction`
        """
        foft = []
        for i in range(self.dim):
            foft.append(np.gradient(self.y[:, i], self.x, **kwargs))
        return self.__class__(self.x, np.stack(foft, axis=1), self._CU)

    @property
    def norm(self):
        f = np.linalg.norm(self.y, axis=1)
        return TimeSeries(self.x, f, self._CU)

    def __getitem__(self, n):
        assert self.dim > n
        return TimeSeries(self.x, self.y[:, n], self._CU)

    def plot(self, fig, **kwargs):
        """
        Plot by :py:class:`matplotlib` or :py:class:`plotly.graph_objects.go`.
        """
        if hasattr(fig, 'plot'):
            for i in range(self.dim):
                fig.plot(self.x, self.y[:, i], **kwargs, label='{} component'.format(i))
        elif hasattr(fig, 'add_trace'):
            import plotly.graph_objects as go

            for i in range(self.dim):
                fig.add_trace(go.Scatter(
                    x=self.x, 
                    y=self.y[:, i],
                    name ='{} component'.format(i),
                    **kwargs,
                ))
        else:
            raise RuntimeError("fig must be 'plotly.graph_objects.go' or 'matplotlib.axes._subplots.AxesSubplot'")

    def preview(self, mode='lines+markers'):
        """
        Plot by :py:class:`plotly.graph_objects.go`
        """
        try:
            import plotly.graph_objects as go
        except:
            print("Package plotly not found. Install plotly to get nice plot.")

        fig = go.Figure()
        self.plot(fig, mode=mode)
        fig.show()

    def copy(self):
        """
        Return a deep copy.
        """
        return self.__class__(self.x, self.y, self._CU)

    def save(self, fname):
        """
        Saves into simple ascii format with 2 collumns (t,y) for real valued data.

        :param str fname: File name.
        """
        np.savetxt(fname, np.vstack((self.x, self.y.T)).T)

class DataSeries(VectorSeries):
    """
    Class describing discrete DataFrame
    """

    def __init__(self, t, y, columns, CU=True):
        super().__init__(t, y, CU)
        assert self.dim == len(columns)
        self.columns = columns

    def __getitem__(self, var):
        assert var in self.columns
        idx = self.columns.index(var)
        return TimeSeries(self.t, self.y[:, idx], self._CU)

    def preview(self, mode='lines+markers'):
        """
        Plot by :py:class:`plotly.graph_objects.go`
        """
        try:
            import plotly.graph_objects as go
        except:
            print("Package plotly not found. Install plotly to get nice plot.")
        
        fig = go.Figure()

        for i in range(self.dim):
            fig.add_trace(go.Scatter(
                x=self.x, 
                y=self.y[:, i],
                mode=mode,
                name ='{}'.format(self.columns[i]),
            ))

        fig.show()

    def copy(self):
        """
        Return a deep copy.
        """
        return self.__class__(self.x, self.y, self.columns, self._CU)

    @property
    def DataFrame(self):
        return pd.DataFrame(self.y, index=self.x, columns=self.columns)


class FrequencySeries:
    """
    Class representing a Fourier spectrum.
    """
    def __init__(self, f, y, CU=True):
        """
        :param array_like f: Frequencies.
        :param array_like y: 
        """
        assert len(f) >= 1, "initial_array must contain at least one sample."
        assert len(f) == len(y), "Frequencies and Values length mismatch"
        assert np.diff(f).min() >= 0, "Frequencies is not monotonically increasing"
        self.f = f
        self.y = y

        self._CU = CU

    @cached_property
    def fmin(self):
        """
        Return the minimum frequency.
        """
        return self.f[0]

    @cached_property
    def fmax(self):
        """
        Return the maximum frequency.
        """
        return self.f[-1]

    @cached_property
    def frange(self):
        """
        Return the range of frequencies.
        """
        return self.fmax - self.fmin

    @cached_property
    def df(self):
        df = np.diff(self.f)
        assert np.allclose(df, df[0], atol=1e-14), "Timeseries is not is regularly sampled"
        return df[0]

    def __len__(self):
        """The number of data points."""
        return len(self.f)

    def __iter__(self):
        for f, y in zip(self.f, self.y):
            yield f, y

    def low_pass(self, f):
        """Remove frequencies higher or equal than the given.

        :param float f: Frequency above which the series will be zeroed.
        """
        msk = np.abs(self.f) <= f
        return self.__class__(self.f[msk], self.y[msk], self._CU)

    def high_pass(self, f):
        """Remove all the frequencies smaller than f

        :param float f: Frequency below which series will be zeroed.
        """
        msk = np.abs(self.f) >= f
        return self.__class__(self.f[msk], self.y[msk], self._CU)

    def band_pass(self, fmin, fmax):
        """Remove all the frequencies below ``fmin`` and above ``fmax``.

        :param float fmin: Minimum frequency.
        :param float fmax: Maximum frequency.
        """
        ret = self.low_passed(fmax)
        return ret.high_passed(fmin)

#TODO
# def sample_common(series):
#     """Take a list of series and return new ones so that they are all defined on the same points.

#     :param list series: The series to resample or redefine on the common points.
#     """
#     s1, *s_others = series
#     if s1.is_regularly_sampled():
#         for s in s_others:
#             if not (len(s) == len(s1)):
#                 break
#             if not np.allclose(s1.x, s.x, atol=1e-14):
#                 break
#             # This is an else to the for loop
#         else:
#             # We have to copy, otherwise one can accidentally modify input data
#             return [ss.copy() for ss in series]

#     if resample:
#         # Find the series with max xmin
#         s_xmin = max(series, key=lambda x: x.xmin)
#         # Find the series with min xmax
#         s_xmax = min(series, key=lambda x: x.xmax)
#         # Find the series with min number of points
#         s_ns = min(series, key=len)
#         x = np.linspace(s_xmin.xmin, s_xmax.xmax, len(s_ns))
#         return [
#             s.resampled(x, piecewise_constant=piecewise_constant)
#             for s in series
#         ]

#     def float_intersection(array_1, array_2):
#         """Here we find the intersection between the two arrays also
#         considering the floating points.
#         """
#         # https://stackoverflow.com/a/32516182
#         return array_2[
#             np.isclose(array_1[:, None], array_2, atol=1e-14).any(0)
#         ]

#     # Here we find the common intersection between all the x, starting with
#     # the first one
#     x = s1.x
#     for s in s_others:
#         x = float_intersection(x, s.x)
#         if len(x) == 0:
#             raise ValueError("Series do not have any point in common")

#     return [s.resampled(x) for s in series]