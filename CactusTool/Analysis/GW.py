from numpy.core.einsumfunc import _parse_einsum_input
from . import TimeSeries
import numpy as np

class GravitationalWave(TimeSeries):
    """
    Class describing discrete 1d wave functions
    """

    def __init__(self, t, y, CU=True):
        """
        Initialize the `GravitationalWave` class.

        :param t: time series
        :type t: array_like
        :param y: values
        :type y: array_like
        """
        super().__init__(t, y, CU)

    @property
    def frequency(self):
        """Returns the instantenous frequency $\omega$ as a DiscreteFunction object

        arXiv:1006.1632
        """
        return self.phase.gradient()

    # def ShiftToZero(self):
    #     idx = self.abs().y.argmax()
    #     return self.shift(-self.t[idx])
        

class Match:
    """
    Return the match between the two waveforms
    """
    def __init__(self, WF1, WF2, psd, low_frequency_cutoff=None, high_frequency_cutoff=None):
        """AI is creating summary for __init__

        :param WF1: [description]
        :type WF1: [type]
        :param WF2: [description]
        :type WF2: [type]
        :param psd: A power spectral density to weight the overlap.
        :type psd: Frequency Series
        :param low_frequency_cutoff: The frequency to begin the match.
        :type low_frequency_cutoff: {None, float}, optional
        :param high_frequency_cutoff: The frequency to stop the match.
        :type high_frequency_cutoff: {None, float}, optional
        """
        self.WF1 = WF1
        self.WF2 = WF2
        self.psd = psd