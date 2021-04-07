from .DiscreteFunction import TimeSeries, VectorSeries
import numpy as np

class BlackHoleBinary:

    def __init__(self, t, m1, m2, m3, Sx1, Sy1, Sz1, Sx2, Sy2, Sz2, Sx3, Sy3, Sz3):
        self.t = t
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.Sx1 = Sx1
        self.Sy1 = Sy1
        self.Sz1 = Sz1
        self.Sx2 = Sx2
        self.Sy2 = Sy2
        self.Sz2 = Sz2
        self.Sx3 = Sx3
        self.Sy3 = Sy3
        self.Sz3 = Sz3

    # @property
    # def M(self):
    #     return VectorSeries(self.t, np.stack([self.m1, self.m2, self.m3], axis=1))

    @property
    def Spin1(self):
        return VectorSeries(self.t, np.stack([self.Sx1, self.Sy1, self.Sz1], axis=1))

    @property
    def Spin2(self):
        return VectorSeries(self.t, np.stack([self.Sx2, self.Sy2, self.Sz2], axis=1))

    @property
    def Spin3(self):
        return VectorSeries(self.t, np.stack([self.Sx3, self.Sy3, self.Sz3], axis=1))