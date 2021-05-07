from ..utils.log import logger
from .DiscreteFunction import TimeSeries, VectorSeries
import numpy as np

class BlackHoleBinary:

    def __init__(self, t, a1, a2, a3, m1, m2, m3, s1, s2, s3, Sx1, Sy1, Sz1, Sx2, Sy2, Sz2, Sx3, Sy3, Sz3):
        self.t = t
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
        self.Sx1 = Sx1
        self.Sy1 = Sy1
        self.Sz1 = Sz1
        self.Sx2 = Sx2
        self.Sy2 = Sy2
        self.Sz2 = Sz2
        self.Sx3 = Sx3
        self.Sy3 = Sy3
        self.Sz3 = Sz3


    @property
    def area(self):
        return VectorSeries(self.t, np.stack([self.a1, self.a2, self.a3], axis=1))

    @property
    def spin(self):
        return VectorSeries(self.t, np.stack([self.s1, self.s2, self.s3], axis=1))

    @property
    def M(self):
        return VectorSeries(self.t, np.stack([self.m1, self.m2, self.m3], axis=1))

    @property
    def Spin1(self):
        return VectorSeries(self.t, np.stack([self.Sx1, self.Sy1, self.Sz1], axis=1))

    @property
    def Spin2(self):
        return VectorSeries(self.t, np.stack([self.Sx2, self.Sy2, self.Sz2], axis=1))

    @property
    def Spin3(self):
        return VectorSeries(self.t, np.stack([self.Sx3, self.Sy3, self.Sz3], axis=1))

    @property
    def tmerger(self):
        idx = np.where(self.m3 != 0)[0][0]
        return self.t[idx]

    @property
    def mf(self):
        """
        remnant mass.
        """
        tmerger = self.tmerger
        tend = self.t[-1]
        tstart = tmerger + 2/3*(tend - tmerger)
        BHf = TimeSeries(self.t, self.m3).clip(tstart, tend)
        mf, res, rank, singular, threshold = BHf.polyfit(0, full=True)
        if res > 1e-2:
            logger.warning("Residuals of the least-squares fit may be too large: {}",format(res))
        return mf[0]

    @property
    def spinf(self):
        """
        remnant spin.
        """
        tmerger = self.tmerger
        tend = self.t[-1]
        tstart = tmerger + 2/3*(tend - tmerger)
        BHf = self.Spin3.clip(tstart, tend)
        spinf = []
        for i in range(BHf.dim):
            spini, res, rank, singular, threshold = BHf.get_component(i).polyfit(0, full=True)
            if res > 1e-2:
                logger.warning("Residuals of the least-squares fit may be too large: {}",format(res))
            spinf.append(spini[0])
        return spinf

    @property
    def chif(self):
        """
        remnant chi.
        """
        return self.spinf / self.mf**2


class BHDiags:
    """
    This class collects the information BH_Diagnostic files saved by the thorn AHFinderDirect.

    :ivar time:        cctk_time
    :ivar pos_x:       x-position of spherical surface center.
    :ivar pos_y:       y-position of spherical surface center.
    :ivar pos_z:       z-position of spherical surface center.
    :ivar rmean:       Mean coordinate radius.
    :ivar r_circ_xy:   Circumferential radius in xy plane.
    :ivar r_circ_xz:   Circumferential radius in xz plane.
    :ivar r_circ_yz:   Circumferential radius in yz plane.
    :ivar m_irr:       Irreducible BH mass.
    """
    def __init__(self, dsets):
        self.dsets = dsets
    
    def __len__(self):
        return len(self.dsets)

    def __iter__(self):
        return iter(self.dsets)
    
    def position(self, ahidx):
        assert ahidx in self.dsets
        dset = self.dsets[ahidx]
        return VectorSeries(dset[0], np.stack([dset[1], dset[2], dset[3]], axis=1))

    def M_irr(self, ahidx):
        assert ahidx in self.dsets
        dset = self.dsets[ahidx]
        return TimeSeries(dset[0], dset[-1])