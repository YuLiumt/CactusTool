# Copyright (C) 2020-2021 Gabriele Bozzola
from abc import ABC, abstractmethod

import numpy as np


class BaseNumerical(ABC):
    """Base abstract class for generic numerical data.

    This class provides the infrastructure needed to implement mathematical operations in all the derived classes.

    The derived classes have to implement:
    - _apply_unary(self, function) that returns function(self)
    - _apply_binary(self, other, function) that returns function(self, other)
    - _apply_reduction(self, function) that returns function(self)

    """

    @abstractmethod
    def _apply_unary(self, function):
        raise NotImplementedError

    @abstractmethod
    def _apply_binary(self, other, function):
        raise NotImplementedError

    @abstractmethod
    def _apply_reduction(self, reduction):
        raise NotImplementedError

    def __add__(self, other):
        return self._apply_binary(other, np.add)

    __radd__ = __add__

    def __sub__(self, other):
        return self._apply_binary(other, np.subtract)

    def __rsub__(self, other):
        return -self._apply_binary(other, np.subtract)

    def __mul__(self, other):
        return self._apply_binary(other, np.multiply)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if other == 0:
            raise ValueError("Cannot divide by zero")
        return self._apply_binary(other, np.divide)

    def __rtruediv__(self, other):
        # This self._apply_binary(other, np.divide)
        # divives self by other, so, we reverse that
        # with ** -1
        return (self._apply_binary(other, np.divide)) ** -1

    def __pow__(self, other):
        return self._apply_binary(other, np.power)

    def __iadd__(self, other):
        return self + other

    def __isub__(self, other):
        return self - other

    def __imul__(self, other):
        return self * other

    def __itruediv__(self, other):
        return self / other

    def __ipow__(self, other):
        return self ** other

    def __neg__(self):
        return self._apply_unary(np.negative)

    def __abs__(self):
        return self._apply_unary(np.abs)

    def min(self):
        return self._apply_reduction(np.min)

    def max(self):
        return self._apply_reduction(np.max)

    def nanmin(self):
        return self._apply_reduction(np.nanmin)

    def nanmax(self):
        return self._apply_reduction(np.nanmax)

    def abs_min(self):
        """Return the minimum of the absolute value"""
        # skipcq PYL-W0212
        return abs(self)._apply_reduction(np.min)

    def abs_max(self):
        """Return the maximum of the absolute value"""
        # skipcq PYL-W0212
        return abs(self)._apply_reduction(np.max)

    def abs_nanmin(self):
        """Return the minimum of the absolute value ignoring NaNs"""
        # skipcq PYL-W0212
        return abs(self)._apply_reduction(np.nanmin)

    def abs_nanmax(self):
        """Return the maximum of the absolute value ignoring NaNs"""
        # skipcq PYL-W0212
        return abs(self)._apply_reduction(np.nanmax)

    def abs(self):
        return self._apply_unary(np.abs)

    amp = abs

    def real(self):
        return self._apply_unary(np.real)

    def imag(self):
        return self._apply_unary(np.imag)

    def sin(self):
        return self._apply_unary(np.sin)

    def cos(self):
        return self._apply_unary(np.cos)

    def tan(self):
        return self._apply_unary(np.tan)

    def arcsin(self):
        return self._apply_unary(np.arcsin)

    def arccos(self):
        return self._apply_unary(np.arccos)

    def arctan(self):
        return self._apply_unary(np.arctan)

    def sinh(self):
        return self._apply_unary(np.sinh)

    def cosh(self):
        return self._apply_unary(np.cosh)

    def tanh(self):
        return self._apply_unary(np.tanh)

    def arcsinh(self):
        return self._apply_unary(np.arcsinh)

    def arccosh(self):
        return self._apply_unary(np.arccosh)

    def arctanh(self):
        return self._apply_unary(np.arctanh)

    def sqrt(self):
        return self._apply_unary(np.sqrt)

    def exp(self):
        return self._apply_unary(np.exp)

    def log(self):
        return self._apply_unary(np.log)

    def log2(self):
        return self._apply_unary(np.log2)

    def log10(self):
        return self._apply_unary(np.log10)

    def conjugate(self):
        return self._apply_unary(np.conjugate)

    def phase(self):
        def unwrap_phase(*args, **kw):
            return np.unwrap(np.angle(*args, **kw))
        return self._apply_unary(unwrap_phase)


class BaseSeries(BaseNumerical):
    """
    Base class (not intended for direct use) for generic series data.
    """
    def __init__(self, x, y):
        """Initialise a generic series.

        :param array_like x: need to be strictly increasing.
        :param array_like y: values
        """
        assert len(x) >= 1, "initial_array must contain at least one sample."
        assert len(x) == len(y), "Times and Values length mismatch"
        assert np.diff(x).min() >= 0, "Timeseries is not monotonically increasing"

        self.x = np.array(x)
        self.y = np.array(y)

        self.is_complex = issubclass(self.y.dtype.type, complex)

    @property
    def xmin(self):
        """Return the minimum of the independent variable x.

        :rvalue: Minimum of x.
        :rtype: float
        """
        return self.x[0]

    @property
    def xmax(self):
        """Return the maximum of the independent variable x.

        :rvalue: Maximum of x
        :rtype: float
        """
        return self.x[-1]

    def __len__(self):
        """The number of data points."""
        return len(self.x)

    def __iter__(self):
        for x, y in zip(self.x, self.y):
            yield x, y

    def is_regularly_sampled(self):
        """Return whether the series is regularly sampled.

        If the series is only one point, an error is raised.

        :returns:  Is the series regularly sampled?
        :rtype:    bool
        """
        if len(self) == 1:
            raise RuntimeError(
                "Series is only one point, "
                "it does not make sense to compute dx"
            )

        dx = self.x[1:] - self.x[:-1]

        return np.allclose(dx, dx[0], atol=1e-14)

    def _apply_binary(self, other, function):
        """
        This is an abstract function that is used to implement mathematical operations with other series (if they have the same t).

        :param other: Other TimeSeries.
        :param callable function: Function to apply to the series.
        """
        # If the other object is of the same type
        if isinstance(other, type(self)):
            if (len(self.x) != len(other.x)) or (not np.allclose(other.x, self.x, atol=1e-14)):
                raise ValueError("The objects do not have the same t!")
            return self.__class__(self.x, function(self.y, other.y))

        # If it is a number
        if isinstance(other, (int, float, complex)):
            return self.__class__(self.t, function(self.y, other))

        # If we are here, it is because we cannot add the two objects
        raise TypeError("I don't know how to combine these objects")

    def __eq__(self, other):
        """Check for equality up to numerical precision."""
        if isinstance(other, type(self)):
            return np.allclose(self.x, other.x, atol=1e-14) and np.allclose(self.y, other.y, atol=1e-14)
        return False

    def _apply_unary(self, function):
        """Apply a unary function to the data.

        :param callable function: Function to apply to the series.
        """
        return self.__class__(self.x, function(self.y), self._CU)

    def _apply_reduction(self, reduction):
        """Apply a reduction to the data.

        :param callable function: Function to apply to the series.
        """
        return reduction(self.y)

    def shift(self, xshift):
        return self.__class__(self.x + xshift, self.y)  

    def crop(self, init=None, end=None):
        """
        Remove data outside the the interval ``[init, end]``. If ``init`` or ``end`` are not specified or None, it does not remove anything from this side.

        :param float init: Left boundary cut interval
        :param float end: Right boundary cut interval
        """
        x = self.x
        index = np.where((x >= init) & (x <= end))
        return self.__class__(x[index], self.y[index])
    
    # Define aliases
    clip = crop

    def nans_remove(self):
        """Filter out nans/infinite values.
        Return a new series with finite values only.

        :returns: A new series with only finite values.
        """
        msk = np.isfinite(self.y)
        return self.__class__(self.x[msk], self.y[msk])