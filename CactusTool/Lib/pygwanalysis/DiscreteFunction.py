#----------------------------------------------------------
#
# Classes and routines for representing and manipulating 
# discrete functions.
#
# Released under the MIT License.
# (C) Christian Reisswig 2009-2011
#
#----------------------------------------------------------


import sys

from math import *
from numpy import *
from .Utils import *

debug = 0

class DiscreteFunction:
    """Class describing discrete functions which transform 1D data"""

    def HasUniformSpacing(self):
        """Checks if this function has a unifrom grid spacing"""
        if (len(self.x) > 1):
            # check if we have uniform spacing
            dx = self.x[1]-self.x[0]
            for ii in range(1, len(self.x)-1):
                if abs(dx-abs(self.x[ii]-self.x[ii+1])) > 1e-7:
                    return False
            return True
        return False
            

    def __init__(self, x_, f_, g_=None):
        """Initialise a 1D transformation function"""
        self.x = array(x_)
        if g_ != None:
            self.f = array(f_) + 1j*array(g_)
        else:
            self.f = array(f_)
        if self.HasUniformSpacing() == True: 
            self.dx = self.x[1]-self.x[0]
        else:
            self.dx = 0.0

    def __del__(self):
        self.x = []
        self.f = []
    
    def Load(self, filename, col1=1, col2=2, col3=3):
        """load a 1d function from file"""
        
        fp = open(filename, 'r')
        x_ = []
        f_ = []
        for line in fp:
            if line[0]=='#' or line[0]=='"' or len(line.strip())==0:
                continue
            data = []
            #data = map(float, line.split())
            data = line.split()
            x_.append(float(data[col1-1]))
            if (len(data) == 2 or col3<=0):
                f_.append(float(data[col2-1]))
            else:
                f_.append(complex(float(data[col2-1]), float(data[col3-1])))
        fp.close()
        
        self.x = array(x_)
        self.f = array(f_)
        if self.HasUniformSpacing() == True: 
            self.dx = self.x[1]-self.x[0]
        else:
            self.dx = 0.0
    
    def MakeMonotonicCoordinates(self, maxvariance=1e10):
        """Given a function whose coordinates are not monotonic, we remove all entries that violate monotonicity.
           maxvariance is the maximal change that may occur from one point to the next"""
        x_ = []
        f_ = []
                
        jj = 0
        for ii in range(0, len(self.x)-1):
            if (self.x[jj] >= self.x[ii]):
                continue
            if (self.x[ii] > self.x[jj]+maxvariance):
                continue
            jj = ii
            x_.append(self.x[ii])
            f_.append(self.f[ii])
        
        return self.__class__(x_, f_)
    
    
    def Restrict(self, x0, x1):
        """Restrict 1d function to intervall between [x0, x1]"""
                
        x_ = []
        f_ = []
        
        for ii in range(0, len(self.x)-1):
            if (self.x[ii] < x0 or self.x[ii] > x1):
                continue
            x_.append(self.x[ii])
            f_.append(self.f[ii])
        
        return self.__class__(x_, f_)
    
    def Invert(self, x0, x1):
        """Invert function on a specific interval"""
        DF = self.Restrict(x0, x1)
        DF = self.__class__(DF.f, DF.x)
        DF = DF.MakeMonotonicCoordinates()
        return DF
        
    
    def Extend(self, n, dx, g, after=1):
        """Extend the function by n function values g separated with spacing dx"""
        x_ = []
        f_ = []
        
        # add new values at beginning
        if after==0:
            for ii in range(0, n):
                x_.append(self.x[0]+(ii-n)*dx)
                f_.append(g(x_[-1]))
        # copy old values
        for ii in range(0, len(self.x)):
            x_.append(self.x[ii])
            f_.append(self.f[ii])
        # add new values at the end
        if after==1:
            for ii in range(0, n):
                x_.append(self.x[-1]+(ii+1)*dx)
                f_.append(g(x_[-1]))
        
        return self.__class__(x_, f_)
    
    def Append(self, DF):
        """Appends another discrete function to this function"""
        x_ = []
        f_ = []
        
        # copy old values
        for ii in range(0, len(self.x)):
            x_.append(self.x[ii])
            f_.append(self.f[ii])
        # add new values at the end
        for ii in range(0, DF.Length()):
            x_.append(self.x[-1]+self.dx+DF.x[ii]-DF.x[0])
            f_.append(DF.f[ii])
        
        return self.__class__(x_, f_)
    
    def FindNearestNeighbor(self, x, BruteForceSearch=False):
        """Find nearest neighboring point (index) to x"""
        n = 0
        if (self.dx > 0 and BruteForceSearch == False):
            # We have uniform data and can compute the nearest point
            n = int(round((x-self.x[0])/self.dx))
        if (self.dx <= 0 and BruteForceSearch == False):
            # Use bijection to find our point
            if x <= self.x[0]: return 0
            if x >= self.x[-1]: return self.Length()-1
            found = False
            a = 0
            b = self.Length()-1
            i = int(round((b-a)*0.5))
            while (found == False):
                if (b-a == 1):
                    #print("(x, a, b, x[a], x[b])=", x, a, b, self.x[a], self.x[b])
                    if (abs(x-self.x[b]) < abs(x-self.x[a])): return b
                    else: return a
                if (x == self.x[i]): return i
                if (x < self.x[i]): 
                    b = i
                else:
                    a = i
                i = a+int(round((b-a)*0.5))
            #print("Error: FindNearestNeighbor: bijection search not implemented")
        if (BruteForceSearch == True):
            for ii in range(0, len(self.x)-1):
                if (self.x[ii] <= x and x <= self.x[ii+1]):
                    if (abs(x-self.x[ii]) < abs(x-self.x[ii+1])):
                        n = ii
                    else:
                        n = ii+1
                    break
        return n
    
    def Interpolate(self, xi):
        """interpolate to point x"""
        
        interp_order = 5
        
        # Don't do anything if we are outside of domain, just copy value from boundary points
        if (xi < self.x[0]): return self.f[0]
        if (xi > self.x[-1]): return self.f[-1]
        
        # find closest point to interpolation point
        n = self.FindNearestNeighbor(xi)
        
        #if (abs(xi-self.x[n]) < 1e-8): print("Skipping interpolation: point too close to gridpoint."; return self.f[n])
        
        npoints = interp_order+1
        if (self.x[n] < xi and npoints % 2 == 0):
            nstart = n - int(npoints/2) + 1
        else:
            nstart = n - int(npoints/2)
        if (nstart < 0): nstart = 0
        if (nstart+npoints >= len(self.x)): nstart = len(self.x)-npoints
        
        res = 0.0
        for ii in range(nstart, nstart+npoints):
            l = 1.0
            for jj in range(nstart, nstart+npoints):
                if (ii != jj):
                    l = l * (xi-self.x[jj]) / (self.x[ii]-self.x[jj])
            res = res + l*self.f[ii]
        
        return res

    def QuadraticSpline(self, x, z0=0):
        """Use a quadratic spline to interpolate at point x"""
        
        # Don't do anything if we are outside of domain, just copy value from boundary points
        if (xi < self.x[0]): return self.f[0]
        if (xi > self.x[-1]): return self.f[-1]
        
        # find closest point to interpolation point
        n = 0
        for ii in range(0, len(self.x)-1):
            if (self.x[ii] <= xi and xi <= self.x[ii+1]):
                n = ii
                break
        
        if (abs(xi-self.x[n]) < 1e-8): return self.f[n]
        
        nstart = n
        
        #self.f[ii] + 
        
        
        return res


    def Resample(self, dx, intervall=0):
        """Resample field to resolution dx"""
        
        if intervall == 0:
            xstart = self.x[0]
            npoints = int((self.x[-1]-self.x[0]) / dx)+1
        else:
            xstart = intervall[0]
            npoints = int((intervall[1]-intervall[0]) / dx)+1
        
        res = self.__class__(zeros(npoints, float), zeros(npoints, type(self.f[0])))
        res.x[0] = xstart
        res.f[0] = self.Interpolate(xstart)
        res.dx = dx
        
        for ii in range(1, npoints):
            res.x[ii] = xstart+dx*ii
            res.f[ii] = self.Interpolate(res.x[ii])
        
        return res


    def TransformCoord(self, g):
        """ Transform coordinate according to function g """
        self.x = map(g, self.x)
        return self
            
    def Length(self):
        """Returns number of datapoints"""
        return len(self.x)


    def Union(self, f2):
        """ Returns the union of points between this and function f"""
        x_ = []
        jj = 0
        npoints = 0

        # find union number of points
        for ii in range(0, self.Length()):
            while self.x[ii] > f2.x[jj] and jj != len(f2.x)-1:
                x_.append(f2.x[jj])
                npoints += 1
                if (jj < len(f2.x)-1):
                    jj += 1
            if (self.x[ii] == f2.x[jj]):
                if (jj < len(f2.x)-1):
                    jj += 1
            x_.append(self.x[ii])
            npoints += 1
        
        res = self.__class__(x_, zeros(npoints, type(self.f)))
        
        return res
        

    def NormL1(self):
        """Calculates L1-norm of function"""
        norm = 0
        for ii in range(0, self.Length()):
            norm += abs(self.f[ii])
        return norm/self.Length()

    def NormL2(self):
        """Calculates L2-norm of function"""
        norm = 0
        for ii in range(0, self.Length()):
            norm += abs(self.f[ii])**2
        return sqrt(norm)/self.Length()
        
        return norm


    def NormSup(self):
        """Calculates infinity-norm of function"""
        
        norm = 0
        for ii in range(0, self.Length()):
            if (abs(self.f[ii]) > norm): norm = abs(self.f[ii])
        
        return norm


    def FindMax(self):
        """Finds grid index location of discrete maximum"""
        
        fmax = -1e50
        imax = 0
        for ii in range(0, self.Length()):
            if (self.f[ii] > fmax):
                imax = ii
                fmax = self.f[ii]
        
        return imax

    def FindAbsMax(self):
        """Finds grid index location of absolute discrete maximum"""
        
        fmax = -1e50
        imax = 0
        for ii in range(0, self.Length()):
            if (abs(self.f[ii]) > fmax):
                imax = ii
                fmax = abs(self.f[ii])
        
        return imax

    def FindAbsMaxInterpolated(self, precision=1e-4):
        """Finds location of continuous maximum via interpolation"""
        
        ii = self.FindAbsMax()  # find discrete maximum
        #print("Discrete Maximum at", self.x[ii])
        
        DF2 = self.Restrict(self.x[ii-1], self.x[ii+1])
        
        # Preselect intervall based on first derivative of 2nd order Lagrange polynomial (its derivative has only one root!)
        a = self.x[ii-1]
        b = self.x[ii+1]
        DF2 = self.Restrict(a, b)
        if Lagrange_derivative(self.x[ii], DF2.x, DF2.f, 2)*Lagrange_derivative(self.x[ii+1], DF2.x, DF2.f, 2) < 0:
            a = self.x[ii]
        else:
            b = self.x[ii]

        #xroot = find_root([a, b], lambda x: Lagrange_derivative(x, DF2.x, DF2.f, order), 0, precision)

        # Find maximum by interpolating between a and b with given precision
        prev = self.Interpolate(a)
        a += precision
        while (a != b):
            cur = self.Interpolate(a)
            if cur < prev:
                break
            prev = cur
            a += precision
        xroot = a

        #print("Maximum at", xroot)
        return xroot    


    def DifferentiateFunction(self):
        """ Differentiate this function and return a function"""
        
        npoints = self.Length()
        res = self.__class__(self.x, self.f*0.0)
        
        # TODO: fail if not uniform spacing
        self.dx = self.x[1]-self.x[0]
        res.dx = self.dx
        
        # Simple second order centred differentiation
        
        for ii in range (1, self.Length()-1):
            res.f[ii] = (self.f[ii+1] - self.f[ii-1])/(2*res.dx)
        
        return res
    
    
    def IntegrateFunction(self):
        """ Integrate this function and return a function"""
        
        npoints = self.Length()
        res = self.__class__(self.x, self.f*0.0)
        
        # TODO: fail if not uniform spacing
        self.dx = self.x[1]-self.x[0]
        res.dx = self.dx
        
        # 4th-order Adams-Moulton coeffs
        beta = array([9.0/24.0, 19.0/24.0, -5.0/24.0, 1.0/24.0])
        
        for ii in range (3, self.Length()):
            s = 0
            for k in range(0, 4):
                s += beta[k]*self.f[ii-k]
            s *= self.dx
            res.f[ii] = s+res.f[ii-1]
        
        return res
    
    
    def Integrate(self, a, b, Resample=True, secondOrder=False):
        """ Integrate this function in [a,b] and return a number. If a,b do not line up with gridpoints then resample"""
        
        res = 0
        self.dx = self.x[1]-self.x[0]
        
        if (a<self.x[0]): a = self.x[0]
        if (b>self.x[-1]): b = self.x[-1]
        
        nstart = int((a-self.x[0])/self.dx)
        nend = int((b-self.x[0])/self.dx)
        if (Resample == True):
            ntot = nend-nstart
            dx = (b-a) / ntot
            dummy = self.__class__(zeros(ntot, float), zeros(ntot, type(self.f[0])))
            for ii in range(0, dummy.Length()):
                dummy.x[ii] = a+ii*dx
                dummy.f[ii] = self.Interpolate(dummy.x[ii])
            nstart = 0
            nend = dummy.Length()-1
        else:
            dummy = self.__class__(self.x, self.f)
            for ii in range(nstart, nend+1):
                dummy.x[ii] = self.x[ii]
                dummy.f[ii] = self.f[ii]
        
        if (secondOrder):
            for i in range(nstart, nend+1):
                res = res + coeff2nd(i-nstart, nend-nstart)*dummy.f[i]
        else:
            for i in range(nstart, nend+1):
                res = res + Simpson_coeff(i-nstart, nend-nstart)*dummy.f[i]

        res = self.dx * res

        return res



    def FirstDerivative(self, order=4, dx=-1):
        """ Take first derivative of this function"""
        
        if (self.HasUniformSpacing() == False and dx==-1): print("FirstDerivative: f has no uniform spacing!")
        
        npoints = len(self.x)   
        res = self.__class__(zeros(npoints, float), zeros(npoints, type(self.f[0])))
        self.dx = self.x[1]-self.x[0]
        
        # Use dx given by parameter
        if (dx > 0):
            for ii in range(0, res.Length()):
                res.x[ii] = self.x[ii]
                res.f[ii] = (self.Interpolate(res.x[ii]+dx) - self.Interpolate(res.x[ii]-dx)) / (2.0*dx)
        else: # use internal dx
            if (order==2):
                res.x[0] = self.x[0]
                for ii in range(1, res.Length()-1):
                    res.x[ii] = self.x[ii]
                    res.f[ii] = (self.f[ii+1] - self.f[ii-1]) / (2.0*self.dx)
                res.x[-1] = self.x[-1]
            if (order==4):
                res.x[0] = self.x[0]
                res.x[1] = self.x[1]
                for ii in range(2, res.Length()-2):
                    res.x[ii] = self.x[ii]
                    res.f[ii] = (-self.f[ii+2]/12. + 2.*self.f[ii+1]/3. - 2.*self.f[ii-1]/3. + self.f[ii-2]/12.) / (self.dx)
                res.x[-1] = self.x[-1]
                res.x[-2] = self.x[-2]
        return res


    def SecondDerivative(self, order=2, dx=-1):
        """ Take second derivative of this function"""
    
        npoints = len(self.x)   
        res = self.__class__(zeros(npoints, float), zeros(npoints, type(self.f[0])))
        self.dx = self.x[1]-self.x[0]
        
        # Use dx given by parameter
        if (dx > 0):
            for ii in range(0, res.Length()):
                if order==2:
                    res.f[ii] = (self.Interpolate(res.x[ii]+dx) - 2*self.f[ii] + self.Interpolate(res.x[ii]-dx)) / (dx**2)
        else: # use internal dx
            if order==2:
                res.x[0] = self.x[0]
                for ii in range(1, res.Length()-1):
                    res.x[ii] = self.x[ii]
                    res.f[ii] = (self.f[ii+1] - 2*self.f[ii] + self.f[ii-1]) / (self.dx**2)
                res.x[-1] = self.x[-1]
            if order==4:
                res.x[0] = self.x[0]
                res.x[1] = self.x[1]
                for ii in range(2, res.Length()-2):
                    res.x[ii] = self.x[ii]
                    res.f[ii] = (-self.f[ii+2]/12 + 4*self.f[ii+1]/3 - 5*self.f[ii]/2 + 4*self.f[ii-1]/3 - self.f[ii-2]/12) / (self.dx**2)
                res.x[-1] = self.x[-1]
                res.x[-2] = self.x[-2]
        return res


    def FourierTransform(self, cutNegFreq=False):
        """ Fourier-transform of this function"""

        npoints = len(self.x)
        res = self.__class__(zeros(npoints, float), zeros(npoints, type(self.f[0])))

        dt = self.dx
        res.f = dt * fft.fft(self.f) #/2  # divide by two because we're going to cut the negative frequencies!
        
        # get the freq
        
        #for i in range(0, self.Length()):
        #    res.x[i] = 2*pi * i*1.0/(self.Length()*dt)
        res.x = fft.fftfreq(self.Length(), dt)
        res.dx = res.x[1]-res.x[0]
        
        if cutNegFreq==False:
            # Undo standard FFT ordering where negative frequencies are stored beyond i=n/2 (n=length of array)
            # so that we have a continuous parametrization
            # This needs to be undone later!
            res.x = fft.fftshift(res.x)
            res.f = fft.fftshift(res.f)
        else:
            # rather delete negative frequency components
            res.x = res.x[0:npoints/2]
            res.f = res.f[0:npoints/2]
        
        return res


    def InverseFourierTransform(self, shiftNegFreq=False, noNegFreq=False):
        """ Fourier-transform of this function"""

        npoints = len(self.x)
        res = self.__class__(zeros(npoints, float), zeros(npoints, type(self.f[0])))
        
        if noNegFreq==False:
            # assuming we have a continuous parameterization, we have to put
            # this to FFT standard form where negative frequencies are stored beyond n/2 (n=length of array)
            self.x = fft.ifftshift(self.x)
            self.f = fft.ifftshift(self.f)

        df = self.dx
        res.f = fft.ifft(self.f)
        
        # get the time
        res.x = fft.fftfreq(self.Length(), df)
        res.dx = res.x[1]-res.x[0]
        if noNegFreq==False:
            res.x = fft.fftshift(res.x)
        if shiftNegFreq==True:    
            res.f = fft.fftshift(res.f)
        res.f = res.f / (res.dx)
        
        return res



    def Smooth(self):
        """Smooth this function """
    
        res = self.__class__(zeros(npoints, float), zeros(npoints, type(self.f[0])))
        
        return res


    def __add__(self, rval):
        res = self.Union(rval) 
        for ii in range(0, res.Length()):
            res.f[ii] = self.Interpolate(res.x[ii]) + rval.Interpolate(res.x[ii])
        return res

    def __sub__(self, rval):
        res = self.Union(rval) 
        for ii in range(0, res.Length()):
            res.f[ii] = self.Interpolate(res.x[ii]) - rval.Interpolate(res.x[ii])
        return res

    # implement pointwise multiplication (to apply windowing to a wavefunction)
    def __mul__(self, rval):
        res = self.Union(rval)
        for ii in range(0, res.Length()):
            res.f[ii] = self.Interpolate(res.x[ii]) * rval.Interpolate(res.x[ii])
        return res

    def __div__(self, rval):
        res = self.Union(rval)
        for ii in range(0, res.Length()):
            res.f[ii] = self.Interpolate(res.x[ii]) / rval.Interpolate(res.x[ii])
        return res

    def __abs__(self):
        npoints = len(self.x)
        res = self.__class__(self.x, zeros(npoints, float))
        for ii in range(0, Length()):
            res.f[ii] = abs(self.f[ii])
        return res

    def __pos__(self):
        npoints = len(self.x)
        res = self.__class__(self.x, self.f)
        return res

    def __neg__(self):
        npoints = len(self.x)
        res = self.__class__(self.x, -self.f)
        return res

    def real(self):
        npoints = len(self.x)
        res = self.__class__(self.x, zeros(npoints, float))
        for ii in range(0, self.Length()):
            res.f[ii] = self.f[ii].real
        return res

    def imag(self):
        npoints = len(self.x)
        res = self.__class__(self.x, zeros(npoints, float))
        for ii in range(0, self.Length()):
            res.f[ii] = self.f[ii].imag
        return res

    def conjugate(self):
        npoints = len(self.x)
        res = self.__class__(self.x, zeros(npoints, type(self.f[0])))
        for ii in range(0, self.Length()):
            res.f[ii] = self.f[ii].conjugate()
        return res


    # implement iterator
    def __iter__(self):
        self.index = 0
        return self
    
    def next(self):
        if self.index == self.Length()-1:
            raise StopIteration
        self.index = self.index + 1
        return [self.x[self.index], self.f[self.index]]
    
    
    def __repr__(self):
        s = "# x    f(x) \n"
        #if type(self.f[0])==float:
        for ii in range(0, self.Length()):
            s += str(self.x[ii])+" "+str(self.f[ii])+"\n"
        return s
        #if type(self.f[0])==complex:
        #    for ii in range(0, self.Length()):
        #       s += str(self.x[ii])+" "+str(self.f[ii].real)+" "+str(self.f[ii].imag)+"\n"
        #    return s
    
        



#---------------------------------------------------------------
#
# Here come functions that operate on DiscreteFunction objects
#
#---------------------------------------------------------------



def ScalProd(f1, f2, kernel):
    """ Calculates scalar-product of functions f1, f2 over kernel """
    
    res = 0.0
    
    
    return res


def MakeConsistent(DF1, DF2, x0, x1, dx):
    """ Interpolates functions to same interval on same points """
    
    if (x0 < max(DF1.x[0], DF2.x[0])): x0 = max(DF1.x[0], DF2.x[0])
    if (x1 > min(DF1.x[-1], DF2.x[-1])): x1 = min(DF1.x[-1], DF2.x[-1])

    if ((x1-x0) < 0):
        print("Warning: MakeConsistent: Empty overlap.")
        npoints = 1
    else:
        npoints = int ((x1-x0)/dx)+1
    
    newDF1 = DF2.__class__(zeros(npoints, float), zeros(npoints, type(DF1.f[0])))
    newDF2 = DF2.__class__(zeros(npoints, float), zeros(npoints, type(DF2.f[0])))
    
    for ii in range(0, npoints):
        newDF1.x[ii] = x0+ii*dx
        newDF2.x[ii] = newDF1.x[ii]
        newDF1.f[ii] = DF1.Interpolate(newDF1.x[ii])
        newDF2.f[ii] = DF2.Interpolate(newDF2.x[ii])

    return newDF1, newDF2


def Blend(DF1, DF2, x0, x1, a=lambda x: x):
    """Blend two functions DF1->Df2 over an intervall [x0,x1] using function a: [0,1]->[0,1].
       The blending is applied according to the convex combination (1-a)*DF1 + a*DF2
       where a takes values in [0,1]. Usually, this is just the identity, i.e. a(x)=x."""

    if (x1 <= x0): print("Blend error: x1<=x0 !")

    x = []
    f = []

    # normalize
    l = x1-x0

    for ii in range(0, DF1.Length()):
        if (DF1.x[ii] < x0):
            x.append(DF1.x[ii])
            f.append(DF1.f[ii])
        if (DF1.x[ii] >= x0 and DF1.x[ii] <= x1):
            # project current coordinate x in intervall c\in[0,1]
            c = (DF1.x[ii]-x0)/(x1-x0)
            x.append(DF1.x[ii])
            # blend
            f.append((1-a(c))*DF1.f[ii] + a(c)*DF2.Interpolate(DF1.x[ii]))
    
    # append remaining part of DF2
    for ii in range(0, DF2.Length()):
        if (DF2.x[ii] > x1):
            x.append(DF2.x[ii])
            f.append(DF2.f[ii])

    return DF2.__class__(x, f)



def AlignOverRange(DF1, DF2, x0, x1, start_xoff, delta, dx, accuracy=1e-2, monotonic=False):
    """Aligns two functions by minimizing their difference in the L2 norm and 
       returns the x-shift that needs to be applied.
       [x0,x1]: The minimization interval.
       start_xoff: An initial guess to speed up computation.
       delta: the shift will be applied in the range [start_xoff-delta, start_xoff+delta]. 
              The smaller delta, the faster (but keep in mind that you need to be sure that the minimum is located in that range!
       dx: The delta spacing of interpolated function values on which the pointwise difference is taken.
       accuracy: The accuracy in finding the shift that minimizes the L2 norm of the error ('the delta spacing of the shift').
       monotonic: Set to 'True' if function is monotonic (this will greatly speed up computation)"""

    # Keep DF1 fixed. Find offset x_off that needs to be applied to DF2.x
    # in order to make the L2-norm of the difference over a range [x0,x1]
    # as small as possible.

    x_off = start_xoff-delta
    minErr = 10e10
    minX = x_off
    trueMinX = x_off
    firstTime = True

    # If function is monotonic, we can use a faster algorithm
    if (monotonic==True):
        factor = 2.0*delta/accuracy*0.1
        maxX = start_xoff+delta
        while (factor >= 1.0 or firstTime == True):
            x_off = minX
            minErr = 10e10
            while (x_off < maxX):
                DF2shift = DiscreteFunction(DF2.x+x_off, DF2.f)
                subDF1, subDF2 = MakeConsistent(DF1, DF2shift, x0, x1, dx)
                err = 0
                for ii in range (0, subDF1.Length()):
                    err += (subDF1.f[ii]-subDF2.f[ii])**2
                err = sqrt(err)
                #print(minErr, err)
                if (err < minErr):
                    minErr = err
                    minX = x_off-accuracy*factor
                    trueMinX = x_off
                else:
                    maxX = x_off
                    break
                x_off += accuracy*factor
            factor *= 0.1
            firstTime = False
        return trueMinX
        
        print("ERROR! ...using brute force method!")
    
    # Proceed with brute force algorithm...

    x_off = start_xoff-delta

    minErr = 10e10
    minX = x_off


    # We shift DF2 and compute the error over range [x0,x1]
    while (x_off < start_xoff+delta):
        DF2shift = DiscreteFunction(DF2.x+x_off, DF2.f)
        subDF1, subDF2 = MakeConsistent(DF1, DF2shift, x0, x1, dx)
        err = 0
        for ii in range (0, subDF1.Length()):
            err += (subDF1.f[ii]-subDF2.f[ii])**2
        err = sqrt(err)
        #print(err)
        if (err < minErr):
            minErr = err
            minX = x_off
        x_off += accuracy

    return minX


