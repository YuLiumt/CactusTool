#!/usr/bin/python

#----------------------------------------------------------
#
# This module computes differences in two given waveforms.
# 
#
# This module depends on DataAnalysis and DiscreteFunction
#
# Released under the MIT License.
# (C) Christian Reisswig 2010-2011
#
#----------------------------------------------------------


import sys
import os
import glob
from .DiscreteFunction import *
from .WaveFunction import *
from .harmonics import *
from .DataAnalysis import *



class WaveDiff:
    """A class holding two wavefunction objects 
       and their differences"""
       
    def __init__(self, WF1, WF2, dx, x0, x1):
	"""Initializes this class by loading two wavefunction.
	   The wavefunctions are assumed to be in complex form (_not_ decomposed into amplitude and phase!).
	   The waveforms will be resampled to new dx at common points and restricted to interval defined by [x0,x1]."""
	self.WF1 = WaveFunction(WF1.x, WF1.f)
	self.WF2 = WaveFunction(WF2.x, WF2.f)
	self.align = ""
	self.phi = 0
	self.omega = 0
	self.dx = dx
	self.x0 = x0
	self.x1 = x1
	self.isAligned = False
	# we cache the aligned waveforms
	self.aWF1 = WaveFunction([],[])
	self.aWF2 = WaveFunction([],[])
	
	
    def AlignInPhase(self, phi):
	"""Align given waveforms at same phase value phi"""
	self.align = "phase"
	self.phi = phi
	self.isAligned = False
	
    def AlignInFreq(self, omega):
	"""Align given waveforms at same frequency omega"""
	self.align = "freq"
	self.omega = omega
	self.isAligned = False
    
    def AlignAtAmpMax(self):
	"""Align two waveforms at their amplitude maximum"""
	self.align = "max"
	self.isAligned = False
    
    def NoAlign(self):
	self.align = ""
	self.isAligned = False
    
    def Align(self):
	"""Returns aligned waveforms according to settings"""
	
	if self.isAligned: 
	    # return copy of cached waveforms
	    return WaveFunction(self.aWF1.x, self.aWF1.f), WaveFunction(self.aWF2.x, self.aWF2.f)
	
	if (self.align == ""):
	    aWF1 = WaveFunction(self.WF1.x, self.WF1.f)
	    aWF2 = WaveFunction(self.WF2.x, self.WF2.f)
	if (self.align == "max"):
	    # align at amplitude maximum
	    tmax1 = self.WF1.Amplitude().FindAbsMaxInterpolated(1e-6)
	    tmax2 = self.WF2.Amplitude().FindAbsMaxInterpolated(1e-6)
	    aWF1 = WaveFunction(self.WF1.x-tmax1, self.WF1.f)
	    aWF2 = WaveFunction(self.WF2.x-tmax2, self.WF2.f)
	if (self.align == "freq"):
	    # get omega(t)
	    freq1 = self.WF1.Frequency()
	    freq2 = self.WF2.Frequency()
	    # invert omega(t) -> t(omega)
	    t1_vs_omega1 = DiscreteFunction(abs(freq1.f), freq1.x)
	    t2_vs_omega2 = DiscreteFunction(abs(freq2.f), freq2.x)
	    t1_vs_omega1 = t1_vs_omega1.MakeMonotonicCoordinates()
	    t2_vs_omega2 = t2_vs_omega2.MakeMonotonicCoordinates()
	    # find time t at which WF has frequency omega
	    t1 = t1_vs_omega1.Interpolate(self.omega)
	    t2 = t2_vs_omega2.Interpolate(self.omega)
	    #print t1, t2
	    # align waveforms at that frequency
	    aWF1 = WaveFunction(self.WF1.x-t1, self.WF1.f)
	    aWF2 = WaveFunction(self.WF2.x-t2, self.WF2.f)
	    
	    #self.aWF1 = WaveFunction(freq1.x-t1, freq1.f)
	    #self.aWF2 = WaveFunction(freq2.x-t2, freq2.f)
	    
	# interpolate waveforms to same points in same intervall
	self.aWF1, self.aWF2 = MakeConsistent(aWF1, aWF2, self.x0, self.x1, self.dx)
	
	self.isAligned = True
	
	return WaveFunction(self.aWF1.x, self.aWF1.f), WaveFunction(self.aWF2.x, self.aWF2.f)
    

    def AlignedPhases(self, reparam=""):
	"""Compute aligned/reparametrized phase as function of time, phase or frequency"""
	
	aWF1, aWF2 = self.Align()
	
	phase1 = aWF1.Phase()
	phase2 = aWF2.Phase()
	
	if (self.align == "max"):
	    # shift phase to zero at t=0 (where the maximum is)
	    shift1 = phase1.Interpolate(0.0)
	    shift2 = phase2.Interpolate(0.0)
	    
	    phase1.f -= shift1
	    phase2.f -= shift2
	
	return phase1, phase2


    def AlignedAmps(self, reparam=""):
	"""Computes aligned/reparametrized amplitude as function of time, phase or frequency"""
	
	aWF1, aWF2 = self.Align()
	
	amp1 = aWF1.Amplitude()
	amp2 = aWF2.Amplitude()
	
	if (reparam == "phase"):
	    # Reparametrize in terms of phase...
	    # ...and align, if necessary
	    phase1 = aWF1.Phase()
	    phase2 = aWF2.Phase()
	    if (self.align == "max"):
		# shift phase to zero at t=0 (where the maximum is)
		shift1 = phase1.Interpolate(0.0)
		shift2 = phase2.Interpolate(0.0)
	    
		phase1.f -= shift1
		phase2.f -= shift2
	    
	    # we take the absolute value of the phase since the phase can be negative and then our (positive) monotonicity does not hold
	    amp1_vs_phase1 = DiscreteFunction(abs(phase1.f), amp1.f)
	    amp1_vs_phase1 = amp1_vs_phase1.MakeMonotonicCoordinates()
	    
	    amp2_vs_phase2 = DiscreteFunction(abs(phase2.f), amp2.f)
	    amp2_vs_phase2 = amp2_vs_phase2.MakeMonotonicCoordinates()
	
	    amp1_vs_phase1, amp2_vs_phase2 = MakeConsistent(amp1_vs_phase1, amp2_vs_phase2, self.x0, self.x1, self.dx)
	
	    return amp1_vs_phase1, amp2_vs_phase2
	
	# if no reparametrization is requested, just return difference in time
	return amp1, amp2

    
    def DiffAmp(self, reparam=""):
	"""Computes difference in amplitude as function of time, phase or frequency"""
	
	a, b = self.AlignedAmps(reparam)
	
	return a-b

    def DiffPhase(self, reparam=""):
	"""Computes difference in phase as function of time, phase or frequency"""
	
	a, b = self.AlignedPhases(reparam)
	
	return a-b


    def Diff(self, reparam=""):
	"""Computes the complex wave difference as function of time, phase or frequency"""
	
	aWF1, aWF2 = self.Align()
	
	return diffWF


    def Mismatch(self, mass, ftype="h", Sn=SensitivityAdLIGO, frange=[10,8000]):
	"""Computes the mismatch for a given mass (in Msol), assuming the wave-function to be any of ftype=(h,news,psi4).
	   Sn is the detector sensitivity function, frange the (detector) frequency range in Hz."""
	
	aWF1, aWF2 = self.Align()
	
	MM = -1
	
	if ftype == "h":
	    aWF1SI = StrainToSI(aWF1, mass)
	    aWF2SI = StrainToSI(aWF2, mass)
	    MM = MisMatch(aWF1SI, aWF2SI, Sn, frange[0], frange[1])
	if ftype == "news":
	    aWF1SI = NewsToSI(aWF1, mass)
	    aWF2SI = NewsToSI(aWF2, mass)
	    MM = MisMatch(aWF1SI, aWF2SI, Sn, frange[0], frange[1], 'news')
	if ftype == "psi4":
	    aWF1SI = Psi4ToSI(aWF1, mass)
	    aWF2SI = Psi4ToSI(aWF2, mass)
	    MM = MisMatch(aWF1SI, aWF2SI, Sn, frange[0], frange[1], 'psi4')
	
	return MM
	
	
	