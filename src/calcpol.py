#!/usr/bin/env python

# DeVoe.py: calculation of UV and ECD spectra by DeVoe approach.
# Copyright (C) 2015  Daniele Padula dpadula85@yahoo.it
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import numpy as np
from scipy.constants import *

import util as u

# Calculation of some used constants
log10 = np.log(10)
pi_sq = np.pi**2

# Definition of some useful constants in CGS units
CGS_h = h * 1e7
CGS_c = c * 1e2

# Conversion factor for Dipolar Strength from D**2 to CGS (statC * cm)
squareD_to_CGS = 1e-36

# Conversion factor for Polarizability from cm**3 (CGS) to A**3 (human readable)
CGS_to_Acube = 1e24

# Prefactor for absorption spectra in CGS units 
CGS_CNST = 3000 * log10 * CGS_h * CGS_c / (32 * np.pi**3 * N_A)
CGS_CNST2 = 6909 / (8 * pi_sq * N_A)
CGS_CNST3 = 24 * pi_sq * N_A / (3298 * CGS_c**2)

#
# Lorentzian and Gaussian Bandwidth from Chirality, 2013, 25, 243
# Lorentzian gamma : HWHM
# Gaussian sigma : HW at 1/e height
#
def uv_spec_lorentzian(SpectralRange, DipS, ExcFreq, gamma):

    # Working in CGS units
    DipS = DipS * squareD_to_CGS
    cnst = 1 / (4 * CGS_CNST * np.pi**0.5)
    spec = cnst * DipS * ExcFreq * gamma / (gamma**2 + (SpectralRange - ExcFreq)**2)
    spec = np.c_[SpectralRange, spec]

    return spec


def uv_spec_gaussian(SpectralRange, DipS, ExcFreq, sigma):

    # Working in CGS units
    DipS = DipS * squareD_to_CGS
    cnst = 1 / (4 * CGS_CNST * np.pi**0.5)
    exponent = -1 * ((SpectralRange - ExcFreq) / sigma)**2
    spec = cnst * DipS * ExcFreq * np.exp(exponent) / sigma
    spec = np.c_[SpectralRange, spec]

    return spec


def pol_im(uvspec):

    # print SpectralRange
    SpectralRange = uvspec[:,0]
    epsilon = uvspec[:,1]

    # Working in CGS units, using SpectralRange in Hz
    # in fact speed of light at the numerator is missing
    pol_im_part =  CGS_CNST2 * epsilon / SpectralRange

    # Conversion from cm**3 to A**3
    # pol_im_part *= CGS_to_Acube  
    pol_im_part = np.c_[SpectralRange, pol_im_part]

    return pol_im_part


def pol_re(pol_im_part):

    SpectralRange = pol_im_part[:,0]
    pol_im = pol_im_part[:,1]
    pol_real_part = KK(SpectralRange, pol_im)
    pol_real_part = np.c_[SpectralRange, pol_real_part]

    return pol_real_part


def KK(e, im):

    re = np.zeros(len(im))
    prefactor = 2*(e[1]-e[0])/np.pi
    ee = e*e

    for i in range(len(im)):
        mask = e != e[i]
        re[i] = sum(e[mask]*im[mask]/(ee[mask]-e[i]*e[i]))

    return prefactor * re #+1 is needed only for epsilon1 calculation for xx, yy and zz projections


def KK_1(e, re):

    im = np.zeros(len(re))
    prefactor = 2*(e[1]-e[0])/np.pi
    ee = e*e

    for i in range(len(im)):
        mask = e != e[i]
        im[i] = sum(e[i]*(re[mask]-1)/(e[i]*e[i]-ee[mask]))

    return prefactor * im 


if __name__ == '__main__':
    pass
