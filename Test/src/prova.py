#!/usr/bin/env python

import numpy as np
from calcpol import *

# TUA = 24.
# TUE = 21.9
# TUD = 0.05
# TBM = 0.

# UK = 23.1
# DS = 88.3
# GK = 1.2

# ZZ = (TUE - TUA) / TUD
# print ZZ

# # if ZZ < 0:
# #     ZZ = ZZ - 1.0e-10
# # elif ZZ > 0:
# #     ZZ = ZZ + 1.0e-10

# IZ = ZZ
# IAB = abs(IZ) + 1

# print IAB

# for J in range(1, int(IAB) + 1):

#     ZZ = J - 1
#     UDEL = ZZ * TUD
#     if IZ < 0:
#         UDEL = -1 * UDEL

#     U = TUA + UDEL
#     UKU = UK**2
#     UU = U**2
#     UMU = UKU - UU
#     UU = 10.061*DS*UK/(UMU**2+UU*GK**2)
#     print UMU, UU
#     TBM = UMU*UU
#     print TBM

# Kk
# rr = np.arange(30, 15.99, -0.1)
# excfreq = 23.100
# dips = 40
# bw = 0.800

# cm-1
rr = np.arange(16000, 30001, 100)
excfreq = 23100
dips = 20
sigma = 2000.
# gamma = sigma * 2 * np.sqrt(np.log(2))
gamma = 1200.

uv_lor = uv_spec_lorentzian(rr, dips, excfreq, gamma)
np.savetxt('lor.txt', uv_lor)

uv_gau = uv_spec_gaussian(rr, dips, excfreq, sigma)
np.savetxt('gau.txt', uv_gau)

pol_im_part = pol_im(uv_lor)
np.savetxt('im.txt', pol_im_part, fmt='%.5e')
print pol_im_part

pol_real_part = pol_re(pol_im_part)
np.savetxt('re.txt', pol_real_part, fmt='%.5e')
print pol_real_part

# im = np.loadtxt('boh.txt')

# real = KK(rr, im)

# print im
# print real
