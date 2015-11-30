#!/usr/bin/env python

import numpy as np
from calcpol import *

rr = np.arange(24000, 21999, -50)
excfreq = 23100
dips = 88.3
bw = 1200

uv = uv_spec_lorentzian(rr, dips, excfreq, bw)
# np.savetxt('uv.txt', uv)
pol_im_part = pol_im(uv)
pol_real_part = pol_real(pol_im_part)
print len(pol_im_part)
print len(pol_real_part)
# np.savetxt('real2.txt', pol_real_part)
print pol_im_part
print pol_real_part
