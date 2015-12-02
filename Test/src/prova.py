#!/usr/bin/env python

import numpy as np
from calcpol import *

# cm-1
rr = np.arange(16000, 30001, 100)
excfreq = 23100
dips = 80
sigma = 2000.
gamma = sigma

# uv_lor = uv_spec_lorentzian(rr, dips, excfreq, gamma)
# np.savetxt('lor.txt', uv_lor)

# uv_gau = uv_spec_gaussian(rr, dips, excfreq, sigma)
# np.savetxt('gau.txt', uv_gau)

# pol_im_part = pol_im(uv_lor)
# np.savetxt('im.txt', pol_im_part, fmt='%.5e')
# print pol_im_part

# pol_real_part = pol_re(pol_im_part)
# np.savetxt('re.txt', pol_real_part, fmt='%.5e')
# print pol_real_part
