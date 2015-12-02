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
import argparse as arg

import parseinput as pi
import calcpol as cp
import makecoords as m
import util as u

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''
    DeVoe calculation of absorption and ECD spectra of multichromophoric
    systems.''')

    # Optional arguments
    parser.add_argument('-i', '--infile', default='devoe.in', help='''
    Input file containing the definitions of dipole types, application
    points and directions.''')

    parser.add_argument('-s', '--structure', default='structure.pdb', help='''
    File containing the structure of the system.''')

    parser.add_argument('-sf', '--sf', default=None, choices=['mol2', 'pdb', 'xyz'],
    help='''Format of the structure file.''')

    parser.add_argument('-ls', '--lineshape', default='lor', choices=['lor', 'gau'],
    help='''Spectral lineshape.''')

    parser.add_argument('--min', default=16000, help='''Low energy limit in wavenumbers.''')

    parser.add_argument('--max', default=35000, help='''High energy limit in wavenumbers.''')

    parser.add_argument('-v', '--verbosity', action='store_true', default=False,
    help='''Verbosity of the output.''')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    title = '''    
###########         ____     __     __                               ########### 
###########        |  _ \  __\ \   / /__   ___   _ __  _   _         ###########
###########        | | | |/ _ \ \ / / _ \ / _ \ | '_ \| | | |        ###########
###########        | |_| |  __/\ V / (_) |  __/_| |_) | |_| |        ###########
###########        |____/ \___| \_/ \___/ \___(_) .__/ \__, |        ###########
###########                                     |_|    |___/         ###########
'''

    #
    # Get script arguments
    #
    args = options()

    infile = args.infile
    structure = args.structure
    sf = args.sf
    lineshape = args.lineshape
    u.verbosity = args.verbosity

    SpecRange = np.arange(args.min, args.max + 1, (args.max - args.min)/500.)

    if u.verbosity:
        print
        print u.banner(ch='#', length=80)
        print title
        print u.banner(ch='#', length=80)
        print

    #
    # Recapitulation of input files
    #
    if u.verbosity:
        print(" Input files:")
        print(" Dipole orientation file: %s" % infile)
        print(" Structure file: %s" % structure)
        print

    #
    # Try to guess structure's format from extension if not specified
    #
    if not sf:
        sf = args.structure.split('.')[-1]

    #
    # Get structure
    #
    if sf == 'mol2':
        coords = pi.parse_MOL2(structure)

    elif sf == 'pdb':
        coords = pi.parse_PDB(structure)

    elif sf == 'xyz':
        coords = pi.parse_XYZ(structure)

    else:
        print u.banner(text='ERROR', ch='#', length=80)
        print("%s file format not recognized." % structure)
        print
        sys.exit()

    #
    # Get dipole types, application points and orientations in terms of atomic
    # indexes from the input file
    #
    coords = np.array(coords)
    dipole_types, dipoles_info = pi.parse_input(infile)

    #
    # Calculate  imaginary and real parts of the desired type of polarizability
    # for each dipole type
    #

    if lineshape == 'lor':
        ls_funct = cp.uv_spec_lorentzian

    elif lineshape == 'gau':
        ls_funct = cp.uv_spec_gaussian

    for dip_type, params in dipole_types.iteritems():

        info = params[0]
        pol_type = params[1]

        DipStrength = info[0]
        ExcFreq = info[1]
        damping = info[2]

        uvspec = ls_funct(SpecRange, DipStrength, ExcFreq, damping)
        pol_im = cp.pol_im(uvspec)
        pol_re = cp.pol_re(pol_im)

        # np.savetxt('uv.txt', uvspec, fmt='%.6e')
        # np.savetxt('im.txt', pol_im, fmt='%.6e')
        # np.savetxt('re.txt', pol_re, fmt='%.6e')

    #
    # Transform the atomic index representation in coordinates
    #
    for cent, weight, dipoles in dipoles_info:
        
        appl_point = m.make_cent(cent, weight, coords)

        for dipole in dipoles:

            dipole = m.make_dipo(dipole, dipole_types, coords)

            # Translate the dipole in its application point
            # dipole += appl_point

        if u.verbosity:
            print

    if u.verbosity:
        print u.banner(ch='#', length=80)
        print
