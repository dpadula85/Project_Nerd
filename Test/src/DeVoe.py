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

    parser.add_argument('-v', '--verbosity', default=0, action='count', help='''Verbosity of the output.''')

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

    if u.verbosity >= 1:
        print
        print(u.banner(ch='#', length=80))
        print(title)
        print(u.banner(ch='#', length=80))
        print

    #
    # Recapitulation of input files
    #
    if u.verbosity >= 1:
        print(" > READING INPUT FILES...")
        print

    if u.verbosity >= 2:
        print("   Dipole orientation file : %s" % infile)
        print("   Structure file          : %s" % structure)
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
        print(u.banner(text='ERROR', ch='#', length=80))
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
    # Transform the atomic index representation in coordinates and store the
    # application point coordinates and the direction vector in arrays
    #
    centers = np.array([]).reshape(0,3)
    orientations = np.array([]).reshape(0,3)
    # pol_types_list = []

    if u.verbosity >= 1:
        print(" > ORIENTING DIPOLES...")
        print

    for cent, weight, dipoles in dipoles_info:

        # Get application point coordinates
        appl_point = m.make_cent(cent, weight, coords)

        for dipole in dipoles:

            # Get the unit vector for the orientation
            e = m.make_dipo(dipole, dipole_types, coords)
            
            # Translate the dipole in its application point
            # e += appl_point

            centers = np.vstack((centers, appl_point))
            orientations = np.vstack((orientations, e))
            # pol_types_list.append(dipole[0])

        if u.verbosity >= 2:
            print

    #
    # Build Interaction Matrix
    #
    if u.verbosity >= 1:
        print(" > BUILDING INTERACTION MATRIX...")
        print

    G = np.eye(len(centers))

    # We only need to build the G matrix elements with i != j
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):

            # Comparison element by element of the coordinates of the center
            if (centers[i] == centers[j]).all():

                # set the interaction matrix element to 0
                G[i,j] = 0.0
                G[j,i] = 0.0

            else:

                # Calculate the interaction
                e_i = orientations[i]
                e_j = orientations[j]
                r_ij = centers[j] - centers[i]
                e_ij = r_ij / np.linalg.norm(r_ij)
                G_ij = (np.dot(e_i, e_j) - 3 * np.dot(e_i, e_ij) * np.dot(e_i, e_ij)) / np.linalg.norm(r_ij)**3
                G[i,j] = G_ij
                G[j,i] = G[i,j]

    #
    # Calculate imaginary and real parts of the desired type of polarizability
    # for each dipole type
    #
    if lineshape == 'lor':
        ls_funct = cp.uv_spec_lorentzian

    elif lineshape == 'gau':
        ls_funct = cp.uv_spec_gaussian

    if u.verbosity >= 1:
        print(" > CALCULATING POLARIZABILITIES...")
        print

    if u.verbosity >= 2:
        print("   Lineshape: %s" % lineshape)
        print

    # pol_types_dict = {}
    for dip_type, params in dipole_types.iteritems():

        info = params[0]
        pol_type = params[1]

        if pol_type == 'mag':
            bj = params[2]

        else:
            bj = 0.0

        DipStrength = info[0]
        ExcFreq = info[1]
        damping = info[2]

        uvspec = ls_funct(SpecRange, DipStrength, ExcFreq, damping)
        pol_im = cp.pol_im(uvspec)
        pol_re = cp.pol_re(pol_im)

        # pol_types_dict[dip_type] = (pol_re, pol_im)

        # np.savetxt('uv.txt', uvspec, fmt='%.6e')
        # np.savetxt('im.txt', pol_im, fmt='%.6e')
        # np.savetxt('re.txt', pol_re, fmt='%.6e')

    if u.verbosity >= 1:
        print(" > DONE!")
        print
        print(u.banner(ch='#', length=80))
        print
