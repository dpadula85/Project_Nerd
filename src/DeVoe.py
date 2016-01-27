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

    parser.add_argument('-o', '--out', default='spec', type=str, help='''Prefix for output files.''')

    parser.add_argument('--min', default=16000, type=float, help='''Low energy limit in wavenumbers.''')

    parser.add_argument('--max', default=35000, type=float, help='''High energy limit in wavenumbers.''')

    parser.add_argument('--step', default=100, type=float, help='''Step for the spectral range.''')

    parser.add_argument('--nsteps', default=None, type=float, help='''Number of steps for the calculation.''')

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
    uv_outfile = args.out + '.OD.dat'
    cd_outfile = args.out + '.CD.dat'
    u.verbosity = args.verbosity

    nsteps = args.nsteps

    if nsteps:
        step = (args.max - args.min) / nsteps
        SpecRange = np.arange(args.min, args.max + 1, step)

    else:
        step = args.step
        SpecRange = np.arange(args.min, args.max + 1, step)
        nsteps = len(SpecRange)

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

    # if u.verbosity >= 2:
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
        print(" %s file format not recognized." % structure)
        print
        sys.exit()

    #
    # Get dipole types, application points and orientations in terms of atomic
    # indexes from the input file
    #
    A2cm = 1e-8
    coords = np.array(coords) * A2cm
    dipole_types, dipoles_info = pi.parse_input(infile)

    #
    # Transform the atomic index representation in coordinates and store the
    # application point coordinates and the direction vector in arrays
    #
    centers = np.array([]).reshape(0,3)
    orientations = np.array([]).reshape(0,3)
    pol_types = []

    if u.verbosity >= 1:
        print(" > ORIENTING DIPOLES...")
        print

    for cent, weight, dipoles in dipoles_info:

        #
        # Get application point coordinates
        #
        appl_point = m.make_cent(cent, weight, coords)

        for dipole in dipoles:

            #
            # Get the unit vector for the orientation
            #
            e = m.make_dipo(dipole, dipole_types, coords)
            
            # Translate the dipole in its application point
            # e += appl_point

            centers = np.vstack((centers, appl_point))
            orientations = np.vstack((orientations, e))
            pol_types.append(dipole[0])

        if u.verbosity >= 2:
            print

    # np.savetxt('locs.txt', centers, fmt='%10.4f')
    # np.savetxt('ori.txt', orientations, fmt='%10.4f')
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

    # if u.verbosity >= 2:
        print("   Lineshape                : %10s" % lineshape)
        print("   Low energy limit (cm-1)  : %10.2f" % args.min)
        print("   High energy limit (cm-1) : %10.2f" % args.max)
        print("   Step (cm-1)              : %10.2f" % step)
        print("   Number of Steps          : %10d" % nsteps)
        print

    pol_types_dict = {}
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

        # Save Polarizability in its complex form
        pol = np.c_[SpecRange, pol_re[:,1] + 1j * pol_im[:,1]]

        pol_types_dict[dip_type] = pol
        # np.savetxt('uv.txt', uvspec, fmt='%.6e')
        # np.savetxt('im.txt', pol_im * 1e24, fmt='%.6e')
        # np.savetxt('re.txt', pol_re * 1e24, fmt='%.6e')

    #
    # Build Interaction Matrix
    #
    if u.verbosity >= 1:
        print(" > BUILDING INTERACTION MATRIX...")

    #
    # We only need to build the G matrix elements with i != j and to iterate
    # on half of the matrix, diagonal excluded. Each element of the G matrix
    # is a 3x3 matrix.
    #
    Upper_half = []
    I = np.eye(3)
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):

            # Calculate the interaction
            # These distances are in Angstrom
            r_ij = centers[j] - centers[i]
            R_ij = np.linalg.norm(r_ij)

            if R_ij < 1e-25:

                G_ij = np.zeros((3,3))
                Upper_half.append(G_ij)

                if u.verbosity >= 2:
                    print("   Interaction between dipoles %d-%d turned off." % (i + 1, j + 1))

                continue

            else:


                # r_ij in dyadic form
                r_tensor = np.kron(r_ij, r_ij).reshape(3,3)

                G_ij = (I / R_ij**3) - (3 / R_ij**5) * r_tensor
                Upper_half.append(G_ij)

                if u.verbosity >= 2:
                   print("   Distance between dipoles %d-%d (Ang): %10.4f" % (i + 1, j + 1 , R_ij * 1e8))

        if u.verbosity == 1:
            u.progbar(i, len(centers))

    #
    # Generate a 3Nx3N matrix and reshape it to a matrix of 3x3 matrices. Then,
    # fill the upper triangle of G matrix with the 3x3 elements contained in Upper_half
    #
    G = np.zeros((len(centers)*3, len(centers)*3)).reshape(len(centers), len(centers), 3, 3)
    Uh = np.array(Upper_half)
    G[np.triu_indices(len(centers), 1)] = Uh

    #
    # Calculate optical spectra from Pols and G matrix
    #
    if u.verbosity >= 1:
        print
        print
        print(" > CALCULATING OPTICAL PROPERTIES...")

    G = G.astype(complex)
    uv_system = np.array([])
    cd_system = np.array([])

    #
    # For each frequency, put the polarizability tensors at that frequency on
    # the diagonal of G matrix
    #
    for k, freq in enumerate(SpecRange):

        if u.verbosity >= 1:
            u.progbar(k, len(SpecRange))

        for i, pol_type in enumerate(pol_types):

            pol_complex = pol_types_dict[pol_type][k][1]
            G[i,i] = I / pol_complex

        #
        # Reshape G matrix from a matrix of matrices to a regular 3Nx3N matrix
        # and symmetrize it
        #
        G_tmp = G.swapaxes(1,2).reshape(len(centers)*3,-1)
        G_tmp = G_tmp + G_tmp.T - np.diag(G_tmp.diagonal())
        # np.savetxt('G_im.dat', G_tmp.imag, fmt='%12.2e')
        # np.savetxt('G_re.dat', G_tmp.real, fmt='%12.2e')
        # sys.exit()

        #
        # Invert G matrix. Handle errors for inversion
        #
        try:
            A = np.linalg.inv(G_tmp)

        except np.linalg.linalg.LinAlgError as err:
            if 'Singular matrix' in err.message:
                print(u.banner(text='ERROR', ch='#', length=80))
                print(" Matrix Inversion did not work properly.")
                print(" Contact the author if this problem persists.")
                sys.exit()

        #
        # Check that the product between G and A gives the Identity matrix
        #
        if not np.allclose(np.dot(G_tmp, A).real, np.eye(len(centers)*3)):
            print(u.banner(text='ERROR', ch='#', length=80))
            print(" Matrix Inversion did not work properly.")
            print(" Contact the author if this problem persists.")
            sys.exit()

        #
        # Reshape the 3N*3N A matrix in a matrix of 3x3 matrices
        #
        # np.savetxt('A.dat', A.real, fmt='%8.2e')
        # A = np.array([ A[a*3:(a+1)*3, b*3:(b+1)*3] for (a, b) in np.ndindex(len(centers), len(centers)) ]).reshape(len(centers), len(centers), 3, 3)

        # blocks = map(lambda x : np.split(x, A.shape[1]/dim, 1), np.split(A, A.shape[0]/dim, 0))
        # A = np.array(blocks)
        # print A.real

        #
        # Calculate the contributes to the spectra for each matrix element of A
        #
        uv_freq = 0 
        cd_freq = 0 

        for m in range(A.shape[0]):

            for n in range(m, A.shape[1]):

                m_idx = m / 3
                e_m = orientations[m_idx]
                n_idx = n / 3
                e_n = orientations[n_idx]

                # Calculate ECD spectrum value
                r_mn = centers[m_idx] - centers[n_idx]
                C_mn = np.dot(r_mn, np.cross(e_m, e_n)) # - 4 * bj * np.dot(e_m, e_nmag) 

                # Calculate Absorption spectrum value
                if m == n:
                    uv_freq += A[m,n].imag * np.dot(e_m, e_n)
                    cd_freq += A[m,n].imag * C_mn

                else:
                    uv_freq += 2 * A[m,n].imag * np.dot(e_m, e_n)
                    cd_freq += 2 * A[m,n].imag * C_mn
                
        # CHECK THE CONSTANTS FOR UV and CD
        uv_freq = uv_freq * freq * 3 / cp.CGS_CNST2
        uv_system = np.r_[uv_system, uv_freq]

        cd_freq = cd_freq * freq**2 * cp.CGS_CNST3
        cd_system = np.r_[cd_system, cd_freq]

    uv_system = np.c_[SpecRange, uv_system]
    cd_system = np.c_[SpecRange, cd_system]

    # Save calculated spectra
    line = '%10.2f %10.6e'

    with open(uv_outfile, 'w') as uvfile:
        uvfile.write("# UV SPECTRUM\n")
        uvfile.write("# Generated with DeVoe.py\n")
        uvfile.write("# Lineshape : %s\n" % lineshape)
        uvfile.write("# cm-1     epsilon\n")
        np.savetxt(uvfile, uv_system, fmt=line)

    with open(cd_outfile, 'w') as cdfile:
        cdfile.write("# ECD SPECTRUM\n")
        cdfile.write("# Generated with DeVoe.py\n")
        cdfile.write("# Lineshape : %s\n" % lineshape)
        cdfile.write("# cm-1     Delta epsilon\n")
        np.savetxt(cdfile, cd_system, fmt=line)

    if u.verbosity >= 1:
        print
        print
        print(" > DONE!")
        print
        print(u.banner(ch='#', length=80))
        print
