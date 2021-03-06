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
import util as u

def parse_input(infile):
    '''Returns a Dictionary with dipole types,
    and a list of tuples for points of application of the dipoles.
    The keys of dipole types dictionary are the flags set ny the user in the input.
    Each value is a tuple. The first element of the tuple is a list containing the
    Dipolar Strength, the Excitation energy and the damping. The second element of
    the tuple is a flag defining the type of polarizability to calculate (electric or
    magnetic). The third element of the tuple if present only for magnetic polarizabilities
    and it is the bj coefficient.

    For application points, the data are tuples of 3 elements.
    Each tuple contains, as first element, the center expressed in terms
    of atom indexes and, as second element, the weight to be applied on the
    first element of the list. As third element, we have a list of dipoles to be applied
    in that point. The dipole is described through a list containing the dipole
    type and its orientation, in terms of atomic indexes.'''

    u.checkfile(infile)

    type_counter = 0
    center_counter = 0

    dipole_types = {}
    centers = []

    with open(infile) as f:

        while True:

            line = f.readline()

            if not line:
                break

            if line.strip():

                if line.startswith('#'):
                    continue
                
                # Create Dipole Types Dictionary
                elif line.split()[0].lower() == 'type':
                    
                    type_counter += 1
                    type_flag = line.split()[1]

                    # Try to read TYPE definition without optional keyword
                    try:
                        dipole = map(float, line.split()[2:])
                        pol_type = 'ele'
                        dipole_types[type_flag] = (map(float, line.split()[2:]), pol_type)

                    # Read TYPE definition with optional keyword
                    except ValueError:
                        pol_type = line.split()[-2]
                        bj = float(line.split()[-1])
                        dipole_types[type_flag] = (map(float, line.split()[2:-2]), pol_type, bj)

                    type_counter += 1
    
                # Centers
                elif line.split()[0].lower() == 'center':
    
                    center_counter += 1
                    weight = float(line.split()[1])
                    atom_idx = u.extend_compact_list(line.split()[2:])

                    # Adjust indices to Python's numeration
                    atom_idx = map(lambda x: x - 1, atom_idx)

                    data = f.readline().strip()

                    # Retrieve dipoles applied in this center
                    dipoles = []
                    while True:

                        if not data or data.split()[0].lower() != 'dipole':
                            break 
                        
                        data = data.split()
                        dip_type = [data[1]]
                        dip_ori = data[2:]
                        
                        if len(dip_ori) == 2:
                            dip_ori = u.extend_compact_list(dip_ori)

                        elif len(dip_ori) > 2:
                            theta = float(dip_ori[-1])
                            dip_ori = u.extend_compact_list(dip_ori[:-1])
                            dip_ori = dip_ori + [theta]


                        dipole = dip_type + dip_ori
                        dipoles.append(dipole)

                        data = f.readline().strip()

                    centers.append((atom_idx, weight, dipoles))

                # Dipoles
                elif line.split()[0].lower() == 'dipole':
    
                    if center_counter == 0:
                        print(u.banner(text='ERROR', ch='#', length=80))
                        print
                        print(" You defined a DIPOLE before assigning it a CENTER in %s" % infile)
                        print
                        sys.exit()
    
                    # dipole_counter += 1
    
    if type_counter == 0:
        print(u.banner(text='ERROR', ch='#', length=80))
        print
        print(" No dipole TYPE has been defined in %s" % infile)
        print
        sys.exit()


    return dipole_types, centers


def parse_PDB(pdb_file):

    u.checkfile(pdb_file)

    atom_names = []
    res_names = []
    res_ids = []
    atom_coord = []

    with open(pdb_file) as f:

        for line in f:

            # Read Atoms
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':

                atom_names.append(line[12:15])
                res_names.append(line[17:19])

                try:
                    res_ids.append(int(line[22:25]))

                except ValueError:
                    res_ids.append(None)

                x = float(line[30:37])
                y = float(line[38:45])
                z = float(line[46:53])
                atom_coord.append([x, y, z])

    return atom_coord


def parse_XYZ(xyz_file):

    u.checkfile(xyz_file)
    struct = np.loadtxt(xyz_file, skiprows=2, usecols=[1, 2, 3])

    return struct.tolist()


def parse_MOL2(mol2_file):

    u.checkfile(mol2_file)

    atom_names = []
    atom_types = []
    res_names = []
    res_ids = []
    atom_coord = []

    with open(mol2_file) as f:

        while True:

            line = f.readline()

            if not line:
                break

            # Read initial lines
            if line[0:17] == '@<TRIPOS>MOLECULE':
                f.readline()
                info = f.readline().split()
                NAtoms = int(info[0])
                try:
                    NRes = int(info[2])
                except:
                    NRes = 1

            # Read Atoms
            elif  line[0:13] == '@<TRIPOS>ATOM':
                for i in range(NAtoms):

                    data = f.readline().split()
                    
                    # save data for the old one
                    atom_names.append(data[1])
                    atom_types.append(data[5])
                    atom_coord.append(map(float, data[2:5]))
                    
    return atom_coord


if __name__ == '__main__':
    pass
