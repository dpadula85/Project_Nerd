#!/usr/bin/env python

import sys
import numpy as np
import util as u


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    # idxs = idxs.split()

    for idx in idxs:

        to_extend = idx.split('-')

        if len(to_extend) > 1:

            sel =  map(int, to_extend)
            extended += range(sel[0],sel[1]+1,1)

        else:
        
            extended.append(int(idx))
    
    return extended


def parse_input(infile):
    '''Returns a Dictionary with dipole types,
    and a list of tuples for points of application of the dipoles.
    Each tuple contains, as first element, the center expressed in terms
    of atom indexes and, ss second element, the weight to be applied on the
    first element of the list. As third element, we have a list of dipoles to be applied
    in that point. The dipole is described through a list containing the dipole
    type and its orientation, in terms of atomic indexes.'''


    u.checkfile(infile)

    type_counter = 0
    center_counter = 0
    dipole_counter = 0

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
                    dipole_types[type_flag] = map(float, line.split()[2:])
    
    
                # Centers
                elif line.split()[0].lower() == 'center':
    
                    center_counter += 1
                    weight = float(line.split()[1])
                    atom_idx = extend_compact_list(line.split()[2:])

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
                            dip_ori = extend_compact_list(dip_ori)

                        elif len(dip_ori) > 2:
                            theta = float(dip_ori[-1])
                            dip_ori = extend_compact_list(dip_ori[:-1])
                            dip_ori = dip_ori + [theta]


                        dipole = dip_type + dip_ori
                        dipoles.append(dipole)

                        data = f.readline().strip()

                    centers.append((atom_idx, weight, dipoles))

                # Dipoles
                elif line.split()[0].lower() == 'dipole':
    
                    if center_counter == 0:
                        print u.banner(text='ERROR', ch='#', length=80)
                        print
                        print("You defined a DIPOLE before assigning it a CENTER in %s" % infile)
                        print
                        sys.exit()
    
                    # dipole_counter += 1
    
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
