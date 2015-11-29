#!/usr/bin/env python

import numpy as np

import elements
import util as u
import parseinput as pi


def make_cent(cent, weight, coords):

    points = np.array([]).reshape(0,3)

    if weight > 0 and weight <= 1:
        weights = [weight] + [(1 - weight)/ (len(cent) - 1)] * (len(cent) - 1)
    else:
        weights = None

    for idx in cent:
        points = np.vstack((points, coords[idx]))

    return np.average(points, axis=0, weights=weights)


def make_dipo(dipole, dipole_types, coords):

    dip_type = dipole[0]
    dip_s = np.sqrt(dipole_types[dip_type][0])

    dip_ori = dipole[1:]

    # Orient dipole along an interatomic axis
    if len(dip_ori) == 2:
        
        a = coords[dip_ori[0]]
        b = coords[dip_ori[1]]
        direction = (b - a) / np.linalg.norm(b - a)

        dip = dip_s * direction
        
    # Orient dipole along an interatomic axis and forming and angle
    # theta with the plane defined by 3 atoms
    elif len(dip_ori) == 4:

        a = coords[dip_ori[0]]
        b = coords[dip_ori[1]]
        c = coords[dip_ori[2]]
        theta = dip_ori[3]

        # Generate a reference frame from point a, b, c
        # Our desired direction is b - a, i.e. the x axis
        # of this new ref frame.
        ref = u.refframe(a, b, c)
        cartesian = np.eye(3)

        # Define transformation matrices
        T = np.dot(cartesian, ref.T)
        Ry = u.rot_mat_y(theta)

        # Ry is 4x4, reduce it to 3x3
        Ry = Ry[:3,:3]

        x = cartesian[0]

        direction = np.dot(x, Ry)
        direction = np.dot(direction, T.T)
        
        dip = dip_s * direction


    else:

        print u.banner(text='ERROR', ch='#', length=80)
        print("Unrecognized number of parameters for dipole orientation.")
        print
        sys.exit()

    return dip


if __name__ == '__main__':

    # coords = pi.parse_PDB('struct.pdb')
    # coords = np.array(coords)
    # dipole_types, cent_dip_couples = pi.parse_input('devoe.in')

    # for cent, weight, dipoles in cent_dip_couples:
    #     appl_point = make_cent(cent, weight, coords)

    #     for dipole in dipoles:
    #         print make_dipo(dipole, dipole_types, coords)

    pass
