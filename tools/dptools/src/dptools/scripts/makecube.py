#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Potential/charge data conversion to cube format.'''

import sys
import argparse
import numpy as np
import dptools.gridsio as gridsio
import dptools.grids as grids
import dptools.common as common

USAGE = """useage: makecube [options] INPUT OUTPUT

Converts the potential or charge data file from dftb+ transport
calculation into a VTK or Gaussian Cube format file.

Depending on the suffix of OUTPUT, the format can be guessed.
"""

def main(cmdlineargs=None):
    '''Main driver for makecube.
    
    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed. (Default: None)
    '''
    options = parse_cmdline_args(cmdlineargs)
    makecube(options)

def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.
    
    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed. (Default: None)
    '''
    
    parser = argparse.ArgumentParser(description=USAGE)
    
    helpstring = 'file containing X vector coordinates (default Xvector.dat)'
    parser.add_argument('--xvec', default='Xvector.dat',
                        help=helpstring, dest='xvec')
    
    helpstring = 'file containing Y vector coordinates (default Yvector.dat)'
    parser.add_argument('--yvec', default='Yvector.dat',
                        help=helpstring, dest='yvec')
    
    helpstring = 'file containing Z vector coordinates(default Zvector.dat)'
    parser.add_argument('--zvec', default='Zvector.dat',
                        help=helpstring, dest='zvec')
    
    helpstring = "Reference, If a present, the difference between the input" \
                 "file and the reference is printed"
    parser.add_argument('--reference', help=helpstring, dest='reference')
    
    helpstring = 'Output format (cube or vtk). If not specified, guessed' \
                 'from file extension'
    parser.add_argument('--format', help=helpstring, dest='format')

    parser.add_argument('infile', action="store")

    parser.add_argument('outfile', action="store")
    
    args = parser.parse_args()
    
    return args

def makecube(options):
    '''Converts a grided density or potential to VTK or cube format
    
    Args:
        infile: File containing the grided data
        outfile: Resulting file to write
        options: Options (e.g. as returned by the command line parser)
    '''
    
    # Build the grid
    with open(options.xvec, 'r') as xvecfile:
        xvec = np.array(xvecfile.read().split(), dtype=float)
    with open(options.yvec, 'r') as yvecfile:
        yvec = np.array(yvecfile.read().split(), dtype=float)
    with open(options.zvec, 'r') as zvecfile:
        zvec = np.array(zvecfile.read().split(), dtype=float)
    
    origin = (xvec[0], yvec[0], zvec[0])
    xres = xvec[1] - xvec[0]
    yres = yvec[1] - yvec[0]
    zres = zvec[1] - zvec[0]
    basis = ((xres, 0.0, 0.0), (0.0, yres, 0.0), (0.0, 0.0, zres))
    ranges = ((0, len(xvec)), (0, len(yvec)), (0, len(zvec)))
    grid = grids.Grid(origin, basis, ranges)
    
    # Build the data object
    sourcefile = common.openfile(options.infile)
    vec = np.array(sourcefile.read().split(), dtype=float)
    if options.reference is not None:
        referencefile = common.openfile(options.reference)
        ref = np.array(referencefile.read().split(), dtype=float)
        vec -= ref
    # Note: potential.dat is stored in the order x0y0z0...x0y0zn, x0y1z0...
    # which correspond to a row major order (right index changing faster)
    # in a A(x,y,z) notation
    griddata = grids.GridData(grid, vec.reshape(grid.shape, order='C'))
    
    if options.format is None:
        if options.outfile.endswith('.cub') or options.outfile.endswith('.cube'):
            outformat = 'cube'
        elif options.outfile.endswith('.vtk'):
            outformat = 'vtk'
    
    if outformat == 'vtk':
        gridsio.scalarvtk(options.outfile, griddata,
                          varname=options.outfile.split('.')[0])
    elif outformat == 'cube':
        gridsio.cube(options.outfile, griddata)
    else:
        raise ValueError('Cannot determine output format in '
                         'makecube:unknown file extension')
