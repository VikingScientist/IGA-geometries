from splipy import *
from splipy.io import G2, STL
import splipy.curve_factory as cf
import splipy.surface_factory as sf
import splipy.volume_factory as vf
import numpy as np
import re
import click

 
@click.command()
@click.option('--h',            default=.01,      help='Element size in all directions')
@click.option('--p',            default=2,        help='The polynomial degree')
@click.option('--layers',       default=3,        help='number of printable layers')
@click.option('--layer-height', default=.02,      help='layer height [m]')
@click.option('--print-width',  default=.04,      help='print width [m]')
@click.option('--print-length', default=0.5,      help='print length [m]')
@click.option('--filename',     default='out.g2', help='Output GoTools (g2) filename')
@click.option('--stl',          is_flag=True,     help='Create a stl output file')

def straight_wall(h, p, layers, layer_height, print_width, print_length, filename, stl):
    """ Creates a single straight wall piece (matematically a box) """
    nx = print_length / h # total number of length elements
    ny = print_width  / h # total number of width elements
    nz = layer_height / h # number of height elements per layer

    assert abs(nz-round(nz)) < 1e-10, 'Element size (h) does not divide layer-height'
    nz *= layers          # total number of height elements
    nx = round(nx)
    ny = round(ny)
    nz = round(nz)

    model = vf.cube((print_length, print_width, layer_height * layers))
    model.raise_order(p-1, p-1, p-1)
    model.refine(nx-1, ny-1, nz-1)

    with G2(filename) as f:
        f.write(model)
    print(f'\"{filename}\" written successfully')

    if stl:
        stl_filename = re.sub(r'\..*$', '.stl', filename)
        with STL(stl_filename) as f:
            f.write(model)
        print(f'\"{stl_filename}\" written successfully')

if __name__ == '__main__':
    straight_wall()
