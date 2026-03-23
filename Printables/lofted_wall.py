from splipy import *
from splipy.io import G2, STL
import splipy.curve_factory as cf
import splipy.surface_factory as sf
import splipy.volume_factory as vf
import numpy as np
import re
import click

 
@click.command()
@click.option('--h',                   default=.01,      help='Element size in all directions')
@click.option('--p',                   default=3,        help='The polynomial degree')
@click.option('--layers',              default=11,       help='number of printable layers')
@click.option('--layer-height',        default=.02,      help='layer height [m]')
@click.option('--print-width',         default=.04,      help='print width [m]')
@click.option('--print-length',        default=0.5,      help='print length [m]')
@click.option('--horizontal-periods',  default=0.4,      help='number of horizontal periods sin(x) in length direction')
@click.option('--vertical-periods',    default=0.7,      help='number of veritcal periods cos(z) in height direction')
@click.option('--amplitude',           default=0.10,     help='amplitude of swings in the y-direction')
@click.option('--filename',            default='out.g2', help='Output GoTools (g2) filename')
@click.option('--stl',                 is_flag=True,     help='Create a stl output file')

def lofted_wall(h, p, layers, layer_height, print_width, print_length, horizontal_periods, vertical_periods, amplitude, filename, stl):
    """ Creates a loftede wall piece ( matematically sin(x)*cos(z) ) """
    nx = print_length / h # total number of length elements
    ny = print_width  / h # total number of width elements
    nz = layer_height / h # number of height elements per layer

    assert abs(nz-round(nz)) < 1e-10, 'Element size (h) does not divide layer-height'
    nz *= layers          # total number of height elements
    nx = round(nx)
    ny = round(ny)
    nz = round(nz)

    print_height = layer_height * layers

    sx = np.pi * 2 * horizontal_periods / print_length
    sz = np.pi * 2 * vertical_periods   / print_height

    surfs = []
    for l in range(layers):
        def math_curve(u):
            z = l*layer_height * np.ones(u.shape)
            return np.array([u, amplitude * np.sin(sx * u) * np.cos(sz * z), z]).T
        c1 = cf.fit(math_curve, 0, print_length)
        c2 = c1 + (0, print_width, 0)
        surfs.append(sf.edge_curves(c1,c2))

    model = vf.loft(surfs)
    model = model.rebuild(n=(nx+p,ny+p,nz+p), p=(p+1,p+1,p+1))

    with G2(filename) as f:
        f.write(model)
    print(f'\"{filename}\" written successfully')

    if stl:
        stl_filename = re.sub(r'\..*$', '.stl', filename)
        with STL(stl_filename, 2) as f:
            f.write(model)
        print(f'\"{stl_filename}\" written successfully')

if __name__ == '__main__':
    lofted_wall()
