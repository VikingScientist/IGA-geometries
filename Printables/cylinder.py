from splipy import *
from splipy.io import G2, STL
import splipy.curve_factory as cf
import splipy.surface_factory as sf
import splipy.volume_factory as vf
from splipy.utils import rotation_matrix
import numpy as np
import re
import click

 
@click.command()
@click.option('--nx',                    default=20,       help='Number of elements in lengthwise direction')
@click.option('--ny',                    default=1,        help='Number of elements in width direction')
@click.option('--nz',                    default=1,        help='Number of elements in height direction')
@click.option('--p',                     default=3,        help='The polynomial degree')
@click.option('--layers',                default=7,        help='number of printable layers')
@click.option('--layer-height',          default=.02,      help='layer height [m]')
@click.option('--print-width',           default=.02,      help='print width [m]')
@click.option('--print-length',          default=1.8,      help='print length = 2pi*radius [m]')
@click.option('--flat-floor',            default=True,     help='Flatten the bottom layer (degenerate face on patch 0)')
@click.option('--bulge-layers',          default=0,        help='Bulge each layer to get 3d-printed texture, ranging from 0 [%] (no bulge) to 100 [%] (degenerate elements)')
@click.option('--filename',              default='out.g2', help='Output GoTools (g2) filename')
@click.option('--stl',                   is_flag=True,     help='Create a stl output file')

def cylinder(nx, ny, nz, p, layers, layer_height, print_width, print_length, flat_floor, bulge_layers, filename, stl):
    """ Creates a cylinder wall piece """
    # nx = print_length / h # total number of length elements
    # ny = print_width  / h # total number of width elements
    # nz = layer_height / h # number of height elements per layer

    # assert abs(nz-round(nz)) < 1e-10, 'Element size (h) does not divide layer-height'
    # nz *= layers          # total number of height elements
    # nx = round(nx)
    # ny = round(ny)
    # nz = round(nz)

    radius = print_length / 2 / np.pi
    height = layer_height * layers

    # gives a B-spline approximation to the helix with arclength parametrization
    # unlike curve_factory.circle which is exact, but not arclength
    def arclength_circle(t):
        return np.array( [radius * np.cos(t), radius * np.sin(t), t/2/np.pi*layer_height] ).T
    helix = cf.fit(arclength_circle, 0, 2*np.pi)
    helix = helix.rebuild(p=p+1, n=nx+p+1)

    crossection = sf.square(size=(print_width, layer_height), lower_left=(radius,0)).rotate(np.pi/2, (1,0,0))

    assert 0 <= bulge_layers <= 100, 'Bulge layers need to be between 0% and 100%'
    if not bulge_layers == 0:
        crossection = Surface()
        crossection.raise_order(0,2)
        t = np.pi/2/100 * bulge_layers
        s = 1/3
        crossection[0,1,:] = [ -s*np.sin(t),   s*np.cos(t)]
        crossection[0,2,:] = [ -s*np.sin(t), 1-s*np.cos(t)]
        crossection[1,1,:] = [1+s*np.sin(t),   s*np.cos(t)]
        crossection[1,2,:] = [1+s*np.sin(t), 1-s*np.cos(t)]
        crossection *= (print_width, layer_height)
        crossection.rotate(np.pi/2, (1,0,0))
        crossection += (radius,0,0)

    one_layer = vf.revolve(crossection)
    one_layer.swap('w','v')
    one_layer.swap('u','v')
    one_layer.reverse('v')
    one_layer = one_layer.split( 0, 'u')
    for i in range(9):
        weights = one_layer[i,:,:,3]
        one_layer[i,:,:,2] += weights * layer_height*i/8
    one_layer = one_layer.rebuild((p+1,p+1,p+1), (nx+p, ny+p, nz+p))
    model = []
    for l in range(layers):
        model.append(one_layer + [0,0,l * layer_height])

    
    if flat_floor:
        for vol in model: vol -= [0,0,layer_height]
        bottom = model[1].faces()[4]
        flat = bottom.clone().project('xy')
        first_layer = vf.edge_surfaces(flat, bottom)
        first_layer.raise_order(0,0,p-1)
        first_layer.refine(0,0,nz-1)
        model[0] = first_layer

        
    with G2(filename) as f:
        f.write(model)
    print(f'\"{filename}\" written successfully')

    if stl:
        stl_filename = re.sub(r'\..*$', '.stl', filename)
        with STL(stl_filename, 2) as f:
            for vol in model:
                f.write(vol)
        print(f'\"{stl_filename}\" written successfully')

if __name__ == '__main__':
    cylinder()
