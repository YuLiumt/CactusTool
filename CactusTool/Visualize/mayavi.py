from mayavi import mlab
import numpy as np

# TODO:
def surf(dsets, rl):
    dset_init = dsets[0]

    for c in sorted(dset_init[rl]):
        origin = dset_init[rl][c]['origin']
        delta = dset_init[rl][c]['delta']
        mesh = dset_init[rl][c]['data']
        size = mesh.shape
        dim = 2
        coord = tuple(np.arange(0,size[(dim-1)-i])*delta[i]+origin[i] for i in range(dim))
        grid = np.meshgrid(*coord)
        mlab.surf(grid[1], grid[0], mesh, warp_scale='auto')   
    mlab.show()