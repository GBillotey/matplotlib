'''
Created on 2 avr. 2013

@author: Geoffroy
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import triinterpolate
import trirefine
import tristream

from numpy.testing import assert_array_equal, assert_array_almost_equal,\
    assert_array_less

def test_walk_stream_linear():
    dx, dy = (1.5, 0.6)
    nx, ny = (10, 10)
    x, y = np.meshgrid(np.linspace(0.0, dx * 1.0, nx),
                       np.linspace(0.0, dy * 1.0, ny))
    x = x.ravel()
    y = y.ravel()
#    print 'x', x
#    print 'y', y
    triang = mtri.Triangulation(x.ravel(), y.ravel())

    ax, ay = 1.5 * dx, 0.5 * dy
    ux = x*0. + ax
    uy = (y-0.75)*x*2.5# + ay

    vec_interp = triinterpolate._StreamIntegrator(triang, u=(ux, uy))

    x0 = 0.1 * np.ones(30) * dx #- dx * 0.5
    y0 = np.linspace(0.1, 0.8, 30) * dy#- dy * 0.5
    dt0 =  0.01 * np.ones(30)
    dt0 = None
    
    x1, y1 = x0, y0 # initialisaing

    plt.scatter(x0, y0, color='red')


    tri_index0 = vec_interp._trifinder(x0, y0)
    tris_pts0 = vec_interp._tris_pts[tri_index0]
    alpha0 = vec_interp._get_alpha_vec(x0, y0, tri_index0)
    u0=None
    du0=None
    for i in range(200):
        print "**** loop i:", i
        (tri_index1, alpha1, u1, du1, dt1, flag_stagnate, flag_frontier
         ) = vec_interp.walk_streams_local(tri_index0, alpha0, u0, du0, lmax=0.1,
                           dt0=dt0)
#        print 'flag_frontier', flag_frontier
#        print 'flag_stagnate', flag_stagnate
        if np.any(flag_stagnate): 
            print "stagnate"
        
        ru = ~flag_frontier & ~flag_stagnate
        x1, y1 = vec_interp.get_xy_vec(alpha1, tri_index1)
        dx, dy = x1 - x0, y1 - y0,
        
#        assert_array_almost_equal(dx, ax * dt1)
#        assert_array_almost_equal(dy, ay * dt1)
        plt.scatter(x1, y1)

 
        if np.all(~ru):
            #raise
            break

        alpha0 = alpha1[ru, :]
        tri_index0 = tri_index1[ru]
        u0 = u1[ru,:, :]
        du0 = du1[ru, :, :]
        x0, y0 = x1[ru], y1[ru]
        dt0 = dt1[ru]
    print "x1:\n", x1

    plt.triplot(triang, color='blue')
    plt.show()




def test_walk_stream_constant():
    dx, dy = (1.5, 0.6)
    nx, ny = (10, 10)
    x, y = np.meshgrid(np.linspace(0.0, dx * 1.0, nx),
                       np.linspace(0.0, dy * 1.0, ny))
    x = x.ravel()
    y = y.ravel()
#    print 'x', x
#    print 'y', y
    triang = mtri.Triangulation(x.ravel(), y.ravel())

    ax, ay = 1.5 * dx, 0.5 * dy
    ux = x*0.0 + ax
    uy = y*0.0 + ay

    vec_interp = triinterpolate._StreamIntegrator(triang, u=(ux, uy))

    x0 = np.array([0.1, 0.1, 0.1]) * dx #- dx * 0.5
    y0 = np.array([0.1, 0.5, 0.8]) * dy#- dy * 0.5
    dt0 = np.array([0.01, 0.02, 0.4]) * 0.1
    x1, y1 = x0, y0 # initialisaing

#    plt.scatter(x0, y0, color='red')


    tri_index0 = vec_interp._trifinder(x0, y0)
    tris_pts0 = vec_interp._tris_pts[tri_index0]
    alpha0 = vec_interp._get_alpha_vec(x0, y0, tri_index0)
    u0=None
    du0=None
    for i in range(200):
#        print "**** loop i:", i
        (tri_index1, alpha1, u1, du1, dt1, flag_stagnate, flag_frontier
         ) = vec_interp.walk_streams_local(tri_index0, alpha0, u0, du0, lmax=0.1,
                           dt0=dt0)
#        print 'flag_frontier', flag_frontier
#        print 'flag_stagnate', flag_stagnate
        if np.any(flag_stagnate): 
            raise
        
        ru = ~flag_frontier & ~flag_stagnate
        x1, y1 = vec_interp.get_xy_vec(alpha1, tri_index1)
        dx, dy = x1 - x0, y1 - y0,
        
        assert_array_almost_equal(dx, ax * dt1)
        assert_array_almost_equal(dy, ay * dt1)
#        plt.scatter(x1, y1)

 
        if np.all(~ru):
            #raise
            break

        alpha0 = alpha1[ru, :]
        tri_index0 = tri_index1[ru]
        u0 = u1[ru,:, :]
        du0 = du1[ru, :, :]
        x0, y0 = x1[ru], y1[ru]
        dt0 = dt0[ru]
    print "x1:\n", x1

#    plt.triplot(triang, color='blue')
#    plt.show()


def test_walk_to_border():
    
    for i_tri_config in range(12):
        itri1 = i_tri_config%3
        itri2 = i_tri_config//3

        tri1 = np.roll([0, 1, 3], itri1)
        tri2 = np.roll([1, 2, 3], itri2)
        unit_x, unit_y = 2.0, 3.0
        x = np.array([0., unit_x, unit_x,  0.]) - 1.0
        y = np.array([0., 0., unit_y, unit_y]) - 1.0

        if i_tri_config < 9:
            triangles = [tri1, tri2]
        else:
            triangles = [tri1]

        triang = mtri.Triangulation(x, y, triangles)
        ux = x*0.0 + unit_x
        uy = y*0.0 + unit_y
        vec_interp = triinterpolate._StreamIntegrator(triang, u=(ux, uy))

        x0 = np.array([0.1, 0.4, 0.7]) * unit_x - 1.0
        y0 = np.array([0.7, 0.4, 0.1]) * unit_y - 1.0
        dt0 = np.array([0.2, 0.1, 0.3])
        #plt.scatter(x0, y0, color='red')

        tri_index0 = vec_interp._trifinder(x0, y0)
        alpha0 = vec_interp._get_alpha_vec(x0, y0, tri_index0)

        # Computing dalpha
        u0 = vec_interp.interpolate_local(tri_index0, alpha0, return_key=('u'))
        dx = (u0 * dt0) / dt0
        dalpha = vec_interp.get_dalpha_vec(dx, tri_index0)
        (dt, alpha_new, tri_index_new, has_neigh
         ) = vec_interp.walk_to_border(alpha0, dalpha, tri_index0)

        x1, y1 = vec_interp.get_xy_vec(alpha_new, tri_index_new)

        assert_array_almost_equal(u0, 1.) # note that u0 is rescaled.
        assert_array_almost_equal(x1, np.array([0.2, 0.5, 0.8]) * unit_x - 1.0)
        assert_array_almost_equal(y1, np.array([0.8, 0.5, 0.2]) * unit_y - 1.0)
        assert_array_almost_equal(dt, [0.1, 0.1, 0.1])

#        plt.scatter(x1, y1)
#        plt.triplot(triang, color='blue')
#    plt.show()


def test_growing_array():
    growing = tristream._GrowingArray()
    n = 55
    tab = []
    for i in range(n):
        growing.push(np.ones(i+1) * (i+1))
        tab += [(i+1)*1.0] * (i+1)
    assert_array_equal(growing[:], tab)


    #print growing[:100:2]
#    print growing.array[:100]

def vectriinterpolate():
    # Test points within triangles of masked triangulation.
    x, y = np.meshgrid(np.arange(4), np.arange(4))
    x = x.ravel() * 1.33
    y = y.ravel() * 2.
    u = (1.23*x - 4.79*y, -2.56*x + 3.79*y)
    triangles = [[0, 1, 4], [1, 5, 4], [1, 2, 5], [2, 6, 5], [2, 3, 6],
                 [3, 7, 6], [4, 5, 8], [5, 9, 8], [5, 6, 9], [6, 10, 9],
                 [6, 7, 10], [7, 11, 10], [8, 9, 12], [9, 13, 12], [9, 10, 13],
                 [10, 14, 13], [10, 11, 14], [11, 15, 14]]
    mask = np.zeros(len(triangles))
    mask[8:10] = 1
    triang = mtri.Triangulation(x, y, triangles, mask)

    vec_interp = triinterpolate._StreamIntegrator(triang, u)
    xs = np.linspace(0.25, 2.75, 6)
    ys = [0.25, 0.75, 2.25, 2.75]
    xs, ys = np.meshgrid(xs, ys)
    xs = xs.ravel() * 1.33
    ys = ys.ravel() * 2.
    us, dus = vec_interp.interpolate(xs, ys, return_key='u_du')
    npts = us.shape[0]

    # Testing computed u
    u0 = np.empty([npts, 2, 1], dtype=np.float64)
    u0[:, 0, 0] = 1.23*xs - 4.79*ys
    u0[:, 1, 0] = -2.56*xs + 3.79*ys
    assert_array_almost_equal(us, u0)

    # Testing computed du
    du0 = np.empty([npts, 2, 2], dtype=np.float64)
    du0[:, 0, 0] = 1.23
    du0[:, 0, 1] = -4.79
    du0[:, 1, 0] = -2.56
    du0[:, 1, 1] = 3.79
    assert_array_almost_equal(dus, du0)

    # Testing points outside triangulation.
    xs = np.array([-0.25, 1.25, 1.75, 3.25])
    ys = xs
    xs, ys = np.meshgrid(xs * 1.33, ys * 2.)
    xs = xs.ravel()
    ys = ys.ravel()
    us = vec_interp.interpolate(xs, ys, return_key='u')
    assert_array_equal(us.mask, [[[True]*1]*2]*16)

    # Testing mixed configuration (outside / inside).
    xs = np.linspace(0.25, 1.75, 6)
    ys = [0.25, 0.75, 1.25, 1.75]
    xs, ys = np.meshgrid(xs * 1.33, ys * 2)
    xs = xs.ravel()
    ys = ys.ravel()
    us = vec_interp.interpolate(xs, ys, return_key='u')
    npts = us.shape[0]
    u0 = np.empty([npts, 2, 1], dtype=np.float64)
    u0[:, 0, 0] = 1.23*xs - 4.79*ys
    u0[:, 1, 0] = -2.56*xs + 3.79*ys
    assert_array_almost_equal(us, u0)
    mask = (xs >= 1.33) * (xs <= 2.66) * (ys >= 2.) * (ys <= 4.)
    assert_array_equal(us.mask[:, 0, 0], mask)
    assert_array_equal(us.mask[:, 1, 0], mask)


if __name__ == '__main__':
    test_growing_array()
    vectriinterpolate()
    test_walk_to_border()
    test_walk_stream_constant()
    #test_walk_stream_linear()
