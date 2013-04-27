'''
Created on 29 mars 2013

Streamplots for triangular grids

http://reference.wolfram.com/mathematica/guide/VectorVisualization.html

note : numpy fromiter 
Create a new 1-dimensional array from an iterable object.
http://stackoverflow.com/questions/5559888/numpy-arrays-filling-and-extracting-data-quickly

lib/matplotlib/streamplot.py

def streamplot(axes, x, y, u, v, density=1, linewidth=None, color=None,
               cmap=None, arrowsize=1, arrowstyle='-|>', minlength=0.1): 

a Stream is a PolyCollection
a StreamSet is a silent_list of PolyCollection


cbook.silent_list('Line2D', self.lines)

class ContourSet(cm.ScalarMappable, ContourLabeler):

or maybe TriStreamSet(TriContour)

@author: Geoffroy
'''
from __future__ import division

import matplotlib.cbook as cbook
import matplotlib.cm as cm

import numpy as np
import math

from matplotlib.tri.tricontour import TriContourSet
import matplotlib.tri as mtri

from matplotlib.tri.triinterpolate import _prod_vectorized,\
    _scalar_vectorized, _transpose_vectorized



class ArrowSetter(object):
    """
    Mixin to provide arrow adding capacity to a StreamSet
    """
    def locate_arrows(self, stream):
        pass
    
    def add_arrows(self):
        pass
        
    def arrow(self, arrow_density):
        """
        Arrow density: dictionary - {'kind': {'user_loc', 'density'}
                                     'data': XX}
        """
        pass


class _StreamSet(cm.ScalarMappable, ArrowSetter):
    def __init__(self, ax, *args, **kwargs):
#        ContourSet.__init(self, ax, *args, **kwargs)
        pass


class TriStreamSet(_StreamSet):
    """
    Specialized StreamSet class for triangular grids
    """
    def color(self): pass
    def linewidth(self): pass


def tristreamplot(triangulation, u, v, minlength=0.1,
    stream_locator={'kind': 'auto', 'cells_grid': (25, 25), 'pts': None},
    linewidth=None, color=None, cmap=None, norm=None,
    arrow_locator={'kind': 'unique', 'loc': 0.5, 'size': 1., 'style': '->'},
    ):
    """
    Parameters
    ----------
    *triangulation* (npts, ntri)
    *u*, *v* (npts,)


    **kwargs similaire a streamplot
    *density*  -> ICI
    *linewidth* -> passe au constructeur TriStreamSet ? constante ou vecteur de taille npts
    *color* -> color constante ou vecteur de taille npts
    *cmap* -> ok
    *norm* -> ok
    *arrowsize*   ! problematic, to be dealt at TriStreamSet level
    *arrowstyle*  ! problematic, to be dealt at TriStreamSet level ? 
                # Add arrows half way along each trajectory.
    *minlength*  -> ??  Minimum length of streamline in axes coordinates.

    Returns
    -------
    TriStreamSet
"""
    # default to data coordinates
    if transform is None:
        transform = axes.transData

    # Ici  on aura des multiples lien collactions, 1 par streamline
    # (TriStreamSet)

    lc = mcollections.LineCollection(streamlines,
                                     transform=transform,
                                     **line_kw)


""" dtype to store a stream point """
stream_point = np.dtype([('pt', np.float64, (2,1)), ('cell', np.int32, 2)])



class _StreamIntegrator(object):
    """
    Integrates vector fields defined on a triangular grid.

    Parameters
    ----------
    *triangulation* (npts, ntri)
    *u*, *v* (npts,)
    *stream_locator* dictionary of keys ['kind', 'spacing_grid']
        - kind can be 'user', 'auto' (default to 'auto')
        - spacing_grid is a tuple (nx, ny) ; defauts to (25, 25) which will
          build a grid of 25 x 25 cells. Only 1 path is allowed to enter in each
          cell of the spacing grid.
    
    Attributes:
    -----------
    starts: np.array of type stream_point
        Points from where starts the stream plot
    cells: grid of mutually excluding cells, to avoid plotting too close

    Private attributes:
    -------------------
    _streams: _GrowingStreamsArray
        keeps in memory the growing streams
    _triangulation: mtri.Triangulation
    _u, _v: mtri.CubicTriInterpolator
    _x_range, _y_range: range of valid data
    _metric: dX = dot(metric, dx)  TODO
    
    """
#    head_state = np.dtype([('pt', np.float64, (2,1)),
#                           ('f', np.float64, (2,1)),
#                           ('df', np.float64, (2,2)),
#                           ('cell', np.int32)])

    def __init__(self, triangulation, u, v, stream_locator):
        self._triangulation = triangulation
        self._u = mtri.TriCubicInterpolator(triangulation, u)
        self._v = mtri.TriCubicInterpolator(triangulation, v)

        # We need to know the points used in the triangulation
        _, _x, _y = mtri.TriAnalyzer(
            triangulation)._get_compressed_triangulation()

        # We define some fictive reference coordinates x_ref, yref which fit a
        # unit square ; the transformation matrix from data coordinates
        # is diagonal (we only store the 2 diag terms).
        self._x_range = [np.min(_x), np.max(_x)]
        self._y_range = [np.min(_y), np.max(_y)]
        self._data_to_ref = [1. / (self.x_range[1] - self.x_range[0]),
                             1. / (self.y_range[1] - self.y_range[0])]

        # Speed grid... not needed as e can play with time scale !

        # Initializing the cells.
        ncells_x, ncells_y = self.stream_locator.get('cells_grid')
        self._cells = _CellsKeeper(
            (self._x_range[0], self._x_range[1], ncells_x),
            (self._y_range[0], self._y_range[1], ncells_y))
        cell_dx, cell_dy  = self._cells.dpts

        # Initializing the starting points.
        loc_kind = self.stream_locator.get('kind', 'auto')
        if loc_kind == "auto":
            start_pts = np.meshgrid(self.cells_x[:-1] + cell_dx * 0.5,
                                    self.cells_y[:-1] + cell_dy * 0.5).ravel()
        elif loc_kind == "user":
            start_pts = self.stream_locator['pts']
        else:
            raise ValueError(
                "In a stream_locator dictionary, `kind` value should be one"
                " of {'auto', 'user'}. Found: {0}".format(loc_kind))
        
        self._head = _HeadKeeper(start_pts)
        
#        np.empty(np.size(start_pts, 0), dtype=stream_point)
#        self.starts['pts'] = start_pts
#        self.starts['cells'] = self._cells.get_cell_index(start_pts)

    def integrate():
        """
        Computes the streampathes in data coordinates.
        """
        # Initialize the head
        self._streams = _GrowingStreamsArray(self.starts)
        n_streams = self._streams.count
        
        head_u, head_v = self.get_data_uv(self.starts['pts'])
        self._streams.store_head(
            runing=np.ones(n_streams, dtype=np.bool), pts=self.starts['pts'],
            cells=self.starts['cells'], data_u=head_u, data_v=head_v)
        
#        self._head = np.empty(n_streams)
#        self._errors = np.empty(n_streams, dtype=np.float64)
#        self._dt = np.ones(self._streams.count, dtype=np.float64)
        dt = self._streams.compute_dt(ds, metric=self._metric)

        # Grows the streams using Runge-Kutta integration, until termination
        # of all (termination occurs for point out of the grid of when
        # crossing another path (in the same cell).
        while np.any(self._streams.alive):
            pts = self.time_cycle(dt)
            eaters = self.get_eaters()
            terminated_grid
            self._streams.grow(pts, eaters, killed)

        # Cleaning the path, removing too short pathes
        self.finalize()

    def time_cycle(self):
        """
        A new time cycle on all still growing pathes

        A time cycle will be finished when all steps are < ds and
        error < error_max
        """
        alive = self._streams.alive

        # Loop on the alive streams until termination
        runing = alive.copy()

        while np.any(runing):
            dpts = self._streams.walk(runing, dt)
            ds = self.pt_norm(dpts)
            ds_ok = (ds > max_l)

            ok = (runing & ds_ok)
            new_pts = self._streams.head_pts[ok] + dpts[ok]
            new_tri_index[ok] = self._triangulation.get_trifinder()(
                new_pts[0, :], new_pts[1, :])
            pts_in[ok] = (new_tri_index[ok] != -1)

            data_u, data_v = get_data_uv(pts_in[ok], pts_in[ok])
            (u_new, _, _) = data_u
            (v_new, _, _) = data_v
            # TODO norm the error
            error_u = self.f_norm(f_new[ok] - (f[ok] + df[ok] * f[ok] * dt[ok]))
            converged = (error_u < error_crit)

            # Still runing : we have to define dt
            dt[~ds_ok] *= (max_l / ds[~ds_ok])
            dt[~converged] *= 0.85 * (error_crit / error_u[~converged])
            runing = ~ok or ~converged

        # Step finished runing: we have to store

        self.streams.grow(alive, )
            
            
            
            
        # Updating self._streams
        self._streams.kill()

    def taylor_2nd_order_step(self, running):
        """
        Since we have a cheap estimation for df, we will rely on it rather
        than using a 2nd order Runge-Kutta step.

        Error estimation based on comparison with f(t+dt).

        Updates:
        running: unconverged streams for this step.
        """
        (u, dudx, dudy) = head_u
        (v, dvdx, dvdy) = head_v
        
        
        ## This error is below that needed to match the RK4 integrator. It
        ## is set for visual reasons -- too low and corners start
        ## appearing ugly and jagged. Can be tuned.
        maxerror = 0.003

        ## This limit is important (for all integrators) to avoid the
        ## trajectory skipping some mask cells. We could relax this
        ## condition if we use the code which is commented out below to
        ## increment the location gradually. However, due to the efficient
        ## nature of the interpolation, this doesn't boost speed by much
        ## for quite a bit of complexity.
        #maxds = min(1. / dmap.mask.nx, 1. / dmap.mask.ny, 0.1)
        max_ds = 1. # length with a 'good' normalization...

        dx = (u + (dudx*u + dudy*v) * 0.5*dt) * dt
        dy = (v + (dvdx*u + dvdy*v) * 0.5*dt) * dt
        
        ds = self.line_length(dx, dy)

        ds = maxds
        stotal = 0
        xi = x0
        yi = y0
        xf_traj = []
        yf_traj = []
    
        while dmap.grid.within_grid(xi, yi):
            # GBY a mettre a la fin..
            #xf_traj.append(xi)
            #yf_traj.append(yi)
            try:
                k1x, k1y = f(xi, yi)
                k2x, k2y = f(xi + ds * k1x,
                             yi + ds * k1y)
            except IndexError:
                # Out of the domain on one of the intermediate integration steps.
                # Take an Euler step to the boundary to improve neatness.
                ds, xf_traj, yf_traj = _euler_step(xf_traj, yf_traj, dmap, f)
                stotal += ds
                break
            except TerminateTrajectory:
                break
    
            dx1 = ds * k1x
            dy1 = ds * k1y
            dx2 = ds * 0.5 * (k1x + k2x)
            dy2 = ds * 0.5 * (k1y + k2y)
    
            nx, ny = dmap.grid.shape
            # Error is normalized to the axes coordinates
            error = np.sqrt(((dx2 - dx1) / nx) ** 2 + ((dy2 - dy1) / ny) ** 2)
    
            # Only save step if within error tolerance
            if error < maxerror:
                xi += dx2
                yi += dy2
                try:
                    dmap.update_trajectory(xi, yi)
                except InvalidIndexError:
                    break
                if (stotal + ds) > 2:
                    break
                stotal += ds
    
            # recalculate stepsize based on step error
            if error == 0:
                ds = maxds
            else:
                ds = min(maxds, 0.85 * ds * (maxerror / error) ** 0.5)
    
        return stotal, xf_traj, yf_traj

    def clean():
        """
        Cleaning unreasonably short pathes.
        """
        pass
        
    def get_data_uv(self, pts, pts_tri_index=None):
        """
        TODO : should return f, df
        """
        if pts_tri_index is None:
            pts_tri_index = self._triangulation.get_trifinder()(
                pts[0, :], pts[1, :])
        data_u = self._u._interpolate_multikeys(
            pts[0, :], pts[1, :], tri_index=head_tri_index,
            return_keys=('z', 'dzdx', 'dzdy'))
        data_v = self._u._interpolate_multikeys(
            pts[0, :], pts[1, :], tri_index=head_tri_index,
            return_keys=('z', 'dzdx', 'dzdy'))
        return data_u, data_v
    
    def pt_norm(self, pts):
        """ Screen norm of a data path """
        pass
    
    def f_norm(self, pts):
        """ Screen norm of a data speed vec (u, v) ie dx/dt , dy/dt """
        pass


class _CellsKeeper(object):
    """ 
    Utility class to store information about locked cells.
    
    Parameters
    ----------
    cells_x, cells_y : 3-tuples
        (start, stop, num) that could be passed to np.linspace
    """
    def __init__(self, cells_x, cells_y):
        (self._xmin, self._xmax, self.nx) = cells_x
        (self._ymin, self._ymax, self.ny) = cells_y
        self.tokens = -np.ones([self._nx - 1, self._ny - 1], dtype = np.int32)
        self.delta_x = float(self._xmax - self._xmin)
        self.delta_y = float(self._ymax - self._ymin)
 
    def get_cell_index(self, pts):
        """ Get the index of the cells containing pts """
        ind_x = np.array(np.floor(
            (pts[:, 0] - self._xmin) / self.delta_x * self.nx,
            dtype=np.int32))
        ind_y = np.array(np.floor(
            (pts[:, 1] - self._ymin) / self.delta_y * self.ny,
            dtype=np.int32))
        return ind_x, ind_y, self.validate_index(ind_x, ind_y)

    def validate_index(self, ind_x, ind_y):
        return ((ind_x >= 0) & (ind_x < self.nx) &
                (ind_y >= 0) & (ind_y < self.ny))

    def set_tokens(self, pts, tokens):
        """ Sets the tokens at the cells containing pts"""
        ind_x, ind_y, valid = self.get_cell_index(pts)
        self.tokens[ind_x[valid], ind_y[valid]] = tokens[valid]

    def get_tokens(self, pts):
        """ Returns the tokens from the cells containing pts"""
        ind_x, ind_y, valid = self.get_cell_index(pts)
        if np.all(valid):
            return self.tokens[ind_x, ind_y]
        else:
            tokens = -np.ones(np.size(ind_x, 0), dtype=np.int32)
            tokens[valid] = self.tokens[ind_x[valid], ind_y[valid]]
            return tokens

class _StreamHead(object):
    """ 
    Utility class to store & update information about the heads of a streamset.
    
    stored info:
    data pts = (x, y)
    data_f = (u, v)
    data_df = 
    """
    def __init__(self, starts, triang, u_interp, v_interp):
        pass
    
    def dt(self, heads):
        """
        step dt to walk at least dl, based on f and df.
        """
        dt_f = dl / norm(self.f) 
        dt_df = 2 * dl**2 / norm_sq(self.df * self.f)
        return np.nanmin(dt_f, dt_df)

    def grow(self, heads):
        """
        
        """
        dx = dt * f + 0.5 * dt**2 * df * f
        
    
    def error(self, f, heads):
        pass

# Utility class implementing a double chained list in arrays, as numpy arrays.
# appending a value to all 
#
# Note that for a stream path implementation the data type should be able to 
# store (x, y, cell)
class _GrowingStreamsArray(object):
    """
    Utility class implementing an array of simultaneously growing streams.

    TODO  : store only index pointing to a 1-d growing array. implementation
    with a double list.
    
    Conventions: 
        - streams grow from their end in the forward direction and are deleted
          from their begining in the forward direction.
        - Once a stream is dead, it will never grow again.
        - Growing streams can 'eat' other streams by their tails

    Parameters:
    -----------
    starts: (x, y, cells) arrays



    The purpose is to optimise the speed of the following operations:
        - appending 1 item to all the active streams.
        - deleting "tails" (= all the points from a stream belonging to its
          first cell.)

    Private attributes:
    - _chuncksize: size of the unit chunck.
    - _n_chuncks: number of chuncks used.
    - _stream_data_array: list of size nchunck of ndarrays of size 
                       (n_streams[ichunck], chuncksize)
        data type of the array chuncks:
        (float64, float64, int32,            int32)
        (      x,      y,   cell,  prev_cell_index)
    - _chuncks_renum: list of size nchunk of ndarrays of size (n_streams,)
    - _streams_span int32 ndarray of size (n_streams, 2)
    - _streams_state bool ndarray of size (n_streams,) (alive: True, dead: False)
    - head properties ???

    Methods:
    - grow(streams, x, y, cells)
    - kill(streams) : kill streams
    - erase_tail(streams): erase the tail
    """
    stream_point = np.dtype([('x', np.float64), ('y', np.float64),
                             ('cell', np.int32)])

    def __init__(self, starts):
        self._streams = [[] for _ in starts]
        self._data = _GrowingArray(dtype=stream_point)
        

    def grow(self, streams, points):
        """ Grows the streams. """
        pos = self._data.push(points) 
        for i, stream in enumerate(streams):
            self._streams[stream].append(pos + i)

    def kill(self, streams):
        """
        Kills streams.
        """
        
    def untail(self, streams):
        """
        Erases the point froms the given streams belonging to the first cell.
        """

    def tails(self, streams):
        """
        Returns a view to the points of the current tail.
        """
        tail_index = self.span[:,0]
        

#    def heads(self):
#        """
#        Returns a view to the points of the current head.
#        """
#        head_index = self.span[:,1]
        
    def store_head(self, *args, **kwargs):
        """ Stores new head data """
        self._heads.store(*args, **kwargs)

    def compute_dt(self, *args, **kwargs):
        """ Computes time step """
        self._heads.compute_dt(*args, **kwargs)


class _StreamsArrayHead(object):
    """
    Heads of a _GrowingStreamsArray.

    Parameters:
    -----------
    n_heads: integer
        number of heads to track.
    """
    def __init__(self, n_heads):
        self.count = n_heads
        # Points
        self.pts = np.empty([n_heads, 2], dtype=np.float64)
        # Index of containing cells
        self.cells = -np.ones([n_heads], dtype=np.int32)
        # Values of f = [u, v] at points
        self.f = np.empty([2, 1], dtype=np.float64)
        # Differential of f at points
        self.df = np.empty([2, 2], dtype=np.float64)
        # Time_cycle
        self.cycles = -np.ones([n_heads], dtype=np.int32)
        # State
        self.alive = np.array([n_heads], dtype=np.bool)

    def store(self, runing, pts, cells, data_u, data_v):
        """
        Updates stream heads.
        """
        self.pts[runing, :] = pts
        self.pts[runing] = cells
        u, dudx, dudy = data_u
        v, dvdx, dvdy = data_v
        self.f[runing, 0] = u
        self.f[runing, 1] = v
        self.df[runing, 0, 0] = dudx
        self.df[runing, 1, 0] = v
        self.df[runing, 0, 1] = v
        self.df[runing, 1, 1] = dvdy

    def compute_dt(self, ds, metric):
        """
        """

    def walk(self, runing, dt):
        """
        Take a dt walk.
        
        Parameters:
        -----------
        runing: boolean array of size self.count
            "runing" heads.
        dt: 1d array (size: runing.sum)
            time steps
        
        Returns:
        --------
        dpts: 2d array of shape (runing.sum, 2)
            computed dx, dy 
        """
        f = self.f[runing, :, :]
        df = self.df[runing, :, :]
        dt2 = dt * dt
        dpts = (_scalar_vectorized(dt, f + _scalar_vectorized(
            0.5*dt, _prod_vectorized(df, f))))
        return dpts[:, :, 0]


class _GrowingArray(object):
    """ Toy implementation of a 1-d growing array """
    def __init__(self, dtype=np.float64):
        self.dtype = dtype
        self._size = 1
        self._pos = 0
        self._array= np.empty(self._size, dtype=self.dtype)

    def push(self, pushed_array):
        """ Pushes array pushed_array at end of self """
        pushed_array = np.asarray(pushed_array, dtype=self.dtype)
        (pushed_len,) = pushed_array.shape
        size = self._size
        while size < (self._pos + pushed_len):
            size *= 2
        if size != self._size:
            extended_array = np.empty(size, dtype=self.dtype)
            extended_array[:self._size] = self._array
            self._size = size
            self._array = extended_array
        self._array[self._pos:(self._pos + pushed_len)] = pushed_array
        self._pos += pushed_len
        return (self._pos - pushed_len)

    def __getitem__(self, key):
        return self._array[:self._pos][key]
