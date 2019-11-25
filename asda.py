import numpy as np
import warnings as wr
import matplotlib.pyplot as plt

from skimage import measure
from itertools import product
from scipy.interpolate import interp2d

__author__ = ['Jiajia Liu', 'Norbert Gyenge']
__email__ = ['jj.liu@sheffield.ac.uk']

class Asda_Calc:
    '''
    import asda
    lo = asda.Lamb_Oseen(vmax=2.0, rmax=50)
    vx, vy = lo.get_vxvy(x_range=[-100, 100], y_range=[-100, 100])
    lo.gamma_values()
    lo.center_edge()
    '''

    def __init__(self, vx, vy):

        '''
        as = asda(vx, vy, vx, vy)

        '''

        self.vx = vx

        self.vy = vy

    def reform2d(self, array, factor=1):
        ''' Missing DOC
        '''

        if factor > 1:

            congridx = interp2d(np.arange(0, array.shape[0]),
                                np.arange(0, array.shape[1]), array.T)

            array = congridx(np.arange(0, array.shape[0], 1/factor), 
                             np.arange(0, array.shape[1], 1/factor)).T

        return array

    def points_in_poly(self):
        """Return polygon as grid of points inside polygon.

        Input : poly (list of lists)
        Output : output (list of lists)
        """

        # Missing comment
        xs, ys = zip(*self.v)

        # Missing comment
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)

        # Missing comment
        newPoly = [(int(x - minx), int(y - miny)) for (x, y) in self.v]

        # Create mask by using skimage's measure
        mask = measure.grid_points_in_poly((round(maxx - minx) + 1,
                                            round(maxy - miny) + 1), newPoly)

        return [(x + minx, y + miny) for (x, y) in zip(*np.nonzero(mask))]

    def gamma_values(self, r=3, factor=1):
        '''
        Purpose: calculate gamma1 and gamma2 values of velocity field vx and vy
                 Fomulars can be found in Graftieaux et al. 2001
        Inputs:
            vx - velocity field in x direction
            vy - velocity field in y direction
            r - maximum distance of neighbour points from target point
            factor - default is 1. Magnify the original data to find sub-grid
                     vortex center and boundary
        Output:
            gamma - tuple in form of (gamma1, gamma2), where gamma1 is useful in finding
                  vortex centers and gamma2 is useful in finding vortex edges'''

        def Gen_Vel(i, j):
            '''   Missing documentation     '''

            vel = np.array([[self.vx[i+im, j+jm], self.vy[i+im, j+jm]]
                            for im in np.arange(-r, r+1)
                            for jm in np.arange(-r, r+1)])


            return np.array([vel, vel - vel.mean(axis=0)])

        def Calc_Gamma(pm, vel, pnorm, N):

            # Missing comment
            cross = np.cross(pm, vel)

            # Missing comment
            vel_norm = np.linalg.norm(vel, axis=2)

            # Missing comment
            sint = cross / (pnorm * vel_norm + 1e-10)

            return np.nansum(sint, axis=1) / N

        # Not sure why we need this
        self.vx = np.array(self.vx, dtype=np.float32).T
        self.vy = np.array(self.vy, dtype=np.float32).T

        # Check the dimensions of the velocity fields
        if self.vx.shape != self.vy.shape:

            raise Exception("Velocity field vx and vy do not match!") 
 
        if isinstance(r,int) is False:

            raise Exception("Parameter 'r' must be an integer!") 

        if isinstance(factor,int) is False:

            raise Exception("Parameter 'factor' must be an integer!") 

        if factor > 1:
            self.vx = reform2d(self.vx, factor=factor)
            self.vy = reform2d(self.vy, factor=factor)

        # Missing comment
        self.gamma = np.array([np.zeros_like(self.vx), np.zeros_like(self.vy)]).T

        # Missing comment
        pm = np.array([[i,j]
                      for i in np.arange(-r, r+1)
                      for j in np.arange(-r, r+1)], dtype=float)

        # Normalising array 'pm' 
        pnorm = np.linalg.norm(pm, axis=1)

        # Missing comment
        N = (2 * r + 1) ** 2

        # Create index array
        index=np.array([[i,j]
                        for i in np.arange(r, self.vx.shape[0]-r)
                        for j in np.arange(r, self.vy.shape[1]-r)]).T

        # Generate velocity field
        vel = Gen_Vel(index[0],index[1])

        # Iterate over the array gamma
        for dim, (i, j) in enumerate(product(np.arange(r, self.vx.shape[0]-r, 1),
                                             np.arange(r, self.vy.shape[1]-r, 1))):

            self.gamma[i, j, 0], self.gamma[i, j, 1] = Calc_Gamma(pm, vel[...,dim], pnorm, N)

        return self.gamma


    def center_edge(self, factor=1, rmin=4, gamma_min=0.89):
        '''
        Find vortices from gamma1, and gamma2
        Output:
            center: center location of vortices
            edge: edge location of vortices
            points: all points within vortices
            peak: maximum/minimum gamma1 value in vortices
            radius: equivalent radius of vortices
            All in pixel coordinates
        '''

        # Initial dictionary setup
        self.edge_prop={'center':0, 'edge':0, 'points':0, 'peak':0, 'radius':0}

        # Turn interactive plotting off
        plt.ioff()

        # Find countours
        cs = plt.contour(self.gamma[...,1].T, levels=[-2 / np.pi, 2 / np.pi])

        for i in range(len(cs.collections)):

            # Extract a contour and iterate over
            for c in cs.collections[i].get_paths():

                self.v = np.rint(c.vertices).tolist()
                ps = self.points_in_poly()

                dust = []

                for p in ps:
                    dust.append(self.gamma[...,0][int(p[0]), int(p[1])])

                if len(dust) > 1:

                    # Missing comment
                    re = np.sqrt(np.array(ps).shape[0] / np.pi) / factor

                    # allow some error around 0.9
                    # vortex with radius less than 4 pixels is not reliable
                    if np.max(np.fabs(dust)) >= gamma_min and re >= rmin:

                        # Extract the index, only first dimension
                        idx = np.where(np.fabs(dust) == np.max(np.fabs(dust)))[0][0]

                        # Update dictionary key 'center'
                        self.edge_prop['center'] += np.array(ps[idx])/factor

                        # Update dictionary key 'edge'
                        self.edge_prop['edge'] += np.array(self.v)/factor

                        # Update dictionary key 'points'
                        self.edge_prop['points'] +=  np.array(ps)/factor

                        # Update dictionary key 'peak'
                        self.edge_prop['peak'] += dust[idx]

                        # Update dictionary key 'radius'
                        self.edge_prop['radius'] += re

        return self.edge_prop

    def vortex_property(self, image=None):
        '''
        Calculate expanding, rotational speed, equivalent radius and average
            intensity of given swirls.
        Output:
            ve: expanding speed, pixel/frame
            vr: rotational speed, pixel/frame
            vc: velocity of the center, pixel/frame
            ia: average the observation values (intensity or magnetic field)
                within the vortices if image is given
        '''
        self.vx = np.array(self.vx)
        self.vy = np.array(self.vy)
        n_swirl = len(centers)
        ve = ()
        vr = ()
        vc = ()
        ia = ()

        for i in range(n_swirl):
            cen = centers[i]
            edg = edges[i]
            pnt = np.array(points[i], dtype=int)
            x0 = int(round(cen[0]))
            y0 = int(round(cen[1]))
            vcen = [self.vx[x0, y0], vy[x0, y0]]
            vc = vc + (vcen, )
            if image is not None:
                image = np.array(image).T
                value = 0
                for pos in pnt:
                    value = value + image[pos[0], pos[1]]
                value = value * 1.0 / pnt.shape[0]
            else:
                value = None
            ia = ia + (value, )
            ve0 = []
            vr0 = []
            for j in range(edg.shape[0]):
                idx = [edg[j][0], edg[j][1]]
                pm = [idx[0]-cen[0], idx[1]-cen[1]]
                tn = [cen[1]-idx[1], idx[0]-cen[0]]
                idx = np.array(idx, dtype=int)
                v = [self.vx[idx[0], idx[1]], self.vy[idx[0], idx[1]]]
                ve0.append(np.dot(v, pm)/np.linalg.norm(pm))
                vr0.append(np.dot(v, tn)/np.linalg.norm(tn))
            ve = ve + (np.nanmean(ve0), )
            vr = vr + (np.nanmean(vr0), )

        return (ve, vr, vc, ia)

    def visual_gamma(self, gamma_mean=True, fname=None, **kwargs):

        if gamma_mean:

            # Select data
            gamma = self.gamma[...,1]

            # Plot title:
            title = r'$\rho'
        else:

            # Select data
            gamma = self.gamma[...,0]

            # Plot title:
            title = r'$\rho'

        fig, ax = plt.subplots()

        # Show the image
        ax.imshow(gamma)

        # Set image title
        ax.set_title(title) 

        # Set axis labesl
        ax.set(xlabel='x', ylabel='y')

        if fname is None:
            plt.show()

        else:
            plt.savefig(fname, **kwargs)


class Lamb_Oseen(Asda_Calc):

    '''Creating an artifactual Lamb Oseen vortex

     Examples
    --------
    >>> import asda
    >>> lo = asda.Lamb_Oseen(vmax=2.0, rmax=50)

    * Generate grid

    >>> xx, yy = lo.get_grid(x_range=[-100, 100], y_range=[-100, 100])

    * Generate vx and vy

    >>> vx, vy = lo.get_vxvy(x_range=[-100, 100], y_range=[-100, 100])
    
    * Create a fancy matplotlib plot

    >>> lo.visual()'''

    def __init__(self, vmax=2.0, rmax=5, gamma=None, rcore=None, ratio_vradial=0):
        ''' Initialization of the Lamb Oseen vortex

        Parameters
        ----------
        vmax: `float`, optional
            maximum value of v_theta, negative value for clockwise vortex

        rmax: `float`, optional
            radius of the position where v_theta reaches vmax

        ratio_vradial: `float`, optional
            ratio between expanding/shrinking speed and rotating speed'''

        # Missing commend
        self.alpha = 1.256430

        # Missing comment   
        self.ratio_vradial = ratio_vradial

        if gamma is None or rcore is None:

            # Check if one of the input parameters is None but the other one is not None
            if (gamma is None) != (rcore is None):

                # Missing input parameter
                wr.warn("One of the input parameters is missing, setting both to 'None'")
                gamma, rcore = None, None
   
            # Radius of the position where v_theta reaches vmax
            self.rmax = rmax

            # Maximum value of v_theta
            self.vmax = vmax

            # Core radius
            self.rcore = self.rmax / np.sqrt(self.alpha)         

            # gamma is the circulation of the vortex in real world, gamma has unit of m^2 s^{-1}
            self.gamma = 2 * np.pi * self.vmax * self.rmax * (1 + 1/(2*self.alpha))

        else:

            # Radius of the position where v_theta reaches vmax
            self.rmax = self.rcore * np.sqrt(self.alpha)

            # Maximum value of v_theta
            self.vmax = self.gamma / (2 * np.pi * self.rmax * (1 + 1/(2*self.alpha)))

            # Missing comment
            self.rcore = rcore

            # Missing comment
            self.gamma = gamma


        # Missing comment
        self.vcore = (1 - np.exp(-1.0)) * self.gamma / (2 * np.pi * self.rcore)

    def get_grid(self, x_range, y_range):
        ''' Return meshgrid of the coordinate of the vortex

        Parameters
        ----------
        variable_name: `variable type`, optional or not optional?
            comment

        Return
        ------
        same here '''

        # Missing comment
        self.xx, self.yy = np.meshgrid(np.arange(x_range[0], x_range[1]),
                                       np.arange(y_range[0], y_range[1]))

        return self.xx, self.yy

    def get_vtheta(self, r=0):
        ''' Return v_theta at radius of r 

        Parameters
        ----------
        variable_name: `variable type`, optional or not optional?
            comment

        Return
        ------
        same here '''

        # Missing comment
        r = r + 1e-10

        return self.gamma *(1.0-np.exp(0-np.square(r)/np.square(self.rcore))) /(2*np.pi*r)


    def get_vradial(self, r=0):
        ''' Missing documentation 

        Parameters
        ----------
        variable_name: `variable type`, optional or not optional?
            comment

        Return
        ------
        same here '''

        # Missing comment
        r = r + 1e-10

        return self.get_vtheta(r) * self.ratio_vradial

    def get_vxvy(self, x_range, y_range, x=None, y=None):
        ''' calculate vx and vy value at point (x, y) 

        Parameters
        ----------
        variable_name: `variable type`, optional or not optional?
            comment

        Return
        ------
        same here '''

        # Check the dimensions of x_range
        if len(x_range) != 2:
            self.x_range = [0-self.rmax, self.rmax]

        # Check the dimensions of y_range
        if len(y_range) != 2:
            self.y_range = [0-self.rmax, self.rmax]

        if (x is None) or (y is None):

            # Check if one of the input parameters is None but the other one is not None
            if (x is None) != (y is None):

                # Missing input parameter
                wr.warn("One of the input parameters is missing, setting both to 'None'")
                x, y = None, None

            # Creating mesh grid
            x, y = self.get_grid(x_range=x_range, y_range=y_range)

        #Missing comment
        r = np.sqrt(np.square(x) + np.square(y)) + 1e-10

        #Missing comment
        vector = [0 - self.get_vtheta(r) * y + self.get_vradial(r) * x,
                  self.get_vtheta(r) * x + self.get_vradial(r) * y]

        self.vx = vector[0] / r
        self.vy = vector[1] / r
        

        return self.vx, self.vy

    def visual_vortex(self, fname=None, **kwargs):

        ''' calculate vx and vy value at point (x, y) 

        Parameters
        ----------
        variable_name: `variable type`, optional or not optional?
            comment

        Return
        ------
        same here '''
        
        fig, ax = plt.subplots()


        # Set image title
        ax.set_title('Lamb-Oseen Vortex') 

        # Generate a stream plot
        ax.streamplot(self.xx, self.yy, self.vx.T, self.vy.T, **kwargs)

        # Set axis labesl
        ax.set(xlabel='x', ylabel='y')

        if fname is None:
            plt.show()

        else:
            plt.savefig(fname, **kwargs)
