import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

# N is number of columns, M is number of rows
class LegendreFit:



    def __init__(self,N,M,order):

        """Initialize ZernikeFit object.

        Arguments:
            N -- first dimension of image
            M -- second dimension of image
            order -- Zernike order to fit up to
        """


        self.N = N
        self.M = M

        self.order = order

        # calculate number of terms based on Zernike order
        self.terms = (order+1)**2
        self.P = self.terms-1
        # calculate total size of image
        self.N0 = self.N*self.M
        # define coordinate system
        self.x1 = np.linspace(-1,1,M)
        self.y1 = np.linspace(-1,1,N)
        #self.x1, self.y1 = np.meshgrid(x1, y1)

        # initialize wavefront mapping coordinates
        self.x = self.x1
        self.y = self.y1

        # calculate Zernike polynomials on this grid
        leg_2d = self.get_legendre(self.x1,self.y1)


        #code for plotting Zernikes
        #for i in range(self.terms):
        #    zern2[str(i)] = np.reshape(self.zern[str(i)],(self.N,self.M))

        #for i in range(self.terms):
        #    plt.figure()
        #    plt.imshow(np.flipud(leg_2d[i]),cmap=plt.get_cmap('gray'))
        #    plt.show()

        # initialize Zernike matrices
        self.A = np.zeros((2*self.N0,self.P))
        self.A2 = np.zeros((self.N0,self.P))

        # calculate Zernike derivatives
        LegendreFit.make_A(self,self.x1,self.y1)

        

    
    #def set_zernikes(self,x1,y1,dx):

    #    """Method to re-generate zernikes onto a new grid size.
    #    Arguments:
    #        x1 -- x coordinates (2d grid)
    #        y1 -- y coordinates (2d grid)
    #        dx -- pixel size (m)
    #    """

    #    # rescale coordinates to the unit circle
    #    self.x = x1/dx*2/self.M
    #    self.y = y1/dx*2/self.N

    #    # calculate zernikes on the new grid
    #    self.zern = self.get_zernikes(self.x,self.y)

    def get_legendre(self,x1,y1):

        """Generate 2D Legendre polynomials
        Arguments:
            x1 -- x coordinates (1d grid)
            y1 -- y coordinates (1d grid)

        Returns:
            zern -- dictionary of zernike polynomials defined on the input grid.
        """

        # flatten coordinate arrays
        xf = x1.flatten()
        yf = y1.flatten()

        # initialize zernike dictionary
        leg_x = {}
        leg_y = {}

        leg_x[0] = np.ones(np.size(xf))
        leg_y[0] = np.ones(np.size(yf))
        leg_x[1] = xf
        leg_y[1] = yf

        # iterate through legendres
        if self.order>1:

            # start recurrence relation for Legendre polynomials
            for i in range(self.order-1):
                
                n = i+1

                leg_x[n+1] = ((2*n+1)*xf*leg_x[n] - n*leg_x[n-1])/(n+1)
                leg_y[n+1] = ((2*n+1)*yf*leg_y[n] - n*leg_y[n-1])/(n+1)

        leg_2d = {}

        for i in range(self.terms):
            
            nx = int(np.floor((i)/(self.order+1)))
            ny = int(np.mod((i),self.order+1))

            normx = 2./(2.*nx+1)
            normy = 2./(2.*ny+1)

            leg_2d[i] = (np.tile(leg_x[nx],(np.size(yf),1))*
                    np.tile(np.reshape(leg_y[ny],(np.size(yf),1)),(1,np.size(xf))))
                   #/normx/normy) 

        self.leg_x = leg_x
        self.leg_y = leg_y

        # populate Legendres into a matrix
        #self.l = np.zeros((x1.size*y1.size,self.terms - 1))
        #for i in range(self.terms - 1):
        #    self.l[:,i] = leg_2d[i+1].flatten()

        self.l = np.zeros((x1.size*y1.size,self.P))

        for i in range(self.P):
            self.l[:,i] = leg_2d[i+1].flatten() 

       
        return leg_2d


    
    def make_A(self,x1,y1):

        """Function to generate matrix for fitting wavefront gradient
        onto 2D Legendre basis
        Arguments:
            None
        Returns:
            None
        """

        xf = x1.flatten()
        yf = y1.flatten()

        # initialize zernike gradient dictionaries
        lx = {}
        ly = {}

        lx[0] = np.zeros(np.size(xf))
        ly[0] = np.zeros(np.size(yf))
        lx[1] = np.ones(np.size(xf))
        ly[1] = np.ones(np.size(yf))

        # loop through Legendre orders
        for i in range(self.order - 1):
            # calculate x and y gradient terms based on Legendre recurrence relations

            n = i+1

            lx[n+1] = 0.
            ly[n+1] = 0.

            num = int(np.floor((n+2)/2))
            for j in range(num):
                n2 = n - j*2

                norm = 2./(2.*n2+1.)

                lx[n+1] += 2.*self.leg_x[n2]/norm
                ly[n+1] += 2.*self.leg_y[n2]/norm


        # combine into one dictionary
        lgrad = {}
        for i in range(self.terms):
            i1 = i

            nx = int(np.floor((i)/(self.order+1)))
            ny = int(np.mod((i),self.order+1))
            
            normx = 2./(2.*nx+1)
            normy = 2./(2.*ny+1)

            lx_2d = (np.tile(lx[nx],(np.size(yf),1))*
                    np.tile(np.reshape(self.leg_y[ny],(np.size(yf),1)),(1,np.size(xf))))
                    #/normx/normy
            ly_2d = (np.tile(self.leg_x[nx],(np.size(yf),1))*
                    np.tile(np.reshape(ly[ny],(np.size(yf),1)),(1,np.size(xf))))
                    #/normx/normy)
            lgrad[i1] = np.append(lx_2d.flatten(),ly_2d.flatten())
            
        # temporary matrix A1, containing zernike gradient
        A1 = np.zeros((2*self.N0,self.P))

        # populate A1, columns are each zernike order, rows are pixel number
        for i in range(self.P):

            A1[:,i] = lgrad[i+1]

        # qr decomposition of A1, to get an orthonormal basis for the wavefront gradient
        A,r = np.linalg.qr(A1)

        # matrix for mapping orthonormal basis back onto Zernike basis
        self.mapping = np.linalg.inv(np.dot(np.transpose(A1),A))

        # set A as an object variable
        self.A = A

    def make_B(self,h_grad,v_grad):

        """Function to take gradient data inside the unit circle and make it consistent with basis
        Arguments:
            h_grad: horizontal gradient (2d)
            v_grad: vertical gradient (2d)
        Returns:
            B: 1d vector containing gradient data inside unit circle
        """

        h_flat = h_grad.flatten()
        v_flat = v_grad.flatten()
        B = np.zeros((2*self.N0,1))


        B[0:self.N0,0] = h_flat
        B[self.N0:,0] = v_flat

        return B

    def z_coeff_grad(self,h_grad,v_grad,dx,i_mask):

        """Function to project gradient onto Zernike coefficients
        Arguments:
            h_grad -- horizontal gradient (2d)
            v_grad -- vertical gradient (2d)
            dx -- pixel size
            i_mask -- amplitude-based mask to avoid fitting noise
        Returns:
            W -- zernike coefficients
        """

        # rescale gradient due to the fact that we normalized coordinates onto the unit circle
        h_grad = h_grad*dx*self.M/2
        v_grad = v_grad*dx*self.N/2

        # flatten amplitude mask to 1d array
        i_flat = i_mask.flatten()

        # tile the mask to account for the fact that the gradient has twice as many pixels as the amplitude
        i_flat = np.tile(i_flat,(2))

        # generate gradient vector
        B = self.make_B(h_grad,v_grad)

        # remove any area outside the amplitude mask
        B = B[i_flat,:]

        # remove any area outside the amplitude mask in the basis matrix
        A = self.A[i_flat,:]

        # projection onto orthonormal basis
        W0 = np.dot(np.transpose(A),B)

        # map back onto Zernike basis
        W = np.dot(np.transpose(self.mapping),W0)

        return W

    def wavefront_fit(self,W):

        """Function to calculate wavefront based on Zernike coefficients,
        onto whatever the current grid has been set to using the "set_zernikes" function
        Arguments:
            W -- Zernike coefficients
        Returns:
            wavefront -- wavefront (2d)"""

        # get grid shape
        M1 = np.size(self.x)
        N1 = np.size(self.y)

        # initialize wavefront
        #wavefront = np.zeros(self.x.size)
        #for i in range(self.P):
        #    wavefront += self.zern[str(i+1)]*W[i]

        wavefront0 = np.dot(self.l,W)

        wavefront = np.reshape(wavefront0,(N1,M1))
        return wavefront
