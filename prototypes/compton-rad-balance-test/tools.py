import numpy as np
#import scipy
#from defs import *


# simple trapedzoidal integration of an array
def integrate(xs, ys):

    # trapedzoidal
    #dxs = np.diff(xs)
    dxs = xs[1] - xs[0] # assume uniform grid
    return np.sum(ys*dxs)

    # simpsons composite formula
    #return scipy.integrate.simps(ys, dx=dxs)

    # romberg's method; requires N = 2**k + 1 grid
    #return scipy.integrate.romb(ys, dxs)


# 2D integration of an array
def integrate2d(xs, ys, zs):

    dxs = xs[1] - xs[0] #np.diff(xs)
    dys = ys[1] - ys[0] #np.diff(ys)

    # 2d trapedzoidal
    return dxs*dys*np.sum(zs)
    #return np.sum(dxs*dys*zs[1:-1, 1:-1] ) # only interior points

    # 2D simpsons
    #return scipy.integrate.simps( scipy.integrate.simps(zs, dx=dys), dx=dxs )

    # TODO precalculated weight array; even x
    #nx = Nx-3;
    #wx = ones(1,Nx);
    #wx(2:2:nx-1) = 4;
    #wx(3:2:nx-2) = 2;
    #wx(Nx-2:Nx) = [3 3 1];

    #ny = Ny-3;
    #wy = ones(Ny,1);
    #wy(2:2:ny-1) = 4;
    #wy(3:2:ny-2) = 2;
    #wy(Ny-2:Ny) = [3 3 1];

    #I = sum(sum(w(1:ny,1:nx).*z(1:ny,1:nx)))/9 + \
    #       (sum(sum(w(1:ny,nx:Nx).*z(1:ny,nx:Nx))) + \
    #               sum(sum(w(ny:Ny,1:nx).*z(ny:Ny,1:nx))))/8 + \
    #        sum(sum(w(ny:Ny,nx:Nx).*z(ny:Ny,nx:Nx)))*9/64



# weighted average of x using delta
def bac_average(d, x):
    r = np.zeros_like(x)
    r[:-1] = (1-d[:-1])*x[1:] + d[:-1]*x[:-1]
    return r

def fro_average(d, x):
    r = np.zeros_like(x)
    r[:-1] = d[:-1]*x[1:] + (1-d[:-1])*x[:-1]
    return r

def derivate(xs, ys):

    #--------------------------------------------------
    # second order
    dxs = xs[1] - xs[0] # assume uniform grid
    r = np.zeros_like(ys)
    r[:-1] = ys[1:] - ys[:-1]
    r[-1]  = ys[-1] - ys[-2]
    return r/dxs

    #--------------------------------------------------
    # fourth order
    dxs = xs[1] - xs[0] # assume uniform grid
    dys = np.zeros_like(ys)

    # middle regions
    fm1 = ys[1:-3]
    f0  = ys[2:-2]
    fp1 = ys[3:-1]
    fp2 = ys[4:  ]
    dys[2:-2] = (1*fm1 - 27*f0 + 27*fp1 - 1*fp2)/24

    # end points
    dys[0] = ys[1]-ys[0]
    dys[1] = ys[2]-ys[1]

    dys[-2]= ys[-1]-ys[-2]
    dys[-1]= ys[-1]-ys[-2]

    return dys/dxs


def find_sorted_nearest(val, ys):
    for j in range(1, len(ys)):
        if ys[j] >= val:
            return j-1
    else:
        return j


# read individual poitns from axis and draw a connecting line post plotting
def plot_connect_points(ax, **kwargs):
    d = ax.get_lines()
    x = [] #np.array(len(d))
    y = [] #np.array(len(d))

    #print('data ----------------', d)
    for i,elem in enumerate(d):
        #print(elem.get_xdata, elem.get_ydata)
        dd = elem.get_data(orig=True)
        #print('dd', dd)
        #print(len(dd))
        #print(np.shape(dd))

        if np.shape(dd)[1] == 2: # skip axhline objects
            continue

        #print('xy', elem.get_xdata(orig=True), elem.get_ydata(orig=True) )
        xv = elem.get_xdata(orig=True)
        yv = elem.get_ydata(orig=True)

        for j in range(len(xv)):
            x.append(xv[j])
            y.append(yv[j])
    #print(x)
    #print(y)

    ax.plot(x, y, **kwargs)
    return



#--------------------------------------------------
# storage class
class Storage:
    def __init__(self,):
        self.data = {}
        return

    def append(self, key, val):
        if not(key in self.data): self.data[key] = []
        self.data[key].append(val)
        return
            

    def save_csv(self, fname):
        fname += '.csv'

        #--------------------------------------------------
        # get size
        def_key = list(self.data.keys())[0]
        Nrow = len(self.data[def_key])
        Ncol = len(list(self.data.keys()))

        print('writing txt to file')
        print('row col', Nrow, Ncol)

        #--------------------------------------------------
        # initialize arrays
        wd = np.zeros((Nrow, Ncol))
        header = ""

        #--------------------------------------------------
        # save first these
        key_order = ['lap', 'out_in', 'theta', 'tau', 'ene_heat',]
        i = 0
        for key in key_order:
            if key in self.data:
                #print('adding key', key, 'with data', self.data[key])
                wd[:, i] = self.data[key]
                header += key + ","
                i += 1

        #--------------------------------------------------
        # then store the rest
        for key in self.data.keys():
            if not(key in key_order):
                wd[:, i] = self.data[key]
                header += key + ","
                i += 1

        header = header[:-1] # remove last trailing comma

        #--------------------------------------------------
        # finally save

        print(header)

        np.savetxt(fname, wd, 
                fmt = '%1.4e',
                header = header,
                delimiter = ',',
                )



        return
