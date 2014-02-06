
import numpy as np
import random


from main import linear
import scipy.optimize
from scipy.stats import norm

import pandas as pd

def normal_dist(x,mu,sigma):
    return (1.0/np.sqrt(2*np.pi*sigma**2))*np.exp(-(x-mu)**2/(2*sigma**2))


def normpdf(x, mu, sigma):
    u = (x-mu)/abs(sigma)
    y = (1/(np.sqrt(2*np.pi)*abs(sigma)))*np.exp(-u*u/2)
    return y


class Pivot:
    def __init__(self,D,L):
        self.D = D
        self.L = L



        self.rot = np.array([[0,1],[-1,0]])

        self.state = np.array([np.arange(L),np.zeros(L),np.zeros(L)],dtype=int).T

        self.possible_steps = list(np.load("/home/ronen/PycharmProjects/Polymer_random_walk/matrices.npy"))
        #self.possible_steps = [rot, rot*(-1), ref, ref*(-1)]
        #for s in self.possible_steps:
        #    print s
        #    print np.linalg.det(s)



    def plot_state(self):
        l_min = np.min(self.state,axis=0)
        l_max = np.max(self.state,axis=0)
        #print self.state
        #print l_min, l_max
        dim = (l_max - l_min)
        #print dim
        bias = l_min
        im = np.zeros(dim+1)
        for idx, l in enumerate(self.state):
            #print l-bias
            im[l[0]-bias[0],l[1]-bias[1],l[2]-bias[2]] = idx

        return im


    def self_intersect(self,state):

       new_len = len(set(tuple(i) for i in state.tolist()))
       if new_len == len(state):
           return False #No intersection
       return True


    def rotate(self,pivot, direction):
        rot =  direction
        base_loc = self.state[pivot-1][0]

        temp_state = self.state.copy()
        for idx, val in enumerate(self.state[pivot:]):
            temp_state[idx+pivot] = rot.dot(val- base_loc) + base_loc
        if self.self_intersect(temp_state):
            #print "Reject!"
            return False
        #    pass
        self.state = temp_state
        return True

    def create_polymer(self, steps):
        """
        Creates one realization of a polymer with 'steps' monomers.
        @param steps: The number of monomers
        @return: (final_point, step_number, grid)
         final point - Cartesian coordinates of R_N, the last monomer.
         step_number - The actual number of monomers.
         grid - A grid contains the current realization.
        """
        if steps>2**16:
            print "Too many steps"
            return



        i = 0
        dists = []
        rejection  = []
        while i<steps:
            i += 1

            #Choose Random direction.
            step = random.choice(self.possible_steps)

            pivot = np.random.randint(0,self.L,1)
            intersect = self.rotate(pivot, step)
            dist = np.linalg.norm(self.state[0] - self.state[-1])
            dists.append(dist)
            rejection.append(intersect)
        return np.array(dists), np.array(rejection)






def draw_distribution(D,N,iterations, data_per_N):

    #PDF stuff
    pylab.figure()
    pylab.title("%d" % N)
    param = norm.fit(data_per_N)
    x_hist = np.linspace(data_per_N.min(), data_per_N.max(), 500)
    pdf_fitted = norm.pdf(x_hist, loc=param[0], scale=param[1])
    xs = np.linspace(min(x_hist),max(x_hist),400)
    pdf_fitted2 = normpdf(xs, param[0],param[1])
    pylab.title("P(N=%d, R)" % N)
    pylab.plot(x_hist, pdf_fitted, 'r-', label="Gaussian fit, \n  $\\mu=%.3f, \\sigma=%.3f$" % param)
    pylab.plot(xs, pdf_fitted2)
    pylab.hist(data_per_N, bins=iterations / 20, alpha=.3, normed=1, label="data")
    pylab.xlabel("R")
    pylab.ylabel("Count (normelized)")
    pylab.legend()
    pylab.savefig("histogram_D=%d_N=%d_pivot.pdf" % (D, N))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Simulate polymer")
    parser.add_argument("-D", type=int, nargs='?', default=2)
    parser.add_argument("-L", type=int, nargs='?', default=500)
    parser.add_argument("--loop_avoid", type=bool, nargs='?', default=False)
    parser.add_argument("--reflect", type=bool, nargs='?', default=False)
    parser.add_argument("--iterations", type=int, nargs='?', default=5000)
    parser.add_argument("--discard", type=int, nargs='?', default=500)
    parser.add_argument("--min_size", type=int, nargs='?', default=50)
    parser.add_argument("--max_size", type=int, nargs='?', default=500)
    parser.add_argument("--step_size", type=int, nargs='?', default=50)

    parser.add_argument("--out", type=str, nargs='?', default=".")
    args = parser.parse_args()
    print args
    import os
    try:
        os.makedirs(args.out)
    except:
        print "dir exists"
    os.chdir(args.out)

    from matplotlib import pylab
    R2 = []
    rg = range(args.min_size,args.max_size,args.step_size)

    rejects = []
    dist_list = []
    dfs = []
    for N in rg:
        print N
        p = Pivot(3,N)
        dists, reject= p.create_polymer(args.iterations)
        dists = dists[args.discard:]
        reject = reject[args.discard:]
        df = pd.DataFrame({"dist":dists,"reject":reject})
        df["N"] = N
        df["count"] = N
        dfs.append(df)
        #draw_distribution(2,N,2000,dists)
        #dist_list.append(dists)
        #R2.append(np.mean(dists**2))
        #rejects.append(reject)
        #pylab.figure()
        #pylab.imshow(p.plot_state(), interpolation='None', cmap=pylab.get_cmap("binary"))
        #pylab.savefig("realization_%d.png" % N)


    #pylab.figure()

    #pylab.savefig("pivot_test.png")
    data = pd.concat(dfs)
    data.to_hdf("D=%d_L=%d_pivot_2D.hdf.pd" % (args.D,args.L ), "a")

    data["R^2"] = data["dist"] ** 2
    agg = data.groupby("count").mean()

    x = agg["R^2"].values
    y = agg.index.values

    x = np.log(x)
    y = np.log(y)

    popt, pcov = scipy.optimize.curve_fit(linear, x,y)
    extrapolated_y = linear(x,popt[0],popt[1])

    pylab.figure()
    pylab.plot(x,y, "o", label="data")
    pylab.plot(x,extrapolated_y,"-", label= "fit, %.3f  * x + %.3f" % tuple(popt))
    pylab.legend()
    pylab.savefig("fit.png")
    print popt[0]/2
