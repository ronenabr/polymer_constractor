__author__ = 'ronen'

from Polymer import Polymer

from numpy import *

import scipy.optimize
from scipy.stats import norm
from matplotlib import pylab


def linear(x,a,b):
    return a*x+b
import pandas as pd

collect_mu_N = []
def draw_distribution(D,L,N,iterations,loop_avoid, data_per_N):

    #PDF stuff
    pylab.figure()
    pylab.title("%d" % N)
    param = norm.fit(data_per_N)
    x_hist = linspace(data_per_N.min(), data_per_N.max(), 500)
    pdf_fitted = norm.pdf(x_hist, loc=param[0], scale=param[1])
    pylab.title("P(N=%d, R)" % N)
    pylab.plot(x_hist, pdf_fitted, 'r-', label="Gaussian fit, \n  $\\mu=%.3f, \\sigma=%.3f$" % param)
    pylab.hist(data_per_N, bins=iterations / 20, alpha=.3, normed=1, label="data")
    pylab.xlabel("R")
    pylab.ylabel("Count (normelized)")
    pylab.legend()
    pylab.savefig("histogram_D=%d_N=%d_loop_avoid=%d.pdf" % (D, N, loop_avoid))


    global collect_mu_N
    collect_mu_N.append(param)

def simulate_sequence(D, L, iterations, min_size, max_size, step_size, loop_avoid, reflect):

    poly = Polymer(D, L, loop_avoid, reflect=reflect)
    data = []

    for N in range(min_size,max_size, step_size):
        print N #Poor's man progress bar
        data_per_N = []
        bad_number = 0
        for i in range(iterations):
            (last, count, grid) = poly.create_polymer(N)
            dist = linalg.norm(last - poly.start_point)
            #print [dist,count,N]
            df = pd.DataFrame({"dist":[dist], "N": [N], "count" : [count]})
            data.append(df)

            #Save some example realization.
            if i == 10 and D == 2:
                pylab.figure()
                pylab.imshow(grid, cmap=pylab.get_cmap("binary"), interpolation="none")
                pylab.savefig("realization_N=%d_self_avoid=%d" % (N, loop_avoid))

        bad_ratio = 1.0*bad_number/iterations
        if bad_ratio > 0.02:
            print "To many bads with N=%d. Bad ratio %.3f" % (N, bad_ratio)


        #draw_distribution(D,L,N,iterations,loop_avoid, data_per_N)


        #R2 = (data_per_N**2).mean()
        #data.append((R2, N))

    data = pd.concat(data)
    return data

def fit_result(x,y, D, self_avoid):

    x = log(x)
    y = log(y)


    popt, pcov = scipy.optimize.curve_fit(linear, x,y)

    extrapolated_y = linear(x,popt[0],popt[1])

    pylab.figure()

    pylab.plot(x,y, "o", label="data")
    pylab.plot(x,extrapolated_y,"-", label= "fit, %.3f  * x + %.3f" % tuple(popt))

    pylab.legend(loc="best")

    pylab.xlabel(r"$\log \left\langle R^{2}\right\rangle$")
    pylab.ylabel(r"$\log N$")
    pylab.savefig("N_fit_D=%d_self_avoid=%d.pdf" % (D, self_avoid))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Simulate polymer")
    parser.add_argument("-D", type=int, nargs='?', default=2)
    parser.add_argument("-L", type=int, nargs='?', default=500)

    parser.add_argument('--loop_avoid',dest='loop_avoid',action='store_true')
    parser.add_argument('--no-loop_avoid',dest='loop_avoid',action='store_false')
    parser.set_defaults(loop_avoid=True)

    
    parser.add_argument('--reflect',dest='reflect',action='store_true')
    parser.add_argument('--no-reflect',dest='loop_avoid',action='store_false')
    parser.set_defaults(reflect=True)

    parser.add_argument("--iterations", type=int, nargs='?', default=500)
    parser.add_argument("--min_size", type=int, nargs='?', default=50)
    parser.add_argument("--max_size", type=int, nargs='?', default=500)
    parser.add_argument("--step_size", type=int, nargs='?', default=50)
    parser.add_argument("--out", type=str, nargs='?', default=".")




    args = parser.parse_args()
    args.loop_avoid = True
    print args
    import os
    try:
        os.makedirs(args.out)
    except:
        print "dir exists"
    os.chdir(args.out)

    data = simulate_sequence(args.D, args.L, args.iterations, args.min_size, args.max_size, args.step_size, True, args.reflect)
    data.to_hdf("D=%d_L=%d_loop_avoid=%d_reflect=%d.hdf.pd" % (args.D,args.L,args.loop_avoid,args.reflect ), "a")

    data["R^2"] = data["dist"] ** 2
    agg = data.groupby("count").mean()

    x = agg["R^2"].values
    y = agg.index.values

    fit_result(x,y, args.D, args.loop_avoid)

