__author__ = 'ronen'

from Polymer import Polymer

from numpy import *

import scipy.optimize
from scipy.stats import norm
from matplotlib import pylab


def linear(x,a,b):
    return a*x+b

collect_mu_N = []
def draw_distribution(D,L,N,iterations,self_avoid, data_per_N):

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
    pylab.savefig("histogram_D=%d_N=%d_self_avoid=%d.pdf" % (D, N, self_avoid))


    global collect_mu_N
    collect_mu_N.append(param)

def simulate_sequence(D, L, iterations, min_size, max_size, step_size, self_avoid, reflect):

    poly = Polymer(D, L, self_avoid, reflect=reflect)
    data = []

    for N in range(min_size,max_size, step_size):
        print N #Poor's man progress bar
        data_per_N = []
        bad_number = 0
        for i in range(iterations):
            (last, count, grid) = poly.create_polymer(N)
            dist = linalg.norm(last - poly.start_point)

            if count == N:
                data_per_N.append(dist)
            else:
                bad_number += 1

            #Save some example realization.
            if i == 10 and D == 2:
                pylab.figure()
                pylab.imshow(grid, cmap=pylab.get_cmap("binary"))
                pylab.savefig("realization_N=%d_self_avoid=%d" % (N, self_avoid))

        bad_ratio = 1.0*bad_number/iterations
        if bad_ratio > 0.02:
            print "To many bads with N=%d. Bad ratio %.3f" % (N, bad_ratio)
        data_per_N = array(data_per_N)

        draw_distribution(D,L,N,iterations,self_avoid, data_per_N)


        R2 = (data_per_N**2).mean()
        data.append((R2, N))

    return data

def fit_result(data, D, self_avoid):
    a_data = array(data)
    x = log(a_data.T[1])
    y = log(a_data.T[0])

    print x,y
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
    parser.add_argument("--self_avoid", type=bool, nargs='?', default=False)
    parser.add_argument("--reflect", type=bool, nargs='?', default=False)
    parser.add_argument("--iterations", type=int, nargs='?', default=500)
    parser.add_argument("--min_size", type=int, nargs='?', default=50)
    parser.add_argument("--max_size", type=int, nargs='?', default=500)
    parser.add_argument("--step_size", type=int, nargs='?', default=50)



    args = parser.parse_args()

    data = simulate_sequence(args.D, args.L, args.iterations, args.min_size, args.max_size, args.step_size, args.self_avoid, args.reflect)
    print data
    fit_result(data, args.D, args.self_avoid)

