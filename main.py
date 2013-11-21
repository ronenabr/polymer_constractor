__author__ = 'ronen'

from Polymer import Polymer

from numpy import *

import scipy.optimize

from scipy.stats import norm




def linear(x,a,b):
    return a*x+b


polymer_num = 500

from matplotlib import pylab
D = 3
L = 300
self_avoid = True



poly = Polymer(D, L, self_avoid)
data = []

collect_mu_N = []

def draw_distribution(D,L,self_avoid, data_per_N):



    #PDF stuff
    pylab.figure()
    pylab.title("%d" % N)
    param = norm.fit(data_per_N)
    x_hist = linspace(data_per_N.min(), data_per_N.max(), 500)
    pdf_fitted = norm.pdf(x_hist, loc=param[0], scale=param[1])
    print param
    pylab.title("P(N=%d, R)" % N)
    pylab.plot(x_hist, pdf_fitted, 'r-', label="Gaussian fit, \n  $\\mu=%.3f, \\sigma=%.3f$" % param)
    pylab.hist(data_per_N, bins=polymer_num / 20, alpha=.3, normed=1, label="data")
    pylab.xlabel("R")
    pylab.ylabel("Count (normelized)")
    pylab.legend()
    pylab.savefig("histogram_D=%d_N=%d_self_avoid=%d.pdf" % (D, N, self_avoid))

    global collect_mu_N
    collect_mu_N.append(param)


for N in range(50,500,50):
    print N
    data_per_N = []
    bad_number = 0
    for i in range(polymer_num):
        (last, count, grid) = poly.create_polymer(N)
        dist = linalg.norm(last - poly.start_point)
        #print dist, count, last
        if count == N:
            data_per_N.append(dist)
        else:
            bad_number += 1
        if i == 10 and D == 2:
            pylab.figure()
            pylab.imshow(grid, cmap=pylab.get_cmap("binary"))
            pylab.savefig("realization_N=%d_self_avoid=%d" % (N, self_avoid))

    bad_ratio = 1.0*bad_number/polymer_num
    if bad_ratio > 0.02:
        print "To many bads with N=%d. Bad ratio %.3f" % (N, bad_ratio)
    data_per_N = array(data_per_N)

    draw_distribution(D,L,self_avoid, data_per_N)


    R2 = (data_per_N**2).mean()
    data.append((R2, N))

a_data = array(data)
x = log(a_data.T[1])
y = log(a_data.T[0])

popt, pcov = scipy.optimize.curve_fit(linear, x,y)

extrapolated_y = linear(x,popt[0],popt[1])

pylab.figure()

pylab.plot(x,y, "o", label="data")
pylab.plot(x,extrapolated_y,"-", label= "fit, %.3f  * x + %.3f" % tuple(popt))

pylab.legend(loc="best")

pylab.xlabel(r"$\log \left\langle R^{2}\right\rangle$")
pylab.ylabel(r"$\log N$")
pylab.savefig("N_fit_D=%d_self_avoid=%d.pdf" % (D, self_avoid))


