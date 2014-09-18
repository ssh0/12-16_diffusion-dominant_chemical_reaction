#! /usr/bin/env python
# -*- coding:utf-8 -*-
#
# written by Shotaro Fujimoto, June 2014.

import numpy as np


def main_rw_d1(L, N, pmax):

    lattice = np.identity(L, dtype=int)
    T = [2 ** p for p in range(1, pmax + 1)]
    A = []
    random = np.random.rand
    np_sum = np.sum
    np_roll = np.roll

    for t in xrange(1, max(T) + 1):
        # 1 step ==============================================================
        for n in xrange(N):
            if np_sum(lattice[n]):
                if random() > 0.5:
                    lattice[n] = np_roll(lattice[n], 1)
                else:
                    lattice[n] = np_roll(lattice[n], -1)

        count = np_sum(lattice, axis=0)

        for l in xrange(L):
            if count[l] == 2:
                for n in xrange(N):
                    lattice[n][l] = 0
        # =====================================================================
        if t in T:
            A.append(L / (1. * np_sum(lattice)) - 1)

    return T, A


def plot_graph(x_data, y_data, x_labels, y_labels,
               xscale='linear', yscale='linear', aspect='auto'):
    """ Plot the graph about y_data for each x_data.
    """
    import matplotlib.pyplot as plt

    d = len(y_data)
    if not len(x_data) == len(y_data) == len(x_labels) == len(y_labels):
        raise ValueError("Arguments must have the same dimension.")
    if d == 0:
        raise ValueError("At least one data for plot.")
    if d > 9:
        raise ValueError("""So much data for plot in one figure.
                            Please divide two or more data sets.""")

    fig = plt.figure(figsize=(9, 8))
    subplot_positioning = [
        '11', '21', '22', '22', '32', '32', '33', '33', '33']
    axes = []
    for n in range(d):
        lmn = int(subplot_positioning[d - 1] + str(n + 1))
        axes.append(fig.add_subplot(lmn))

    for i, ax in enumerate(axes):
        ymin, ymax = min(y_data[i]), max(y_data[i])
        ax.set_aspect(aspect)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.set_xlabel(x_labels[i], fontsize=16)
        ax.set_ylabel(y_labels[i], fontsize=16)
        ax.set_ymargin(0.05)
        ax.plot(x_data[i], y_data[i], 'o-')

    fig.subplots_adjust(wspace=0.2, hspace=0.5)
    fig.tight_layout()
    plt.show()


def fitting(fit_func, parameter0, trial=10):
    import scipy.optimize as optimize

    a = np.zeros(trial)
    for i in range(trial):
        T, A = main_rw_d1(L, N, pmax)
        x, y = np.array(T), np.array(A)
        result = optimize.leastsq(fit_func, parameter0, args=(x, y))
        a[i] = result[0][0]

    print 'a =', np.average(a), 'sigma(a) =', np.std(a)

if __name__ == '__main__':

    L = 1000
    N = L
    pmax = 9

    T, A = main_rw_d1(L, N, pmax)
    plot_graph([T], [A], [r'$t$'], [r'$1/A(t)-1$'],
               xscale='log', yscale='log', aspect='equal')

    def fit_func(parameter0, t, y):
        a = parameter0[0]
        b = parameter0[1]
        residual = y - b * (t ** a)
        return residual
    parameter0 = [1, 0]  # a, b
#    fitting(fit_func, parameter0, trial=10)
