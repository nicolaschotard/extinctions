#!/usr/bin/env python

"""Plot the available extinction laws"""

import pylab as P
from Extinction import extinction

if __name__ == "__main__":

    print "\nUse extinction.py to plot the extinction laws.\n"

    EP = extinction.ExtinctionsPlots()

    # Plot all the avalaible extinction laws
    EP.plot_extinctionLaws()

    # Plot other useful figures
    EP.plot_cardelli_law()
    EP.plot_cardelli_law_variability()
    EP.plot_rlbd_variability()
    EP.plot_albd_variability()

    # See also EP.plot_all_figures()

    P.show()
