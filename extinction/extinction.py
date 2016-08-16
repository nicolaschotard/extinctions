######################################################################
# Filename:      Extinction.py
# Version:       $Revision: 1.13 $
# Description:   Simple module providing Cardelli/O'Donnell law
# Author:        $Author: nchotard $
# $Id: Extinction.py,v 1.13 2015/05/21 11:50:37 nchotard Exp $
######################################################################

"""
.. _Extinction:

ToolBox.Astro.Extinction - Cardelli/O'Donnell extinction law
============================================================

References:

* CCM89: Cardelli, Clayton and Mathis
  (`<http://adsabs.harvard.edu/abs/1989ApJ...345..245C>`_)
* OD94: O'Donnell (`<http://adsabs.harvard.edu/abs/1994ApJ...422..158O>`_)
* FM98: Fitzpatrick & Massa (1998)
* G08: Goobar (`<http://adsabs.harvard.edu/abs/2008ApJ...686L.103G>`_)
"""

from __future__ import print_function

__author__ = "Yannick Copin <ycopin@ipnl.in2p3.fr>"
__version__ = "$Id: Extinction.py,v 1.13 2015/05/21 11:50:37 nchotard Exp $"

import numpy as N
from scipy import interpolate as I

# *Indicative* Bessel filter central wavelengths [A]
BesselFilters = dict(U=3590.,
                     B=4320.,
                     V=5410.,
                     R=6360.,
                     I=8000.)


def CCMextinctionParameters(lbda, odonnell=True):
    """
    Extinction law parameters a(lbda),b(lbda), to be used in A(lambda)/A(V) =
    a(lbda) + b(lbda)/Rv for lbda [A] in the 1000 A-33 microns range (for a
    V-filter centered on 5494.5 A). R(V):=A(V)/E(B-V) ranges from 2.75 to 5.3,
    Rv=3.1 being the canonical Milky Way value. Uses expressions from
    `Cardelli, Clayton and Mathis
    <http://adsabs.harvard.edu/abs/1989ApJ...345..245C>`_, with optical/NIR
    part updated from `O'Donnell
    <http://adsabs.harvard.edu/abs/1994ApJ...422..158O>`_.

    Set `odonnel` to False to use the original CCM89 expression.
    """

    def cardelliIR(x):
        """A/Av in InfraRed: 0.3-1.1 micron^-1"""

        assert ((x >= 0.3) & (x <= 1.1)).all()
        a = +0.574 * x ** 1.61
        b = -0.527 * x ** 1.61
        return a, b

    def cardelliOpt(x):
        """A/Av in Optical/Near IR: 1.1-3.3 micron^-1"""

        assert ((x >= 1.1) & (x <= 3.3)).all()
        y = x - 1.82
        pa = [+0.32999, -0.77530, +0.01979, +
              0.72085, -0.02427, -0.50447, +0.17699, 1]
        pb = [-2.09002, +5.30260, -0.62251, -
              5.38434, +1.07233, +2.28305, +1.41338, 0]
        a = N.polyval(pa, y)
        b = N.polyval(pb, y)
        return a, b

    def cardelliUV(x):
        """A/Av in UV: 3.3-8 micron^-1"""

        assert ((x >= 3.3) & (x <= 8)).all()
        y = x - 5.9
        fa = N.where(x >= 5.9, (-0.04473 - 0.009779 * y) * y ** 2, 0)
        fb = N.where(x >= 5.9, (+0.2130 + 0.1207 * y) * y ** 2, 0)
        a = +1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) + fa
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) + fb
        return a, b

    def cardelliFUV(x):
        """A/Av in Far UV: 8-10 micron^-1"""

        assert ((x >= 8) & (x <= 10)).all()
        y = x - 8.
        a = N.polyval([-0.070, +0.137, -0.628, -1.073], y)
        b = N.polyval([+0.374, -0.420, +4.257, +13.67], y)
        return a, b

    x = 1e4 / N.atleast_1d(lbda)          # [micron^-1]
    iIR = (x >= 0.3) & (x < 1.1)          # IR
    iOpt = (x >= 1.1) & (x < 3.3)         # Opt/Near IR
    iUV = (x >= 3.3) & (x < 8)            # UV
    iFUV = (x >= 8) & (x <= 10)           # Far UV

    a = N.zeros(len(x), 'd') * N.nan
    b = N.zeros(len(x), 'd') * N.nan
    a[iIR], b[iIR] = cardelliIR(x[iIR])
    if odonnell:
        a[iOpt], b[iOpt] = odonnellOpt(x[iOpt])  # Updated O'Donnell expression
    else:
        a[iOpt], b[iOpt] = cardelliOpt(x[iOpt])  # Original CCM89 expression
    a[iUV], b[iUV] = cardelliUV(x[iUV])
    a[iFUV], b[iFUV] = cardelliFUV(x[iFUV])

    return a, b

# O'Donnel94 extinction parameters
def odonnellOpt(x):
    """
    Parametrization from `O'Donnel 94
    <http://adsabs.harvard.edu/abs/1994ApJ...422..158O>`_ A/Av in Optical/Near
    IR: 1.1-3.3 micron^-1
    """

    if N.min(x) >= 3000:
        x = 10000. / N.array(x)  # Transform into micron^-1
    assert ((x >= 1.1) & (x <= 3.3)).all()
    y = x - 1.82
    pa = [-0.505, +1.647, -0.827, -1.718, +1.137, +0.701, -0.609, +0.104, 1]
    pb = [+3.347, -10.805, +5.491, +11.102, -7.985, -3.989, +2.908, +1.952, 0]
    a = N.polyval(pa, y)
    b = N.polyval(pb, y)
    return a, b


def extinctionLaw(lbda, law='OD94', Rv=3.1):
    """
    Return extinction law A(lbda)/Av for a given value of Rv and a given
    extinction law:

    * CCM89: `Cardelli, Clayton and Mathis
      <http://adsabs.harvard.edu/abs/1989ApJ...345..245C>`_
    * OD94: `O'Donnell <http://adsabs.harvard.edu/abs/1994ApJ...422..158O>`_
    * FM98: Fitzpatrick & Massa (1998)
    """

    if law == 'CCM89':
        a, b = extinctionParameters(lbda, odonnell=False)
    elif law == 'OD94':
        a, b = extinctionParameters(lbda, odonnell=True)
    elif law == 'FM98':
        a, b = FMExtinction(lbda, Rv=Rv), 0
    else:
        raise ValueError("'law' has to be either CCM89, OD94 or FM98.")
    return a + b / Rv


def extinctionFactor(lbda, ebmv, law='OD94', Rv=3.1):
    """
    Return extinction factor (in [0,1], actually more like a transmission
    factor) for E(B-V), based on an extinction law for Rv. See
    :func:`extinctionLaw` for available extinction laws.
    """

    a = extinctionLaw(lbda, law=law, Rv=Rv)      # A(lambda)/A_V
    return 10 ** (-0.4 * a * Rv * ebmv)   # A_V := R_V * E(B-V)


def FMExtinction(lbda, Rv=3.1, LMC2=False, AVGLMC=False, P=None,
                 k=5, s=None, CCMStyle=True):
    """
    The R-dependent Galactic extinction curve is that of Fitzpatrick & Massa
    (Fitzpatrick, 1999, PASP, 111, 63; astro-ph/9809387).  Parameterization is
    valid from the IR to the far-UV (3.5 microns to 0.1 microns).  UV
    extinction curve is extrapolated down to 912 Angstroms.

    * R_V - scalar specifying the ratio of total to selective extinction R(V) =
      A(V) / E(B - V).  If not specified, then R = 3.1 Extreme values of R(V)
      range from 2.3 to 5.3

    * AVGLMC - if set, then the default fit parameters c1,c2,c3,c4,gamma,x0 are
      set to the average values determined for reddening in the general Large
      Magellanic Cloud (LMC) field by Misselt et al.  (1999, ApJ, 515, 128)

    * LMC2 - if set, then the fit parameters are set to the values determined
      for the LMC2 field (including 30 Dor) by Misselt et al.  Note that
      neither /AVGLMC or /LMC2 will alter the default value of R_V which is
      poorly known for the LMC.

    The following five input keyword parameters allow the user to customize the
    adopted extinction curve.  For example, see Clayton et al. (2003, ApJ, 588,
    871) for examples of these parameters in different interstellar
    environments.

    * x0 - Centroid of 2200 A bump in microns (default = 4.596)
    * gamma - Width of 2200 A bump in microns (default  =0.99)
    * c3 - Strength of the 2200 A bump (default = 3.23)
    * c4 - FUV curvature (default = 0.41)
    * c2 - Slope of the linear UV extinction component
      (default = -0.824 + 4.717/R)
    * c1 - Intercept of the linear UV extinction component
      (default = 2.030 - 3.007*c2)

    .. NOTE::

       1. The following comparisons between the FM curve and that of Cardelli,
          Clayton, & Mathis (1989), (see ccm_unred.pro):

          a. In the UV, the FM and CCM curves are similar for R < 4.0, but
             diverge for larger R
          b. In the optical region, the FM more closely matches the
             monochromatic extinction, especially near the R band.

       2. Many sightlines with peculiar ultraviolet interstellar extinction
          can be represented with the FM curve, if the proper value of
          R(V) is supplied.
       3. Use the 4 parameter calling sequence if you wish to save the
          original flux vector.

    Based on `<http://idlastro.gsfc.nasa.gov/ftp/pro/astro/fm_unred.pro>`_
    """
    if type(lbda) in [int, float]:
        lbda = [lbda]

    def FMUV(x, Rv):
        xuv = (x - 5.9) * N.array((x - 5.9) > 0, dtype=int)
        yuv = P['c1'] + P['c2'] * x + P['c3'] * x ** 2 / \
              ((x ** 2 - P['x0'] ** 2) ** 2 + (x * P['gamma']) ** 2) + \
              P['c4'] * (0.5392 * xuv ** 2 + 0.05644 * xuv ** 3) + Rv
        return yuv

    def FMOP(Rv):
        xop = N.array([1.66667, 1.82815, 2.14132, 2.43309])
        yop = [N.polyval([2.13572e-04, 1.00270, -4.22809e-01], Rv),
               N.polyval([-7.35778e-05, 1.00216, -5.13540e-02], Rv),
               N.polyval([-3.32598e-05, 1.00184, 7.00127e-01], Rv),
               N.polyval([-4.45636e-05, 7.97809e-04, -5.46959e-03,
                          1.01707, 1.19456], Rv)]
        return xop, yop

    def FMIR(Rv):
        xir = N.array([0.37736, 0.81967])
        yir = N.array([0.26469, 0.82925]) * Rv / 3.1
        return xir, yir

    # Define parameters for the differentes cases
    if P != None:
        keys = ['gamma', 'x0', 'c4', 'c3', 'c2', 'c1']
        if not all([key in P for key in keys]):
            raise KeyError(
                "P dictionary must have the following keys: %s" % keys)
        else:
            print("Use the given set of parameters:")
    elif LMC2:
        P = {'gamma': 1.05, 'x0': 4.626, 'c4': 0.42,
             'c3': 1.92, 'c2': 1.31, 'c1': -2.16}
        print("Use the LMC2 set of parameters:")
    elif AVGLMC:
        P = {'gamma': 0.91, 'x0': 4.596, 'c4': 0.64,
             'c3': 2.73, 'c2': 1.11, 'c1': -1.28}
        print("Use the AVGLMC set of parameters:")
    else:
        c2 = -0.824 + 4.717 / Rv
        P = {'gamma': 0.99, 'x0': 4.596, 'c4': 0.41,
             'c3': 3.23, 'c2': c2, 'c1': 2.030 - 3.007 * c2}
        print("Use the standard set of parameters for MW:")
    print('; '.join(["%s=%.2f" % (p, P[p]) for p in P]))

    # Make sure that the law is constructed over a fair range of wavelength
    rlbda = N.arange(2000, 10000, 1)

    # Convert to inverse microns
    x = 10000. / N.array(rlbda)

    # Initialize the extinction curve array
    extcurve = x * 0.

    # Set the limits
    xcutuv = 3.704
    iuv, N_UV = x >= xcutuv, x < xcutuv
    iopir, Nopir = x < xcutuv, x >= xcutuv

    # Compute the NUV extinction
    if len(x[N_UV]) > 0:
        extcurve[iuv] = FMUV(x[iuv], Rv)

    # Compute the OP/IR extinction
    if len(x[Nopir]) >= 0:
        # For the UV (ref)
        xspluv = N.array([3.704, 3.846])
        yspluv = FMUV(xspluv, Rv)

        # For OP
        xsplop, ysplop = FMOP(Rv)

        # For IR
        xsplir, ysplir = FMIR(Rv)

        # Make a spline
        spline = I.UnivariateSpline(N.concatenate([[0], xsplir,
                                                   xsplop, xspluv]),
                                    N.concatenate(
                                        [[0], ysplir, ysplop, yspluv]),
                                    k=k, s=s)

        # Save only OP/IR
        extcurve[iopir] = spline(x[iopir])

    if CCMStyle:
        #normalization to a CCM law (=1 at lbd=5494.5)
        sp = I.UnivariateSpline(rlbda, extcurve/Rv, k=k, s=s)
        return sp(lbda)/sp(5494.5)
    else:
        sp = I.UnivariateSpline(rlbda, extcurve-Rv, k=k, s=s)
        return sp(lbda)


def FMUnreddened(lbda, flux, ebmv, Rv=3.1, LMC2=False,
                 AVGLMC=False, k=5, s=None):
    """
    Deredden a flux vector using the Fitzpatrick (1999) parameterization

    :param float rx1y: the correlation coefficient between variables x1 and y
    :param float rx2y: the correlation coefficient between variables x2 and y
    :param float rx1x2: the correlation coefficient between variables x1 and x2
    :param int n: sample size
    :return: two-tailed p-value

    :param lbda: wavelength vector (Angstroms)
    :param flux: calibrated flux vector, same number of elements as lbda If
      only 3 parameters are supplied, then this vector will updated on output
      to contain the dereddened flux.
    :param ebmv: color excess E(B-V), scalar.  If a negative E(B-V) is
      supplied, then fluxes will be reddened rather than dereddened.
    :return: unreddened flux vector, same units and number of elements
      as `flux`
    """

    extcurve = FMExtinction(lbda, Rv=3.1, LMC2=LMC2, AVGLMC=AVGLMC, k=k, s=s)
    return flux * 10. ** (0.4 * ebmv * extcurve * Rv)


# Goobar08 extinction law
def goobar08_law(lbd, lbdref, a, p):
    """
    From `Goobar08 <http://adsabs.harvard.edu/abs/2008ApJ...686L.103G>`_, *Low
    R_V from Circumstellar Dust around Supernovae*
    """
    return 1. - a + a * (lbd / lbdref) ** p


def ap_to_Rv(a, p, r=0.8):
    """
    a and p are the parameters of the goobar08 extinction law r is the ratio
    between the mean wavelegnths of the B and the V filter
    """
    return 1. / (a * (r ** p - 1.))


class ExtinctionsPlots:

    """
    Visualization of the extinction law and its evolutions with Rv or E(B-V).

    Example:

    >>> EP = ExtinctionsPlots()
    >>> EP.plot_all_figures()

    Or::

      python Extinction.py.
    """

    def __init__(self, dpi=100, wmin=3000, wmax=10000):

        self.pylab = __import__('pylab')  # Do not import in the library per se

        # Set the different parameters used to create the figures
        self.dpi = dpi
        self.wavelength = N.linspace(wmin, wmax, wmax - wmin)

        # Set the reference wavelength to the one of the CCM law (i.e. 5494.5
        # A)
        self.ref_wavelengths = {'B': 4394.43, 'V': 5494.51}
        self.cmap = self.pylab.cm.RdYlBu
        self.num = 100
        self.Rvs = N.linspace(1.1, 4.1, self.num)
        self.Ebmvs = [0.05, 0.2, 0.6]

    def plot_all_figures(self):
        """
        Plot all the available figures.
        Each figure are created by a self.plot_* method.
        See self.plot_functions.
        """

        f = [m for m in dir(self)
             if m.startswith('plot')
             and m != 'plot_all_figures']
        for m in f:
            getattr(self, m)()

    def plot_extinctionLaws(self):
        """
        Plots all the avalaible extinction laws from Extinction.py
        - The Goobar extinction law
        - The Cardelli extinction law
        - The Fitzpatrick extinction law
        - The O'Donnel extinction law
        """

        fig = self.pylab.figure(dpi=100)
        print(("Figure %i: Extinction laws" % (fig.number)).center(80, '='))
        print("""
        Avalaible extinction laws from Extinction.py
        """)
        gl = "Goobar law".center(30, '-')

        # Compute the goobar extinction law
        a, p = 0.9, -1.5  # Set the values for the MW dust
        Gext = goobar08_law(self.wavelength, self.ref_wavelengths['V'], a, p)
        print("""
        %s
        Set the values for the MW dust:
        a = %.2f ; p = %.2f
        Corresponding: Rv =%.2f
        """ % (gl, a, p, ap_to_Rv(a, p)))

        # Compute the CCM extinction law
        cl = "CCM law".center(30, '-')
        print("""
        %s
        CCM law for Rv = 3.1
        CCM law for Rv = %.1f
        """ % (cl, ap_to_Rv(a, p)))
        CCMext = extinctionLaw(self.wavelength, law='CCM89')
        CCMgoob = extinctionLaw(
            self.wavelength, Rv=ap_to_Rv(a, p), law='CCM89')

        # Compute the FM extinction law
        fl = "Fitzpatrick law".center(30, '-')
        print("""
        %s
        Fitzpatrick law for Rv = 3.1
        """ % (fl))
        FMext = FMExtinction(self.wavelength)

        # Compute the O'Donnel extinction law
        ol = "O'Donnel 94 law".center(30, '-')
        print("""
        %s
        O'Donnel 94 law for Rv = 3.1\n
        """ % ol)
        ODext = extinctionLaw(self.wavelength, law='OD94')

        # Make a figure
        ax = fig.add_axes([0.08, 0.09, 0.88, 0.86],
                          xlabel=r'wavelength [$\AA$]',
                          ylabel=r'$A_{\lambda}$/$A_{V}$',
                          title='Extinction laws')
        ax.plot(self.wavelength, CCMext, '-k',
                lw=3, label='CCM 89, Rv=3.1')
        ax.plot(self.wavelength, ODext, '--r',
                lw=3, label="O'Donnel 94, Rv=3.1")
        ax.plot(self.wavelength, FMext, '-r',
                lw=3, label='Fitzpatrick 99, Rv=3.1')
        ax.plot(self.wavelength, Gext, '-.g',
                lw=3, label="Goobar 08 (MW), 'Rv'=%.2f" % ap_to_Rv(a, p))
        ax.plot(self.wavelength, CCMgoob, '-.k',
                lw=3, label='CCM 89, Rv=%.2f' % ap_to_Rv(a, p))
        ax.legend(loc='best').draw_frame(False)
        ax.set_ylim(ymin=0.3, ymax=2.4)

    def plot_cardelli_law(self, Rv=3.1, ebmv=0.3):
        """
        Plot the cardelli extinction law for a given Rv, and its constituants
        a and b. Plot also the interstellar medium transmission for the given
        Rv and E(B-V).
        """
        fig = self.pylab.figure(dpi=self.dpi)
        print(("Figure %i: Cardelli extinction law" %
               (fig.number)).center(80, '='))
        print("""
        Top panel: Cardelli extinction law parammeters a and b (top pannel).
        Bottom panel:
        - CCM extinction law as A(lbd)/Av = a+b/Rv for Rv = %.2f
        - Interstellar medium transmission for Rv = %.2f and E(B-V) = %.2f.\n
        """ % (Rv, Rv, ebmv))

        # Create the figure and axes
        ax1 = fig.add_axes([0.06, 0.51, 0.92, 0.42],
                           title='CCM extinction law '
                           '($R_V=%.1f, E(B-V)=%.1f$)' % (Rv, ebmv))
        ax2 = fig.add_axes([0.06, 0.09, 0.92, 0.42],
                           xlabel=r'$\lambda$ [$\AA$]')

        # Get the extinction parameters, law and factor (transmission)
        a, b = extinctionParameters(self.wavelength)
        extLaw = extinctionLaw(self.wavelength, Rv=Rv)
        extFactor = extinctionFactor(self.wavelength, ebmv, Rv=Rv)

        # Plot them, adn set the labels
        ax1.plot(self.wavelength, a, 'y', lw=1.5, label=r'$a_{\lambda}$')
        ax1.plot(self.wavelength, b, '-.g', lw=1.5, label=r'$b_{\lambda}$')
        ax2.plot(self.wavelength, extLaw, '--b', lw=1.5,
                 label=r'$A_{\lambda}/A_{V}$ ($=a_{\lambda}+b_{\lambda}/R_V$)')
        ax2.plot(self.wavelength, extFactor, 'r', lw=1.5,
                 label=r'$T_{\lambda}$ ($=10^{-0.4 \times R_V \times '
                 'E(B-V) \times A_{\lambda}/A_V }$)')

        # Add lines and annotations about reference filter postions
        ax1.axvline(self.ref_wavelengths['B'], ls=':', color='k', lw=0.8)
        ax1.axvline(self.ref_wavelengths['V'], ls=':', color='k', lw=0.8)
        ax2.axvline(self.ref_wavelengths['B'], ls=':', color='k', lw=0.8)
        ax2.axvline(self.ref_wavelengths['V'], ls=':', color='k', lw=0.8)
        ax1.axhline(0, ls=':', color='k', lw=0.8)
        ax1.axhline(1, ls=':', color='k', lw=0.8)
        ax2.axhline(1, ls=':', color='k', lw=0.8)
        ax1.annotate(r'$\lambda_{B}$',
                     xy=(self.ref_wavelengths['B'] + 10, 3.3),
                     xycoords='data', color='k',
                     horizontalalignment='left', size='large')
        ax1.annotate(r'$\lambda_{V}$', xy=(self.ref_wavelengths['V'] + 10, 3.3),
                     xycoords='data', color='k',
                     horizontalalignment='left', size='large')

        # Set axes limits to clean the figures
        ax1.set_ylim(ymin=-0.9, ymax=3.8)
        ax1.set_xticklabels([])
        ax1.set_yticks(ax1.get_yticks()[1:-1])
        ax2.set_ylim(ymin=0, ymax=2)
        ax2.set_xticks(ax2.get_xticks()[1:-1])
        ax2.set_yticks(ax2.get_yticks()[1:-1])

        # Set the legends on the two axes
        ax1.legend(loc='best').draw_frame(False)
        ax2.legend(loc='best').draw_frame(False)

    def plot_cardelli_law_variability(self):
        """
        Plot the cardelli extinction law for several values of Rv.
        Expressed as Rlbd-Rv as a function of the inverse wavelength.
        """
        fig = self.pylab.figure(dpi=self.dpi)
        print(("Figure %i: CCM law variability (A(lbd)/Av)"
               % (fig.number)).center(80, '='))
        print("""
        Cardelli extinction law for several values of Rv. This figure express
        the extinction law variability as a function af Rv.\n""")

        ax = fig.add_axes([0.10, 0.09, 0.9, 0.89],
                          xlabel=r'$\lambda$ [$\AA$]',
                          ylabel=r'$A_{\lambda}/A_V$')

        ext = []
        colors = self.cmap(N.linspace(0, 1, len(self.Rvs)))
        for i, rv in enumerate(self.Rvs):
            extLaw = extinctionLaw(self.wavelength, Rv=rv)
            ax.plot(self.wavelength, extLaw, color=colors[i])
            ext.append(extLaw[0])
        scat = ax.scatter([self.wavelength[0]] * self.num,
                          ext, c=self.Rvs, cmap=(self.cmap),
                          visible=False, label='__nolgend__')
        cb = fig.colorbar(scat, format='%.2f')
        cb.set_label(r'$R_{V}$')
        ax.axvline(self.ref_wavelengths['V'], ls=':', color='k', lw=0.8)
        ax.axhline(1, ls=':', color='k', lw=0.8)

        ax.annotate(r'$\lambda_{V}$',
                    xy=(self.ref_wavelengths['V'] + 10, 3),
                    xycoords='data', color='k',
                    horizontalalignment='left', size='large')
        ax.annotate(r'$\frac{A_{\lambda}}{A_{V}} = '
                    r'a_{\lambda} + \frac{b_{\lambda}}{R_V}$',
                    xy=(6400, 1.85), xycoords='data',
                    color='k', horizontalalignment='left', size='xx-large')
        ax.set_ylim(ymin=0.05, ymax=4)
        ax.set_xlim(xmin=self.wavelength.min(), xmax=self.wavelength.max())

    def plot_rlbd_variability(self, Rv=3.1, ebmv=0.3):
        """
        Plot the cardelli extinction law for several values of Rv.
        The figure express the extinction law variability as a function af Rv.
        """
        fig = self.pylab.figure(dpi=self.dpi)
        print(("Figure %i: CCM law variability (R(lbd)-Rv)" %
               (fig.number)).center(80, '='))
        print("""
        CCM law expressed as R(lbd)-Rv for several values of Rv.\n""")

        # Create the figure and axe
        ax = fig.add_axes([0.10, 0.09, 0.9, 0.89],
                          xlabel=r'$1/\lambda$ ($\mu m^{-1}$)',
                          ylabel=r'$R_{\lambda}-R_V$')

        ext = []
        wlp = 1e4 / self.wavelength
        refV = 1e4 / self.ref_wavelengths['V']
        refB = 1e4 / self.ref_wavelengths['B']
        lbd = 1e4 / 9500

        colors = self.cmap(N.linspace(0, 1, len(self.Rvs)))
        for i, rv in enumerate(self.Rvs):
            extLaw = rv * (extinctionLaw(self.wavelength, Rv=rv) - 1)
            ax.plot(wlp, extLaw, color=colors[i])
            ext.append(extLaw[0])
        scat = ax.scatter([wlp[0]] * self.num, ext, c=self.Rvs,
                          cmap=(self.cmap), visible=False)
        ax.axvline(refV, ls=':', color='k', lw=0.8)
        ax.axvline(refB, ls=':', color='k', lw=0.8)
        ax.annotate(r'$\lambda_{B}$', xy=(refB, 3),
                    xycoords='data', color='k',
                    horizontalalignment='left', size='large')
        ax.annotate(r'$\lambda_{V}$', xy=(refV, 3),
                    xycoords='data', color='k',
                    horizontalalignment='left', size='large')
        ax.annotate(r'$R_{B}$', xy=(refB, -2.4),
                    xycoords='data', color='k',
                    horizontalalignment='right', size='large')
        ax.annotate(r'$R_{V}$', xy=(refV, -2.4),
                    xycoords='data', color='k',
                    horizontalalignment='right', size='large')
        cb = fig.colorbar(scat, format='%.2f')
        cb.set_label(r'$R_{V}$')
        ax.annotate(r'$R_{\lambda}=R_{V} \times a_{\lambda}+b_{\lambda}$',
                    xy=(2.5, 1), xycoords='data', color='k',
                    horizontalalignment='left', size='large')
        ax.set_ylim(ymin=-2.5, ymax=3.5)
        ax.set_xlim(xmin=wlp.min(), xmax=wlp.max())
        ax.set_xlim(xmax=3.29)

    def plot_albd_variability(self):
        """
        Plot the dust cloud transmission as a function of Rv and E(B-V).
        This figure express the variablitiy of the transmission as function
        of these two parameters. A degeneracy between the Rv and E(B-V)
        variabilities can be seen in the IR.
        """

        fig = self.pylab.figure(dpi=self.dpi)
        print(("Figure %i: Transmission variability" %
               (fig.number)).center(80, '='))
        print("""
        Dust cloud transmission variablitiy as a function of Rv and E(B-V).
        A degeneracy between the Rv and E(B-V) variabilities can be seen
        in the IR.\n""")

        ax = fig.add_axes([0.10, 0.09, 0.9, 0.89],
                          xlabel=r'$\lambda$ [$\AA$]',
                          ylabel=r'$T_{\lambda}$ ($=10^{-0.4 \times '
                          'R_V \times E(B-V) \times A_{\lambda}/A_V }$)')
        ext = []
        colors = self.cmap(N.linspace(0, 1, len(self.Rvs)))

        alphas = [1, 0.9, 0.8]
        letters = ['a', 'b', 'c']
        ls = [':', '--', '-']

        cst = 0
        for i, rv in enumerate(self.Rvs):
            for j, ebmv in enumerate(self.Ebmvs):
                extFact = extinctionFactor(self.wavelength, ebmv, Rv=rv)
                ax.plot(self.wavelength, extFact, ls=ls[j],
                        color=colors[i], alpha=alphas[j])
                if i == int(len(self.Rvs) / 2):
                    ex = extinctionFactor(3100, ebmv, Rv=rv)
                    ax.annotate(r'$(%s)$' % letters[j],
                                xy=(3100, ex), xycoords='data', color='k',
                                horizontalalignment='left',
                                verticalalignment='center', size='medium')
                    ax.annotate(r'$(%s)$ E(B-V)=%.2f' % (letters[j], ebmv),
                                xy=(7400, 0.15 - cst), xycoords='data',
                                color='k', size='medium')
                    cst += 0.05

            ext.append(extFact[0])

        scat = ax.scatter([self.wavelength[0]] * self.num, ext,
                          c=self.Rvs, cmap=(self.cmap), visible=False,
                          label='__nolgend__')
        cb = fig.colorbar(scat, format='%.2f')
        cb.set_label(r'$R_{V}$')
        ax.axhline(1, ls=':', color='k', lw=0.8)

        ax.set_ylim(ymin=0, ymax=1.01)
        ax.set_xlim(xmin=self.wavelength.min(), xmax=self.wavelength.max())

# Alias
extinctionParameters = CCMextinctionParameters

if __name__ == "__main__":

    print("\nUse Extinction.py to plot the extinction laws.\n")

    EP = ExtinctionsPlots()

    # Plot all the avalaible extinction laws
    EP.plot_extinctionLaws()

    # Plot other useful figures
    EP.plot_cardelli_law()
    EP.plot_cardelli_law_variability()
    EP.plot_rlbd_variability()
    EP.plot_albd_variability()

    # See also EP.plot_all_figures()

    EP.pylab.show()
