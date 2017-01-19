.. image:: http://readthedocs.org/projects/extinctions/badge/?version=latest
   :target: http://extinctions.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
      
.. image:: https://landscape.io/github/nicolaschotard/extinctions/master/landscape.svg?style=flat
   :target: https://landscape.io/github/nicolaschotard/extinctions/master
   :alt: Code Health

.. image:: https://travis-ci.org/nicolaschotard/extinctions.svg?branch=master
   :target: https://travis-ci.org/nicolaschotard/extinctions
   :alt: Travis CI build status (Linux)

.. image:: https://codecov.io/gh/nicolaschotard/extinctions/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/nicolaschotard/extinctions	 

.. inclusion-marker-do-not-remove
	 
extinctions
-----------

Python package including different extinction laws and dust maps. Useful to

- get E(B-V) for a set of coordinates (or list of coordinates) from different
  
  - distant sources (IRSA, NED)
  - local maps (SFD98, Schlafly 2014, Planck 2013, Green 2015)

- compute the ISM transmission for different extincion laws
  
  - CCM89: Cardelli, Clayton and Mathis (`<http://adsabs.harvard.edu/abs/1989ApJ...345..245C>`_)
  - OD94: O'Donnell (`<http://adsabs.harvard.edu/abs/1994ApJ...422..158O>`_)
  - FM98: Fitzpatrick & Massa (1998)
  - G08: Goobar (`<http://adsabs.harvard.edu/abs/2008ApJ...686L.103G>`_)
    
Installation
------------

To install::

  git clone https://github.com/nicolaschotard/extinctions.git
  pip install extinctions/

To install in a local directory ``mypath``, use::

  pip install --prefix='mypath' extinctions/

and do not forget to add it to your PYTHONPATH.

To upgrade to a new version (after a ``git pull`` or a local modification), use::

  pip install --upgrade (--prefix='mypath') extinctions/

To install a release version (no release version available yet)::

  pip install http://github.com/nicolaschotard/extinctions/archive/v0.1.tar.gz

Also works with the master::

  pip install (--upgrade) https://github.com/nicolaschotard/extinctions/archive/master.zip

In the future, release versions will be listed at this `location
<http://github.com/nicolaschotard/extinctions/releases>`_.


Dependencies
------------

`extinctions` has a few python dependencies listed in the `requirements
<requirements.txt>`_ file. To install them all, use::

  pip install -r requirements.txt

Requirements are:

- numpy==1.11.2
- scipy==0.17.1
- matplotlib==1.5.1
- seaborn==0.7.0
- pyyaml==3.11
- astropy==1.2.1
- requests==2.10.0
- astroquery==0.3.1
- wget==2.2
- healpy==1.9.1


  
Dust map setup
--------------

You can automatically download the dust maps listed in the `maps.yaml
<extinction/data/maps.yaml>`_ file using the script `get_maps.py`. It
will, by default, put them in $HOME/.extinction/maps, but other
locations are of course possible (option `--outdir`). If you choose to
store them in an other diretory, you must set a $MAPSDIR environment
variable pointing to it. Already existing maps in the output directory
will not be downloaded again.

The available dust maps are:

- `SFD98 <http://lambda.gsfc.nasa.gov/product/foreground/dust_map.cfm>`_, full sky Healpy format
- SFD98 `north <http://www.sdss3.org/svn/repo/catalogs/dust/trunk/maps/SFD_dust_4096_ngp.fits>`_ and `south <http://www.sdss3.org/svn/repo/catalogs/dust/trunk/maps/SFD_dust_4096_sgp.fits>`_ dust maps
- `Planck <http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/previews/HFI_CompMap_ThermalDustModel_2048_R1.20/index.html>`_
- `Schlafly 2014 <http://lambda.gsfc.nasa.gov/product/foreground/fg_ebv_map_info.cfm>`_
- `Green 2015 <http://lambda.gsfc.nasa.gov/product/foreground/fg_ebv_2015_map_info.cfm>`_


Download the maps
-----------------

To download the extinction maps listed above, use the script
`get_maps.py`. The following options are available::

  get_maps.py -h
  usage: get_maps.py [-h] [--outdir OUTDIR] [--update] [--list]
                   [--select SELECT] [--exclude EXCLUDE]

  optional arguments:
     -h, --help         show this help message and exit
     --outdir OUTDIR    Output directory in where to put the dust maps
     --update           Update the maps directory in case of changes of maps.yaml
     --list             List of available maps and exit
     --select SELECT    Select maps to download (coma separated)
     --exclude EXCLUDE  Exclude map(s) (coma separated).If the select option is
                        used, the exclude option will be ignored.


Usage
-----

Get E(B-V)
``````````

Here is an example of how to get the value of E(B-V) for a set of
coordinates (RA,DEC):

.. code:: python
	  
  In [1]: ra, dec = 340.83, -9.59 # could also be arrays or list of ra/dec
  In [2]: from extinctions import reddening
  In [3]: red = reddening.Reddening(ra, dec) # ra dec can also be lists of coordinates
  INFO: Loading the maps from local directory /home/chotard/.extinction/maps/
  - green is loaded
  - schlafly is loaded
  - sfd is loaded
  - planck is loaded

You can then get E(B-V) from different sources:

.. code:: python
	  
  # from the local maps
  In [4]: red.query_local_map(dustmap='sfd')
  Out[4]: 0.047723956233310674
  In [5]: red.query_local_map(dustmap='schlafly')
  Out[5]: 0.062566755984547445
  In [6]: red.query_local_map(dustmap='planck')
  Out[6]: 0.052058338535565914
  In [7]: red.query_local_map(dustmap='green')
  Out[7]: 0.035786284169101318


  # from the SFD98 north/south maps using `sncosmo`
  In [8]: red.from_sncosmo()
  Out[8]: array([ 0.0473752])

  # Using astroquery
  In [9]: red.from_astroquery()
  Downloading http://irsa.ipac.caltech.edu//workspace/TMP_XG1Joz_30445/DUST/340.8300_-9.5900.v0001/extinction.tbl
  |==============================================================================================| 4.3k/4.3k (100.00%)         0s
  Out[9]: [0.047377326565143825]


Make some plots
```````````````

To have a look at the different extinction laws amd dust maps, you can
use the script `extinction_plots.py`.
