.. image:: https://landscape.io/github/nicolaschotard/Extinction/master/landscape.svg?style=flat
   :target: https://landscape.io/github/nicolaschotard/Extinction/master
   :alt: Code Health

.. image:: https://travis-ci.org/nicolaschotard/Extinction.svg?branch=master
   :target: https://travis-ci.org/nicolaschotard/Extinction
   :alt: Travis CI build status (Linux)

Extinction
==========

**Warning**: Package under development

Python package including different extinction laws and dust maps.

Installation
------------

To install::

  git clone https://github.com/nicolaschotard/Extinction.git
  pip install Extinction/

To install in a local directory `mypath`, use::

  pip install --prefix='mypath' Extinction/

and do not forget to add it to your PYTHONPATH.

To upgrade to a new version (after a `git pull` or a local modification), use::

  pip install --upgrade (--prefix='mypath') Extinction/


In the future, release versions will be listed
`here <http://github.com/nicolaschotard/Extinction/releases>`_, and
installed using, e.g.::

  pip install http://github.com/nicolaschotard/Extinction/archive/v0.1.tar.gz


Dust map setup
--------------

The first thing needed to use this package is to download the dust
maps stored in the `maps.yaml <extinction/data/maps.yaml>`_ file. The
script `get_maps.py` automatically downloads these maps in (by
default) $HOME/.extinction/maps. Other locations are of course
possible with the option `--outdir`, as long at this output directory
is correctly added to the PATH environment variable. Already existing
maps in the output directory will not be downloaded again.

Available dust maps are for now:

- `SFD98 <http://lambda.gsfc.nasa.gov/product/foreground/dust_map.cfm>`_, full sky Healpy format
- SFD98 `north <http://www.sdss3.org/svn/repo/catalogs/dust/trunk/maps/SFD_dust_4096_ngp.fits>`_ and `south <http://www.sdss3.org/svn/repo/catalogs/dust/trunk/maps/SFD_dust_4096_sgp.fits>`_ dust maps
- `Planck <http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/previews/HFI_CompMap_ThermalDustModel_2048_R1.20/index.html>`_
- `Schlafly 2014 <http://lambda.gsfc.nasa.gov/product/foreground/fg_ebv_map_info.cfm>`_
- `Green 2015 <http://lambda.gsfc.nasa.gov/product/foreground/fg_ebv_2015_map_info.cfm>`_


Dependencies
------------

`Extinction` has a few python dependencies:

- numpy
- matplotlib
- seaborn
- astropy / astroquery  
- healpy

Usage
-----

You need to get the maps listed above. To do so, use the script `get_maps.py`::

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

