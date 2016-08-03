# Extinction

***Warning***: Package under development

Python package including different extinction laws and dust maps.

## Installation

To install

```
git clone https://github.com/nicolaschotard/Extinction.git
pip install Extinction/
```

To install in a local directory `mypath`, use
```
pip install --prefix='mypath' Extinction/
```
and do not forget to add it to your PYTHONPATH.

To upgrade to a new version (after a `git pull` or a local modification), use
```
pip install --upgrade (--prefix='mypath') Extinction/
```

In the future, release versions will be listed
[here](http://github.com/nicolaschotard/Extinction/releases), and
installed using, e.g.:

```
pip install http://github.com/nicolaschotard/Extinction/archive/v0.1.tar.gz
```

## Dust map setup

The first thing to do would be to download all the needed dust maps:

- [E(B-V) (SFD)](http://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits) from http://lambda.gsfc.nasa.gov/product/foreground/fg_sfd_get.cfm
- [E(B-V) (Schlafly et al 2014)](http://lambda.gsfc.nasa.gov/data/foregrounds/EBV/ps1-ebv-4.5kpc.fits) from http://lambda.gsfc.nasa.gov/product/foreground/fg_ebv_map_get.cfm
- [E(B-V) (Green et al. 2015)](http://lambda.gsfc.nasa.gov/data/foregrounds/EBV/lambda_green_dust_map_2d.fits) from http://lambda.gsfc.nasa.gov/product/foreground/fg_ebv_2015_map_get.cfm
- [Planck astrophysical foregrounds from parametric component separation, 2015 version](http://pla.esac.esa.int/pla/aio/product-action?MAP.MAP_ID=COM_CompMap_ThermalDust-commander_2048_R2.00.fits), from [here](https://wiki.cosmos.esa.int/planckpla2015/index.php/CMB_and_astrophysical_component_maps#Thermal_dust_emission_2) (Question for Celine: is it the right one?)
- other known maps?

I'm working on a way to automatically download these map, and put them in (by default) $HOME/.extinction/maps. Other locations would of course be possible, as long at the PATH is correctly updated.

## Dependencies

`Extinction` has a few python dependencies:

- numpy
- matplotlib
- seaborn
- astropy / astroquery
- healpy

## Usage

TBD