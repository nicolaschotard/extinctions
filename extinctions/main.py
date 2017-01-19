"""Main entry points for scripts."""


import os
import shutil
import wget
from argparse import ArgumentParser
from pkg_resources import resource_filename
import yaml

from extinctions import extinction


# Download the dust maps


def download(url, out, update=False):
    """Download a file in an output directory (force with update)."""
    local_path = out + "/" + os.path.basename(url)
    if os.path.exists(local_path):
        if update:
            os.remove(local_path)
        else:
            print " File", local_path, "already exist. Pass!"
            return
    print " - %s -> %s" % (url, out)
    wget.download(url, out=local_path)
    print ""


def list_maps(mdic):
    """List the available maps from a dictionnary of maps."""
    print "INFO: Available maps are"
    for j, mm in enumerate(sorted(mdic)):
        print "       %i: %s (%s) - %s" % (j, mm, mdic[mm]['size'], mdic[mm]['url'])


def get_maps(argv=None):
    """Download dust maps."""
    description = """Download dust maps."""
    prog = "get_maps.py"

    parser = ArgumentParser(prog=prog, description=description)
    parser.add_argument("--outdir",
                        help="Output directory in where to put the dust maps")
    parser.add_argument("--update", action='store_true', default=False,
                        help="Update the maps directory in case of changes of maps.yaml")
    parser.add_argument("--list", action='store_true', default=False,
                        help="List of available maps and exit")
    parser.add_argument("--select", help="Select maps to download (coma separated)")
    parser.add_argument("--exclude", help="Exclude map(s) (coma separated)."
                        "If the select option is used, the exclude option will be ignored.")
    args = parser.parse_args(argv)

    if args.outdir is not None:
        outdir = os.path.abspath(args.outdir)
        print "WARNING: Please make sure to add '%s' to your path in order for " % \
            outdir + "the code to find the maps. Use:"
        print " export DUSTMAP='%s'" % outdir
        print " or"
        print " setenv DUSTMAP '%s'" % outdir
    else:
        outdir = os.getenv('HOME') + "/.extinction/maps"

    if not os.path.isdir(outdir):
        print "INFO: Creating output directoy:", outdir
        os.makedirs(outdir)
    else:
        print "INFO: output directory used is", outdir

    # Get the maps yaml file
    m = resource_filename('extinctions', 'data/maps.yaml')
    print "INFO: getting the map list from the local github repository: %s" % m
    shutil.copy2(m, outdir)

    # Now get all or part of the maps
    maps = yaml.load(open(outdir + "/maps.yaml"))

    # List the maps
    if args.list:
        list_maps(maps)
        return

    # Selection by the user?
    if args.select is not None:
        select = args.select.split(',')
        for m in maps.keys():
            if m not in select:
                maps.pop(m)
    elif args.exclude is not None:
        exclude = args.exclude.split(',')
        for m in maps.keys():
            if m in exclude:
                maps.pop(m)

    # Check if some maps are already present
    for m in maps.keys():
        path = outdir + "/" + os.path.basename(maps[m]['url'])
        if os.path.exists(path):
            print " File", path, "already exist."
            maps.pop(m)

    print "INFO: Downloading the following maps:"
    for i, m in enumerate(sorted(maps)):
        print "       %i: %s" % (i, m)

    for m in sorted(maps):
        download(maps[m]['url'], outdir)


# Make plots


def extinction_plots(argv=None):
    """Plot the available extinction laws."""
    description = """Plot extinction curves."""
    prog = "extinction_plots.py"

    print "\nUse extinction.py to plot the extinction laws.\n"

    parser = ArgumentParser(prog=prog, description=description)
    parser.add_argument("--hide", action='store_true', default=False,
                        help="Do NOT show the figures")
    args = parser.parse_args(argv)

    eplot = extinction.ExtinctionsPlots()

    # Plot all the avalaible extinction laws
    eplot.plot_extinction_laws()

    # Plot other useful figures
    eplot.plot_cardelli_law()
    eplot.plot_cardelli_law_variability()
    eplot.plot_rlbd_variability()
    eplot.plot_albd_variability()

    # See also eplot.plot_all_figures()
    if not args.hide:
        extinction.pylab.show()
