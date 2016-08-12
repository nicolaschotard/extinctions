#!/usr/bin/env python

import os
import sys
import yaml
import wget
from argparse import ArgumentParser
from pkg_resources import resource_filename
from shutil import copy2

def download(f, out, update=False):
    path = out+"/"+os.path.basename(f)
    if os.path.exists(path):
        if update:
            os.remove(path)
        else:
            print " File", path, "already exist. Pass!"
            return
    print " - %s -> %s" % (f, out)
    wget.download(f, out=path)
    print ""

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--outdir",
                        help="Output directory in where to put the dust maps")
    parser.add_argument("--update", action='store_true', default=False,
                        help="Update the maps directory in case of changes of maps.yaml")
    parser.add_argument("--list", action='store_true', default=False,
                        help="List of available maps and exit")
    parser.add_argument("--select", help="Select maps to download (coma separated)")
    parser.add_argument("--exclude", help="Exclude map(s) (coma separated)."
                        "If the select option is used, the exclude option will be ignored.")
    args = parser.parse_args()

    if args.outdir is not None:
        outdir = args.outdir
        print "Please make sure to add '%s' to your path in order for the code to find the maps using" % outdir
        print "export DUSTMAP='%s'" % outdir
        print " or"
        print "setenv DUSTMAP '%s'" % outdir
    else:
        outdir = os.getenv('HOME') + "/.extinction/maps"
    
    if not os.path.isdir(outdir):
        print "INFO: Creating output directoy:", outdir
        os.makedirs(outdir)
    else:
        print "INFO: output directory used is", outdir

    # Get the maps yaml file
    m = resource_filename('Extinction', 'data/maps.yaml')
    print "INFO: getting the map list from the local github repository: %s" % m
    copy2(m, outdir)

    # Now get all or part of the maps
    maps = yaml.load(open(outdir+"/maps.yaml"))

    # List the maps
    if args.list:
        print "INFO: Available maps are"
        for i, m in enumerate(sorted(maps)):
            print "       %i: %s (%s) - %s" % (i, m, maps[m]['size'], maps[m]['url'])
        sys.exit()
        
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
        path = outdir+"/"+os.path.basename(maps[m]['url'])
        if os.path.exists(path):
            print " File", path, "already exist."
            maps.pop(m)
    
    print "INFO: Downloading the following maps:"
    for i, m in enumerate(sorted(maps)):
        print "       %i: %s" % (i, m)
    print
    for m in sorted(maps):
        download(maps[m]['url'], outdir)
