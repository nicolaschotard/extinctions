#!/usr/bin/env python

import os
import yaml
import wget
from argparse import ArgumentParser

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--outdir",
                        help="Output directory in where to put the dust maps")
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

    print "INFO: downloading the map list from the github repository"
    
    maps = yaml.load(open('maps.yaml'))
