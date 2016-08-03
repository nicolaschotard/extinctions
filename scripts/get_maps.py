#!/usr/bin/env python

import os
import yaml
import wget
from argparse import ArgumentParser

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
                        help="Update the maps.yaml file")
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

    m = "https://raw.githubusercontent.com/nicolaschotard/Extinction/master/maps/maps.yaml"
    print "INFO: downloading the map list from the github repository\n"
    download(m, outdir, update=args.update)
    maps = yaml.load(open(outdir+"/maps.yaml"))
    print "\nINFO: Downloading the following maps:"
    for i, m in enumerate(sorted(maps)):
        print "       %i: %s" % (i, m)
    print
    for m in sorted(maps):
        download(maps[m], outdir)
