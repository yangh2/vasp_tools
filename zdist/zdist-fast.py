#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
from scipy.stats import multivariate_normal
import matplotlib.pyplot  as plt
import subprocess

parser=OptionParser();
parser.add_option('-s', "--smear", dest="smear", default=0.2, type='float',
                  help="Gaussion smearing [default: %default A]");
parser.add_option('-c', "--nbins", dest="bin_num", default=300, type='int',
                  help="number of bins in histogram [default: %default].");
parser.add_option('-p', "--parallel", dest='pvec', default='0,0,1', type='str',
                  help='Projection vector [default: %default]');
parser.add_option("-z", "--depth", dest='depth', default="1,-1", type='str',
                  help='Specify thickness [default: %default]');
parser.add_option("-n", '--samples', dest='nsamples', default=1, type='int',
                  help='number of samples [default: %default]');
(options, args)=parser.parse_args();

fxyz=args[0];

bin_num=options.bin_num;
smear=options.smear;
pvec=np.array([float(i) for i in options.pvec.split(',')]);
depth=np.array([float(i) for i in options.depth.split(',')]);
isdepth=depth[0]<depth[1];
nsamples=options.nsamples;

bin_sz=(depth[1]-depth[0])/bin_num;

fid = open(fxyz, 'r');
lino=0;
nl=1000;
ncount=0;
atom_species=np.array([]);
while True:
    line=fid.readline();
    lino=lino+1;
    if (lino%nl==1):
        continue;
    if (lino%nl==2):
        continue;
    if (lino%nl==3):
        continue;
    if (lino%nl==4):
        na=int(line.split()[0]);

    for atom_id in range(na):
        line=fid.readline();
        lino=lino+1;
        atom_type=int(line.split()[3]);
        if (atom_type not in atom_species):
            atom_species=np.append(atom_species, atom_type);
    break;
fid.close()
atom_species = np.sort(atom_species);

with open('zdist.in', 'w') as fid:
    fid.write("%d\n" % len(atom_species));
    for type_id in range(len(atom_species)):
        fid.write("%d " % atom_species[type_id]);
    fid.write("\n")
    fid.write("%d\n" % bin_num);
    fid.write("%.6f\n" % bin_sz);
    fid.write("%.6f\n" % smear);
    [fid.write("%.6f " % i) for i in pvec];
    fid.write("\n")
    fid.write("{0:.6f} {1:.6f}\n".format(depth[0],depth[1]));
    fid.write("%d\n" % nsamples);
    fid.write(fxyz);
    fid.write("\n")
fid.close()

print('Processing ...');
cmd='zdist';
subprocess.call(cmd, shell=True);

