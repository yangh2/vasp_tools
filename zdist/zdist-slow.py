#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
from scipy.stats import multivariate_normal
import matplotlib.pyplot  as plt

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
parser.add_option('-o', '--figure-only', action='store_true', dest='islabel', default='false',
                  help='Plot figures without labels');
(options, args)=parser.parse_args();

fxyz=args[0];

bin_num=options.bin_num;
smear=options.smear;
pvec=np.array([float(i) for i in options.pvec.split(',')]);
depth=np.array([float(i) for i in options.depth.split(',')]);
isdepth=depth[0]<depth[1];
nsamples=options.nsamples;
islabel=options.islabel;

bin_size=(depth[1]-depth[0])/bin_num;
bin_lmt=int(smear/bin_size*3);
smear=smear**2;
canvas=np.zeros(bin_num);

fid = open(fxyz, 'r');
lino=0;
nl=1000;
ncount=0;
print("Processing");
while True:
    line=fid.readline();
    lino=lino+1;
    if (lino%nl==1):
        a=np.array([float(i) for i in line.split()]);
        continue;
    if (lino%nl==2):
        b=np.array([float(i) for i in line.split()]);
        continue;
    if (lino%nl==3):
        c=np.array([float(i) for i in line.split()]);
        continue;

    if (lino%nl==4):
        ncount+=1;
        print("{0:4d} out of {1:4d}".format(ncount, nsamples), end='\r');
        lat_cont=np.array([a,b,c]);
        na=int(line.split()[0]);
        nl=na+4;
        
        for atom_id in range(na):
            line=fid.readline();
            lino=lino+1;
            atoms_pos=np.array([float(i) for i in line.split()[0:3]]);
            atom_type=int(line.split()[3]);
            atom_pos=np.copy(atoms_pos[:]);

            atom_rpos=np.dot(atom_pos, lat_cont); # real space position
            pos_z = np.dot(atom_rpos, pvec);      # projected position: z

            if (( isdepth ) and ((pos_z < depth[0]) or (pos_z > depth[1]))):
                continue;

            px_nz = int((pos_z-depth[0])/bin_size);
            var = multivariate_normal(mean=pos_z, cov=smear)
            #print(px_nx,px_ny);
            for paint_z in range(-bin_lmt,bin_lmt+1):
                pos_pz=px_nz+paint_z;
                if (( pos_pz<0 ) or (pos_pz >=bin_num)):
                    continue;
                canvas[pos_pz]+=var.pdf((px_nz+paint_z)*bin_size+depth[0]);
        if ( ncount >= nsamples):
            break;

fid.close();

canvas=canvas/np.sum(canvas)/bin_size;
z=np.linspace(depth[0], depth[1], bin_num);
np.savetxt('zdist.dat', np.transpose(np.array([z,canvas])));
