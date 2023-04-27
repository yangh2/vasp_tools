#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
from scipy.stats import multivariate_normal
import matplotlib.pyplot  as plt
import subprocess

parser=OptionParser();
parser.add_option('-s', "--smear", dest="smear", default=0.2, type='float',
                  help="Gaussion smearing [default: %default A]");
parser.add_option('-x', "--pixel-size", dest="px_sz", default=0.2, type='float',
                  help="Pixel size [default: %default A]");
parser.add_option('-c', "--canvas-size", dest="px_num", default=300, type='int',
                  help="Canvas size [default: %default px]. [-px_num*px_sz/2,px_num*px_sz/2]");
parser.add_option('-p', "--parallel", dest='pvec', default='0,0,1', type='str',
                  help='Projection vector [default: %default]');
parser.add_option("--define-x", dest='vec_x', default='1,0,0', type='str',
                  help='Define x vector on canvas [default: %default]');
parser.add_option("--define-y", dest='vec_y', default='0,1,0', type='str',
                  help='Define y vector on canvas [default: %default]');
parser.add_option("-z", "--depth", dest='depth', default="1,-1", type='str',
                  help='Specify thickness [default: %default]');
parser.add_option("-n", '--samples', dest='nsamples', default=1, type='int',
                  help='number of samples [default: %default]');
parser.add_option('-o', '--figure-only', action='store_true', dest='islabel', default='false',
                  help='Plot figures without labels');
(options, args)=parser.parse_args();

fxyz=args[0];

px_num=options.px_num;
px_sz=options.px_sz;
smear=options.smear;
pvec=np.array([float(i) for i in options.pvec.split(',')]);
vec_x=np.array([float(i) for i in options.vec_x.split(',')]);
vec_y=np.array([float(i) for i in options.vec_y.split(',')]);
depth=np.array([float(i) for i in options.depth.split(',')]);
isdepth=depth[0]<depth[1];
nsamples=options.nsamples;
islabel=options.islabel;

px_ct=int(px_num/2);
smear_mat=[[smear**2,0],[0,smear**2]];
smear_sz=int(smear/px_sz*3);
canvas=np.zeros((px_num,px_num));
width=px_num*px_sz;
height=width;

fid = open(fxyz, 'r');
lino=0;
nl=1000;
ncount=0;

with open('HAADF.in', 'w') as fid:
    fid.write("%d\n" % px_num);
    fid.write("%.6f\n" % px_sz);
    fid.write("%.6f\n" % smear);
    [fid.write("%.6f " % i) for i in vec_x];
    fid.write("\n")
    [fid.write("%.6f " % i) for i in vec_y];
    fid.write("\n")
    [fid.write("%.6f " % i) for i in pvec];
    fid.write("\n")
    fid.write("{0:.6f} {1:.6f}\n".format(depth[0],depth[1]));
    fid.write("%d\n" % nsamples);
    fid.write(fxyz);
    fid.write("\n")
fid.close()

print('Processing ...');
cmd='haadf';
subprocess.call(cmd, shell=True);

canvas=np.loadtxt('HAADF.out');
#print(np.shape(canvas));
print('\nPlotting ...');
# Plot 
x=np.arange(px_num);
y=np.arange(px_num);
fig = plt.figure("HAADF");
plt.pcolormesh(x,y,canvas,cmap="Greys");
plt.axis('square');
#plt.axis('tight');
xlocs=np.arange(0,px_num+1,step=int(px_num/10));
xlabels=[str(int(i*px_sz)) for i in xlocs];
plt.margins(x=0)
plt.xlim(0,px_num);
plt.ylim(0,px_num);
if not islabel:
    plt.xticks(xlocs,xlabels);
    plt.yticks(xlocs,xlabels);
    plt.xlabel('x [A]');
    plt.ylabel('y [A]');
else:
    plt.axis('off')
plt.savefig('HAADF.eps',bbox_inches = 'tight',
    pad_inches = 0);
plt.show();
