#!/usr/bin/env python3

import numpy as np
from optparse import OptionParser
from scipy.stats import multivariate_normal
import matplotlib.pyplot  as plt

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
parser.add_option("-k", "--decay-ratio", dest='dratio', default="1.0", type='float',
                  help='Specify exponential decay ratio with distance [default: %default 1/A]');
parser.add_option("-n", '--samples', dest='nsamples', default=1, type='int',
                  help='number of samples [default: %default]');
parser.add_option('-o', '--figure-only', action='store_true', dest='islabel', default='false',
                  help='Plot figures without labels');
(options, args)=parser.parse_args();

fxyz=args[0];

px_num=options.px_num;
px_sz=options.px_sz;
smear=options.smear;
dratio=options.dratio;
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
        
        # atoms_pos=np.array([]);
        # atoms_type=np.array([]);
        # for atom_id in range(na):
        #     line=fid.readline();
        #     lino=lino+1;

        #     atoms_pos=np.append(atoms_pos, np.array([float(i) for i in line.split()[0:3]]));
        #     atoms_type=np.append(atoms_type, int(line.split()[3]));

        # atoms_pos=np.reshape(atoms_pos, (-1,3));
        #atoms_rpos=np.dot(atoms_pos,lat_cont);
        # print(na);
        # print(lat_cont);
        # print(atoms_pos);
        # print(atoms_rpos);

        for atom_id in range(na):
            line=fid.readline();
            lino=lino+1;
            atoms_pos=np.array([float(i) for i in line.split()[0:3]]);
            atom_type=int(line.split()[3]);
            atom_pos=np.copy(atoms_pos[:]);

            for pb_x in range(-1,2,1):
                for pb_y in range(-1,2,1):
                    atom_pos[0]=atoms_pos[0]+pb_x;
                    atom_pos[1]=atoms_pos[1]+pb_y;
                    #print(atoms_pos[atom_id,:], atom_pos);
                    atom_rpos=np.dot(atom_pos, lat_cont); # real space position
                    pos_x = np.dot(atom_rpos, vec_x);     # projected position: x
                    pos_y = np.dot(atom_rpos, vec_y);     # projected position: y
                    pos_z = np.dot(atom_rpos, pvec);      # projected position: z

                    if (( isdepth ) and ((pos_z < depth[0]) or (pos_z > depth[1]))):
                        continue;
                    px_nx=int(pos_x/px_sz);
                    px_ny=int(pos_y/px_sz);
                    var = multivariate_normal(mean=[atom_rpos[0],atom_rpos[1]], cov=smear_mat)
                    #print(px_nx,px_ny);
                    for paint_x in range(-smear_sz,smear_sz+1):
                        for paint_y in range(-smear_sz,smear_sz+1):
                            pos_px=px_nx+paint_x+px_ct;
                            pos_py=px_ny+paint_y+px_ct;
                            if (( pos_px<0 ) or (pos_px >=px_num)):
                                continue;
                            if (( pos_py<0 ) or (pos_py >=px_num)):
                                continue
                            canvas[pos_px, pos_py]+=np.exp(-dratio*pos_z)*var.pdf([(px_nx+paint_x)*px_sz, (px_ny+paint_y)*px_sz]);
        if ( ncount >= nsamples):
            break;

fid.close();
canvas=canvas/np.max(canvas);
canvas=-canvas;
np.savetxt('TEM_canvas.dat', canvas);
print('Plotting ...');
# Plot 
x=np.arange(px_num);
y=np.arange(px_num);
fig = plt.figure("TEM");
plt.pcolormesh(x,y,canvas,cmap="Greys");
plt.axis('square');
xlocs=np.arange(0,px_num+1,step=int(px_num/10));
xlabels=[str(int(i*px_sz)) for i in xlocs];
if not islabel:
    plt.xticks(xlocs,xlabels);
    plt.yticks(xlocs,xlabels);
    plt.xlabel('x [A]');
    plt.ylabel('y [A]');
    plt.xlim(0,px_num);
    plt.ylim(0,px_num);
else:
    plt.axis('off')
plt.savefig('TEM.eps');    
plt.show();
