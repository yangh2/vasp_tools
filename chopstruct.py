#!/usr/bin/env python3

import numpy as np
import optparse as parser;

parser = parser.OptionParser();
(options, args) = parser.parse_args();

fid=open("chopstruct.in","r")
cpos=np.array([float(i) for i in fid.readline().split()])
veca=np.array([float(i) for i in fid.readline().split()])
vecb=np.array([float(i) for i in fid.readline().split()])
vecc=np.array([float(i) for i in fid.readline().split()])
vol = np.dot(veca,np.cross(vecb, vecc));
b1 = np.cross(vecb, vecc)/vol;
b2 = np.cross(vecc, veca)/vol;
b3 = np.cross(veca, vecb)/vol;
print(b1)
print(b2)
print(b3);

fxyz=args[0]
fid =open(fxyz,"r");
line1=fid.readline();
line2=fid.readline();
line3=fid.readline();
natom=int(fid.readline().split()[0]);
print(natom);
inside=np.array([]);
outside=np.array([]);
db=np.zeros(3);
for i in range(natom):
    line = fid.readline().split();
    x = float(line[0]); y = float(line[1]); z = float(line[2]);
    t = int(line[3]);
    an = int(line[4]);
    dd = np.array([x-cpos[0],y-cpos[1],z-cpos[2]]);

    flag = 0;
    for ia in range(-1,2):
        for ib in range(-1,2):
            for ic in range(-1,2):
                db[0]=dd[0]+ia;
                db[1]=dd[1]+ib;
                db[2]=dd[2]+ic;
                a1=np.dot(db,b1);
                a2=np.dot(db,b2);
                a3=np.dot(db,b3);
                #print(a1,a2,a3);
                if ((a1<1)and(a1>=0)and(a2<1)and(a2>=0)and(a3<1)and(a3>=0)):
                    flag=1;
    if (flag == 1):
        inside=np.append(inside,np.array([x,y,z,t,an]));
    else:
        outside=np.append(outside,np.array([x,y,z,t,an]));

inside=np.reshape(inside,(-1,5));
outside=np.reshape(outside,(-1,5));

win,l=np.shape(inside);
wout,l=np.shape(outside);

print(win, wout);

fid = open("chopin.xyz", "w");
fid.write(line1);
fid.write(line2);
fid.write(line3);
fid.write("{0}\n".format(win));
for i in range(win):
    fid.write("{0} {1} {2} {3:d} {4:d}\n".format(inside[i,0],inside[i,1],inside[i,2],int(inside[i,3]),int(inside[i,4])));

fid = open("chopout.xyz", "w");
fid.write(line1);
fid.write(line2);
fid.write(line3);
fid.write("{0}\n".format(wout));
for i in range(wout):
    fid.write("{0} {1} {2} {3:d} {4:d}\n".format(outside[i,0],outside[i,1],outside[i,2],int(outside[i,3]),int(outside[i,4])));
