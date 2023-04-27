#!/usr/bin/env python3

import numpy as np

fid = open("chopbasis.in", "r");
a1=np.array([float(i) for i in fid.readline().split()])
a2=np.array([float(i) for i in fid.readline().split()])
a3=np.array([float(i) for i in fid.readline().split()])
vol = np.dot(a1,np.cross(a2, a3));
b1 = np.cross(a2, a3)/vol;
b2 = np.cross(a3, a1)/vol;
b3 = np.cross(a1, a2)/vol;
rbasis=np.array([b1,b2,b3]);

fid =open("chopin.xyz","r");
veca=np.array([float(i) for i in fid.readline().split()])
vecb=np.array([float(i) for i in fid.readline().split()])
vecc=np.array([float(i) for i in fid.readline().split()])
basis=np.array([veca, vecb, vecc]);
print (basis);
natom=int(fid.readline().split()[0]);
print(natom);
inside=np.array([]);
outside=np.array([]);
for i in range(natom):
    line = fid.readline().split();
    dx = np.array([float(i) for i in line[0:3]]);
    t = int(line[3]);
    an = int(line[4]);
    rx = np.dot(np.transpose(basis),dx);
    #print(rx);
    dxnew= np.dot(np.transpose(rbasis), rx);
    #print(dxnew);
    while (dxnew[0] >0.5):
        dxnew[0]=dxnew[0]-1;
    while (dxnew[1] >0.5):
        dxnew[1]=dxnew[1]-1;
    while (dxnew[2] >0.5):
        dxnew[2]=dxnew[2]-1;
    while (dxnew[0] <-0.5):
        dxnew[0]=dxnew[0]+1;
    while (dxnew[1] <-0.5):
        dxnew[1]=dxnew[1]+1;
    while (dxnew[2] <-0.5):
        dxnew[2]=dxnew[2]+1;

    inside=np.append(inside,dxnew);
    inside=np.append(inside,t);
    inside=np.append(inside,an);
#print(inside);    
inside=np.ravel(inside);    
inside=np.reshape(inside,(-1,5));
win,l=np.shape(inside);
fid = open("chopempty.xyz", "w");
fid.write("{0} {1} {2}\n".format(a1[0],a1[1],a1[2]));
fid.write("{0} {1} {2}\n".format(a2[0],a2[1],a2[2]));
fid.write("{0} {1} {2}\n".format(a3[0],a3[1],a3[2]));
fid.write("{0}\n".format(win));
for i in range(win):
    fid.write("{0} {1} {2} {3:d} {4:d}\n".format(inside[i,0],inside[i,1],inside[i,2],int(inside[i,3]),int(inside[i,4])));
