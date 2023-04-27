#!/usr/bin/env python3

import numpy as np;
import scipy as sp;
from scipy import optimize;
from scipy import interpolate;
import optparse as parser;
import matplotlib.pyplot as plt;

kb = 8.617e-5;

def nelecfunc():
    eigenval = 'EIGENVAL';
    fh = open(eigenval);
    
    [fh.readline() for i in range(5)];
    return int(fh.readline().split()[0]);
def fermifunc(elist, mu, temp):
    f = np.array([]);
    if temp == 0:
        beta = 0;
    else:
        beta = 1/kb/temp;
    for i in range(len(elist)):
        if ( beta == 0 ):
            if ( elist[i] < mu ):
                f = np.append(f, 1);
            else:
                f = np.append(f, 0);
        else:
            if (beta*(elist[i]-mu) > 20):
                f = np.append(f, 0);
            else:
                f = np.append(f, 1/(np.exp(beta*(elist[i]-mu))+1));
    return f;
def numdiff(mu, elist, ndos, temp, nelec):
    f = fermifunc(elist, mu, temp);
    #print(np.trapz(ndos*f, elist));
    return(np.trapz(ndos*f, elist)-nelec);
def get_mu(elist, ndos, temp, nelec, mu0):
    mu = sp.optimize.root(numdiff, mu0, args=(elist,ndos,temp,nelec), tol=1e-5);
    return mu.x;
# def get_mu(elist, ndos, temp, nelec, emin, emax):
#     print(nelec);
#     left=emin; right = emax; mid = (left + right)/2;
#     while (right-left > 1e-8):
#         mid = (left + right)/2;
#         nleft = numdiff(left, elist, ndos, temp, nelec);
#         nright = numdiff(right, elist, ndos, temp, nelec);
#         nmid = numdiff(mid, elist, ndos, temp, nelec);
#         #print(left, right, nleft, nright);
#         if ( nmid * nright < 0 ): left = mid;
#         else: right= mid;
        
#     return mid;

parser = parser.OptionParser();
parser.add_option("--doscar", dest='doscar', type='str', default='DOSCAR', help='');
parser.add_option("--T0", dest='T0', type='int', default=100, help='');
parser.add_option('--T1', dest='T1', type='int', default=2000, help='');
parser.add_option('--dT', dest='dT', type='int', default=10, help='');
parser.add_option('--mu', dest='mu', type='float', help='');
(options, args) = parser.parse_args();

T0 = options.T0;
T1 = options.T1;
dT = options.dT;
doscar = options.doscar;

fh = open(doscar);
natom = int(fh.readline().split()[1]);

[fh.readline() for i in range(4)];
line = np.array([float(i) for i in fh.readline().split()]);
emax = line[0]; emin=line[1];
nline = int(line[2]);
efermi = line[3];
nelec = nelecfunc();

dos = np.array([]);
for i in range(nline):
    line = np.array([float(i) for i in fh.readline().split()]);
    dos = np.append(dos, line);
dos = np.reshape(dos, (nline, -1));
(l, w) = np.shape(dos);
x = np.linspace(max(emin, efermi-0.2), min(emax, efermi+0.2), 1000);
x = np.append(x, dos[:,0]);
x = np.sort(np.unique(x), axis=None);
l = len(x);
dosnew = x;
for i in range(w-1):
    f=interpolate.interp1d(dos[:,0],dos[:,i+1], kind='cubic');
    dosnew = np.append(dosnew, f(x));
dos = np.transpose(np.reshape(dosnew, (w, l)));
#print(dos);
elist = dos[:,0];

#print(l, w);

if ( w == 3 ):
    ndos = dos[:,1];
if ( w == 5 ):
    ndos = dos[:,1] + dos[:,2];

mu0 = get_mu(elist, ndos, 0, nelec, efermi);
f0 = fermifunc(elist, mu0, 0);
U0 = np.trapz(ndos*elist*f0, elist)/natom;

temp = T0;
while temp <= T1:
    mu = get_mu(elist, ndos, temp, nelec, efermi);
    f = fermifunc(elist, mu, temp);
    e = np.trapz(ndos*elist*(f-f0), elist)/natom;
    vs = np.array([]);
    for i in f:
        if (( i==0 ) or( i==1)): vs = np.append(vs, 0);
        else: vs = np.append(vs, -i*np.log(i)-(1-i)*np.log(1-i))
    s = kb*temp*np.trapz(ndos*vs, elist)/natom;
    #plt.figure("vs");
    #plt.plot(elist, ndos* vs);
    print(temp, e-s, e, s);

    temp += dT;

# if ( w == 3):
#     plt.figure("DOS");
#     #print(mu);
#     #print(numdiff(mu, elist, ndos, temp, nelec));
#     mu = get_mu(elist, ndos, 0, nelec, efermi);
#     f0 = fermifunc(elist, mu, 0);
#     mu = get_mu(elist, ndos, 100, nelec, efermi);
#     f = fermifunc(elist, mu, 100);
#     plt.plot(elist, ndos, '.-');
#     plt.plot(elist, ndos*f0, 'b-');
#     plt.plot(elist, ndos*f, '-k');
#     plt.show();
