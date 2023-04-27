#!/usr/bin/env python3

from optparse import OptionParser
import numpy as np

parser=OptionParser();
parser.add_option('-n', "--number", dest="num", default=1, type='int',
                  help="Number of atomic substitutions [default: %default A]");
parser.add_option('-f', "--from", dest="species_from", default=0, type='int',
                  help="Atomic number of old species [default: %default A]");
parser.add_option('-t', "--to", dest="species_to", default=0, type='int',
                  help="Atomic number of old species [default: %default A]");

(options, args)=parser.parse_args();
num=options.num
if ( num > 0 ):
    species_from=options.species_from
    species_to=options.species_to
else:
    num=-num;
    species_to=options.species_from
    species_from=options.species_to

fxyz=args[0];

species_list=np.loadtxt(fxyz, skiprows=4, usecols=3, dtype='int');
na=len(species_list);
num_species_from = na-(np.count_nonzero(species_list-species_from));
num_species_to = na-(np.count_nonzero(species_list-species_to));

if ( num_species_from < num ):
    print("Not enough atoms for substitution!!");
    exit(0);
    
shuffle=np.arange(num_species_from);
np.random.shuffle(shuffle);
selected=shuffle[0:num];

where_species_from=np.argwhere(species_list==species_from);
species_list[where_species_from[selected]]=species_to;

fid = open(fxyz, "r");
for lino in range(na+4):
    line = fid.readline().rstrip();
    if (lino <4 ):
        print(line);
        continue;

    words=line.split();
    words[3]=str(species_list[lino-4]);

    print(" ".join(words))
