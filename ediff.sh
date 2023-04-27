#!/bin/bash

cwd=`pwd`;

echo -e "Al\tCo\tCu\tNEtot\tNE_gap-NE_fermi"
for i in $@;
do
    if ! [ -d $i ]; then continue; fi
    cd $i;
    if ! [ -f RunAt ]; then cd $cwd; continue; fi
    
    ans=`echo $i | tr '/' ' ' | awk '{print $1}' | tr '[A-Z][a-z]' ' '`
    echo -n "$ans "
    smear=`grep "ISMEAR" INCAR | tr '=' ' '| tr ';' ' '| awk '{print $2}'`
    
    if [ $smear -eq 1 ];
    then
	dos=`ls DOSCAR* -ltr | grep -v " -5" | tail -n1 | awk '{print $NF}'`;
    else
	dos=`ls DOSCAR* -ltr | grep -v " -5" | tail -n2 | head -n1| awk '{print $NF}'`;
    fi
    
    ediff.py $dos;
    cd $cwd;
done
