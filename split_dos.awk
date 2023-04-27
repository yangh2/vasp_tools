#!/usr/bin/awk -f

BEGIN{
    natom=0			# number of atom
    emax=0			
    emin=0
    enum=0			# number of bins
    efermi=0			# fermi energy
    dnum=0			# DOS split index
    nmod=0
    nshift=5
}

NR==1{
    natom=$1
}

{
    if (NR==6){
	emax=$1
	emin=$2
	enum=$3
	nmod=enum+1;
	efermi=$4
    }
    if ((NR>5)&&((NR-nshift)%nmod ==1)){
	emax=$1
	emin=$2
	enum=$3
	nmod=enum+1;
	efermi=$4
    }
    
    if ((NR>5)&&((NR-nshift)%nmod !=1)){
	$1=$1-efermi;
	print $0 >> "DOS_"dnum
	if (dnum>0){
	    print $1,$2,$3+$4+$5,$6+$7+$8+$9+$10 >> "DOS_spd_"dnum
	}
    }
    
    if ((NR>5)&&((NR-nshift)%nmod ==0)){
	dnum = dnum +1;
	print dnum;
    }
}

