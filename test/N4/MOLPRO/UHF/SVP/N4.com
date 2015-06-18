*** Molpro input file for N4
memory,20,m
basis=def2-SVP
geomtyp=xyz
charge=1
geometry={
4                        ! number of atoms

 N,       0.5200,      0.5300,      0.5100
 N,      -0.5200,     -0.5200,      0.5200
 N,      -0.5200,      0.5200,     -0.5200
 N,       0.5200,     -0.5200,     -0.5200
}
uhf;save,mo

put,molden,N4_alpha.molden;orb,mo,set=1 ! alpha spins
put,molden,N4_beta.molden ;orb,mo,set=2 ! beta spins
