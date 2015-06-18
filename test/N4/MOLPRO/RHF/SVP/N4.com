*** Molpro input file for N4
memory,20,m
basis=def2-SVP
geomtyp=xyz
geometry={
4                        ! number of atoms

 N,       0.5200,      0.5300,      0.5100
 N,      -0.5200,     -0.5200,      0.5200
 N,      -0.5200,      0.5200,     -0.5200
 N,       0.5200,     -0.5200,     -0.5200
}
hf
put,molden,N4.molden;
