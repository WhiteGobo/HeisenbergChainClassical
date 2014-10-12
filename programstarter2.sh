#$ -S /bin/sh

# options:
#$ -q normal@node011
#$ -wd /home/student/r/rfechner/heisenbergchain
#$ -N rfclassicheischain1


# run:
/home/student/r/rfechner/heisenbergchain/program2 < input9.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input10.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input11.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input12.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input13.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input14.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input15.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input16.dat&

wait

