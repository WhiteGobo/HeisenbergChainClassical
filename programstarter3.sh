#$ -S /bin/sh

# options:
#$ -q normal@node012
#$ -wd /home/student/r/rfechner/heisenbergchain
#$ -N rfclassicheischain1


# run:
/home/student/r/rfechner/heisenbergchain/program2 < input17.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input18.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input19.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input20.dat&

wait

