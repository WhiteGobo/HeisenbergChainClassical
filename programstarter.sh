#$ -S /bin/sh

# options:
#$ -q normal@node006
#$ -wd /home/student/r/rfechner/heisenbergchain
#$ -N rfclassicheischain1


# run:
/home/student/r/rfechner/heisenbergchain/program2 < input1.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input2.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input3.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input4.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input5.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input6.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input7.dat&
/home/student/r/rfechner/heisenbergchain/program2 < input8.dat&

wait

