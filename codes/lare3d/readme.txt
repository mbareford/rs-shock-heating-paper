This folder contains the source code for the variant of LARE3D (based on v2.10) used for this paper,
"Shock Heating in Numerical Simulations of Kink-unstable Coronal Loops".

There are also environment scripts and makefiles for two computer clusters, wardlaw, the UKMHD cluster
at St Andrews University, and dirac, a supercomputer hosted at Durham university. The ssh commands for
these two facilities are given below.

ssh -Y <username>@wardlaw10.mcs.st-and.ac.uk
ssh -Y <username>@login.cosma.dur.ac.uk

Once, you have decided which facility to use, remove the suffix from the env.sh and Makefile files, and,
after setting up suitable directory hierarchies (see Makefile and replace the terms in angle brackets),
do the following.

source env.sh
make JBN=02

The latter command will build LARE3D with ./src/control/control_02.f90, which is the medium resolution run.
Low resolution is JBN=01 and high is JBN=03.

To submit the job, consult the relevant documentation for the facility on which you intend to run the simulation.
Note, you must create the data directory (see the data_dir variable in the control source code file used to
compile the exe) before you run the job.

