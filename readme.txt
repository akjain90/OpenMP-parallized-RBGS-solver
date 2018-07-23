1) The program is callable by ./rbgs nx ny c 

2) Number of threads is set with 

export OMP_NUM_THREADS = 

3) To run the code please enter the following :

likwid-pin -c N:0-31 ./rbgs nx ny c (Thread pinning)
