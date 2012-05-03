hank: hankel_profiling.c
	gcc -lm -lblas -lgsl hankel_profiling.c -o hank

eval: eval_profiling.o
	g++ eval_profiling.o ../vergini/basis.o ../vergini/billiard.o ../vergini/colloc.o ../vergini/grid.o ../vergini/matrix.o ../vergini/nrutil.o ../vergini/useful.o -o eval -I/usr/include/atlas -L/usr/lib64/atlas -L/usr/lib64/lapack -L/usr/lib64/libblas -lgsl -llapack -lblas -lf77blas -latlas -lm

eval_profiling.o: eval_profiling.c
	g++ -c eval_profiling.c -I/usr/include/atlas -L/usr/lib64/atlas -lgsl -lgslcblas -llapack -lblas -lf77blas -latlas -lm

hank103.c: hank103.f
	f2c hank103.f