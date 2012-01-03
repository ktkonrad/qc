CC = g++
CFLAGS = -g -lgsl -lgslcblas -lm

objects = count_nodal_domains.o random_percolation.o stack2.o util.o
dependencies = ../vergini/billiard.o ../vergini/useful.o ../vergini/nrutil.o util_verg.o

count: nodal_domain_driver.o $(objects) $(dependencies)
	$(CC) $(CFLAGS) -o count nodal_domain_driver.o $(objects) $(dependencies)

perc: percolation_driver.o $(objects)
	$(CC) $(CFLAGS) -o perc percolation_driver.o $(objects)

nodal_domain_driver.o: nodal_domain_driver.c
	$(CC) $(CFLAGS) -c nodal_domain_driver.c

percolation_driver.o: percolation_driver.c
	$(CC) $(CFLAGS) -c percolation_driver.c

count_nodal_domains.o: count_nodal_domains.c count_nodal_domains.h
	$(CC) $(CFLAGS) -c count_nodal_domains.c

random_percolation.o: random_percolation.c random_percolation.h
	$(CC) $(CFLAGS) -c random_percolation.c

stack2.o: util/stack2.c util/stack2.h
	$(CC) $(CFLAGS) -c util/stack2.c

util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c

util_verg.o: util_verg.c util_verg.h util.o
	$(CC) $(CFLAGS) -c util_verg.c

.PHONY : clean
clean:
	-rm count $(objects) interp_matrix_k=*.dat *~