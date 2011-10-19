CC = g++
CFLAGS = -Wall -g
OFLAGS = -lm

objects = nodal_domain_driver.o count_nodal_domains.o random_percolation.o stack2.o util.o
dependencies = ../vergini/billiard.o ../vergini/useful.o ../vergini/nrutil.o

count: $(objects)
	$(CC) $(OFLAGS) -o count $(objects) $(dependencies)

nodal_domain_driver.o: nodal_domain_driver.c
	$(CC) $(CFLAGS) -c nodal_domain_driver.c

count_nodal_domains.o: count_nodal_domains.c count_nodal_domains.h
	$(CC) $(CFLAGS) -c count_nodal_domains.c

random_percolation.o: random_percolation.c random_percolation.h
	$(CC) $(CFLAGS) -c random_percolation.c

stack2.o: util/stack2.c util/stack2.h
	$(CC) $(CFLAGS) -c util/stack2.c

util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c

.PHONY : clean
clean:
	-rm count $(objects) *~