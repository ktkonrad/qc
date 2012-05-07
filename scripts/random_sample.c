/*
  script to take a random sample of entries in a file. single line, separated by commas
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHUNK_SIZE (1<<22) // 4MB
#define LEN 128
#define SEPARATOR ','

int main(int argc, char **argv) {
  FILE *infile, *outfile;
  int cutoff;
  char buf[CHUNK_SIZE+1]; // +1 for terminating null char
  char *current, *prev; // pointers to current and previous space in buf
  int bytes_read;

  if (argc < 4) {
    printf("usage: ./random_sample ratio infile outfile");
    exit(-1);
  }
  
  cutoff = atof(argv[1]) * RAND_MAX;
  infile = fopen(argv[2], "r");
  outfile = fopen(argv[3], "w");
  while (!feof(infile)) {
    bytes_read = fread(buf, sizeof(char), CHUNK_SIZE, infile);
    buf[CHUNK_SIZE] = (char)0;
    prev = buf;
    current = buf;
    while (current = strchr(prev, SEPARATOR)) {
      if (rand() < cutoff) {
        fwrite(prev, sizeof(char), current - prev + 1, outfile);
      }
      prev = current+1;
    }
  }
}
