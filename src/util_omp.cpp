#include "../include/genisw.h"
#include <math.h>

using namespace std;

void partition(long int total, int processors, int *myid, long int *begin, long int *end){
  /* Partitions data according to the processor name.

  Parameters
  ----------
  long int total:
    Size of the full data to be partitioned.
  int processors:
    The number of processors being used for the MPI process.
  int myid:
    The current processors ID.
  long int begin:
    The integer for the end of the divided for-loop iterator.
  long int end:
    The integer for the end of the divided for-loop iterator.
  */
  long int partition_size = total/processors, partition_excess = total%processors;

  *begin = partition_size*(*myid);

  if(*myid == processors-1){
    *end = partition_size*(*myid+1)+partition_excess;
  }
  else{
    *end = partition_size*(*myid+1);
  }
}
