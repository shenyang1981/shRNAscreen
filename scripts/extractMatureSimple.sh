#!/usr/bin/env bash

# first arg is raw fastq
# second arg is libname
# third arg is output
FASTAQ=$1
LIB=$2
OUTPUT=$3
zcat ${FASTAQ} | paste - - - -| mawk -v lib=${LIB} 'BEGIN{FS="\t";OFS="\t"} \
  { \
    mature = substr($2,3,19); \
    a[mature]++; \
  } END { \
    for(i in a){ \
	print lib,i,a[i]; \
    } \
}' | pigz -p 32 -c > ${OUTPUT}
