#! /bin/sh
#
# buildLAHCOpt.sh
# Copyright (C) 2016 haroldo <haroldo@soyuz>
#
# Distributed under terms of the MIT license.
#


CC=gcc
CPP=g++

CFLAGS="-O0 -g -flto -fopenmp -std=c99"
CPPFLAGS="-O0 -g -flto -fopenmp -DCPX -I/opt/ibm/ILOG/CPLEX_Studio1263/cplex/include/ilcplex/"
LDFLAGS="-O0 -g -L/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/x86-64_linux/static_pic/ -flto -lm -fopenmp -lpthread -lcplex"
CSOURCES="lslahc.c instance.c solution.c tokenizer.c vec_int.c vec_str.c \
long_compl_path.c str_utils.c vec_char.c node_heap.c top_sort.c memory.c mode_set.c \
neighborhood.c rrusage.c lahc.c learning.c test.c ms_solver_mip.c list_int.c"
CPPSOURCES="lp.cpp"
BINDIR=./bin/dbg/

rm $BINDIR/*

lnkFiles=""
echo building C sources ...
for cs in $CSOURCES;
do
    command="${CC} ${CFLAGS} -c $cs -o ${BINDIR}/`basename $cs .c`.o"
    printf "\t$command\n"
    lnkFiles="${lnkFiles}${BINDIR}/`basename $cs .c`.o "
    $command
done

echo building C++ sources ...
for cs in $CPPSOURCES;
do
    command="${CPP} ${CPPFLAGS} -c $cs -o ${BINDIR}/`basename $cs .cpp`.o"
    printf "\t$command\n"
    lnkFiles="${lnkFiles}${BINDIR}/`basename $cs .cpp`.o "
    $command
done

echo linking ...
command="${CPP}  ${lnkFiles} ${LDFLAGS} -o ${BINDIR}/psplahc"
printf "\t$command\n"
$command

