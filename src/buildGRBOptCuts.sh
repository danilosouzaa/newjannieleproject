#! /bin/sh
#
# buildLAHCOpt.sh
# Copyright (C) 2016 haroldo <haroldo@soyuz>
#
# Distributed under terms of the MIT license.
#


CC=gcc
CPP=g++


CFLAGS="-Ofast -flto -fopenmp -std=c99 -I/opt/gurobi752/linux64/include/"
CPPFLAGS="-Ofast -flto -fopenmp -DGRB -I/opt/gurobi752/linux64/include/ "
LDFLAGS="-Ofast -flto -fopenmp -lpthread -lm -L/opt/gurobi752/linux64/lib/ -lgurobi75"
CSOURCES=" ippsp.c instance.c solution.c tokenizer.c vec_int.c vec_str.c long_compl_path.c str_utils.c vec_char.c node_heap.c top_sort.c memory.c mode_set.c rrusage.c ms_solver_mip.c list_int.c mip_compact.c stack.c vec_double.c cgraph.c clique_enum.c clique_extender.c clique_separation.c clique.c  dict_int.c grasp.c oddhs.c spaths.c vectormgm.c vint_queue.c vint_set.c cut_pool.c results.c clique_elite_set.c parameters.c cut_cg.c cut_clique.c cut_cover.c cut_precedence.c cut_gpu.c solutionGpu.c cut_oddholes.c prepareCPU.c"
CPPSOURCES=" lp.cpp BKGraph.cpp BKVertex.cpp bron_kerbosch.cpp build_cgraph.cpp"
BINDIR=./bin/opt/

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
    command="${CPP}  ${lnkFiles} ${LDFLAGS} -o ${BINDIR}/ippsp"
    printf "\t$command\n"
    $command

