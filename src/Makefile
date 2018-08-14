CC=clang
CXX=clang++

# Build Options
#CFLAGS=-Wall -O0 -g3 -std=c99 -Wall -DDEBUG
#CXXFLAGS=-Wall -O0 -g3 -Wall -DDEBUG
CFLAGS=-Wall -O0 -g3 -std=c99 -Wall -DDEBUG -fsanitize=address
CXXFLAGS=-Wall -O0 -g3 -Wall -DDEBUG -fsanitize=address
LDFLAGS=-O0 -g3 -fsanitize=address
#SOLVER_FLAGS=-DCPX -I/opt/ibm/ILOG/CPLEX_Studio1263/cplex/include/ilcplex/
#SOLVER_LDFLAGS=-L/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/x86-64_linux/static_pic/ -lcplex -lpthread -lm
SOLVER_FLAGS=-DCBC `pkg-config --cflags cbc`
SOLVER_LDFLAGS=`pkg-config --libs cbc`
BINDIR=./bin/dbg/

all:mip

clean:
	rm -f ./bin/dbg/* ; rm -f ./bin/opt/* ; rm vgcore* -f ; rm core* -f

# LINKING
iinfo:$(BINDIR)/iinfo.o
	$(CC)  $(BINDIR)/iinfo.o $(BINDIR)/instance.o $(BINDIR)/tokenizer.o $(BINDIR)/str_utils.o $(BINDIR)/vec_char.o \
	      $(BINDIR)/vec_str.o $(BINDIR)/vec_int.o  $(BINDIR)/memory.o $(LDFLAGS) -o $(BINDIR)/iinfo

dcg:$(BINDIR)/dcg.o
	$(CXX) $(LDFLAGS) $(SOLVER_LDFLAGS) $(BINDIR)/dcg.o $(BINDIR)/instance.o $(BINDIR)/tokenizer.o $(BINDIR)/str_utils.o $(BINDIR)/vec_char.o \
	      $(BINDIR)/vec_str.o $(BINDIR)/vec_int.o  $(BINDIR)/memory.o $(BINDIR)/mip_compact.o  $(BINDIR)/solution.o $(BINDIR)/rrusage.o \
           $(BINDIR)/node_heap.o $(BINDIR)/lp.o  $(BINDIR)/top_sort.o $(BINDIR)/long_compl_path.o $(BINDIR)/vec_double.o  $(BINDIR)/mode_set.o\
		  $(BINDIR)/dict_int.o -lgomp \
		  -o $(BINDIR)/dcg

mip:$(BINDIR)/mip.o
	$(CXX)   $(BINDIR)/mip.o $(BINDIR)/instance.o $(BINDIR)/tokenizer.o $(BINDIR)/str_utils.o $(BINDIR)/vec_char.o \
	      $(BINDIR)/vec_str.o $(BINDIR)/vec_int.o  $(BINDIR)/memory.o $(BINDIR)/mip_compact.o \
          $(BINDIR)/lp.o  $(BINDIR)/top_sort.o $(BINDIR)/long_compl_path.o $(BINDIR)/vec_double.o $(BINDIR)/ms_solver_mip.o\
		  $(BINDIR)/dict_int.o $(BINDIR)/solution.o $(BINDIR)/mode_set.o $(BINDIR)/rrusage.o $(BINDIR)/node_heap.o $(LDFLAGS) $(SOLVER_LDFLAGS) -lgomp \
		  -o $(BINDIR)/mip

# COMPILING
$(BINDIR)/mip.o:mip.c $(BINDIR)/instance.o $(BINDIR)/mip_compact.o
	$(CC) -c $(CFLAGS) mip.c -o $(BINDIR)/mip.o

$(BINDIR)/dcg.o:dcg.c $(BINDIR)/instance.o $(BINDIR)/mip_compact.o
	$(CC) -c $(CFLAGS) dcg.c -o $(BINDIR)/dcg.o

$(BINDIR)/mip_compact.o:mip_compact.c $(BINDIR)/instance.o $(BINDIR)/lp.o $(BINDIR)/solution.o\
         $(BINDIR)/dict_int.o $(BINDIR)/vec_str.o  $(BINDIR)/ms_solver_mip.o
	$(CC) -c $(CFLAGS) mip_compact.c -o $(BINDIR)/mip_compact.o

$(BINDIR)/ms_solver_mip.o:ms_solver_mip.c $(BINDIR)/instance.o $(BINDIR)/lp.o $(BINDIR)/solution.o\
         $(BINDIR)/dict_int.o $(BINDIR)/vec_str.o
	$(CC) -c $(CFLAGS) ms_solver_mip.c -o $(BINDIR)/ms_solver_mip.o

$(BINDIR)/iinfo.o:iinfo.c $(BINDIR)/instance.o
	$(CC) -c $(CFLAGS) iinfo.c -o $(BINDIR)/iinfo.o

$(BINDIR)/lp.o:lp.cpp $(BINDIR)/instance.o
	$(CXX) -c $(CXXFLAGS) -fopenmp $(SOLVER_FLAGS) lp.cpp -o $(BINDIR)/lp.o

$(BINDIR)/instance.o:instance.c instance.h macros.h $(BINDIR)/vec_int.o $(BINDIR)/vec_str.o \
             $(BINDIR)/tokenizer.o $(BINDIR)/str_utils.o $(BINDIR)/long_compl_path.o
	$(CC) -c $(CFLAGS) instance.c -o $(BINDIR)/instance.o

$(BINDIR)/long_compl_path.o:long_compl_path.c long_compl_path.h $(BINDIR)/vec_double.o $(BINDIR)/top_sort.o
	$(CC) -c $(CFLAGS) long_compl_path.c -o $(BINDIR)/long_compl_path.o

$(BINDIR)/top_sort.o:top_sort.c top_sort.h 
	$(CC) -c $(CFLAGS) top_sort.c -o $(BINDIR)/top_sort.o

$(BINDIR)/vec_double.o:vec_double.c vec_double.h 
	$(CC) -c $(CFLAGS) vec_double.c -o $(BINDIR)/vec_double.o

$(BINDIR)/solution.o:solution.c solution.h macros.h $(BINDIR)/instance.o $(BINDIR)/mode_set.o $(BINDIR)/node_heap.o
	$(CC) -c $(CFLAGS) solution.c -o $(BINDIR)/solution.o

$(BINDIR)/node_heap.o:node_heap.c node_heap.h macros.h 
	$(CC) -c $(CFLAGS) node_heap.c -o $(BINDIR)/node_heap.o

$(BINDIR)/mode_set.o:mode_set.c mode_set.h macros.h $(BINDIR)/rrusage.o
	$(CC) -c $(CFLAGS) mode_set.c -o $(BINDIR)/mode_set.o

$(BINDIR)/rrusage.o:rrusage.c rrusage.h macros.h
	$(CC) -c $(CFLAGS) rrusage.c -o $(BINDIR)/rrusage.o

$(BINDIR)/vec_char.o:vec_char.c vec_char.h macros.h 
	$(CC) -c $(CFLAGS) vec_char.c -o $(BINDIR)/vec_char.o

$(BINDIR)/vec_int.o:vec_int.c vec_int.h macros.h 
	$(CC) -c $(CFLAGS) vec_int.c -o $(BINDIR)/vec_int.o

$(BINDIR)/dict_int.o:dict_int.c dict_int.h macros.h $(BINDIR)/vec_str.o
	$(CC) -c $(CFLAGS) dict_int.c -o $(BINDIR)/dict_int.o

$(BINDIR)/str_utils.o:str_utils.c str_utils.h macros.h 
	$(CC) -c $(CFLAGS) str_utils.c -o $(BINDIR)/str_utils.o

$(BINDIR)/vec_str.o:vec_str.c vec_str.h macros.h $(BINDIR)/memory.o $(BINDIR)/str_utils.o
	$(CC) -c $(CFLAGS) vec_str.c -o $(BINDIR)/vec_str.o

$(BINDIR)/tokenizer.o:tokenizer.c tokenizer.h macros.h $(BINDIR)/vec_int.o $(BINDIR)/vec_str.o $(BINDIR)/vec_char.o
	$(CC) -c $(CFLAGS) tokenizer.c -o $(BINDIR)/tokenizer.o

$(BINDIR)/memory.o:memory.c memory.h macros.h 
	$(CC) -c $(CFLAGS) memory.c -o $(BINDIR)/memory.o
