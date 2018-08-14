#/bin/sh
cppcheck --std=c99 --library=/usr/share/cppcheck/std.cfg  -DDEBUG -I/opt/ibm/ILOG/CPLEX_Studio127/cplex/include/ -DCPX *.c*
cppcheck-htmlreport --file=CppCheckResults.xml  --report-dir=cppcheck
