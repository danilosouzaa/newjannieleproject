************************************************************************
file with basedata            : me23_.bas
initial value random generator: 242155063
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  16
horizon                       :  105
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  0   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     14      0       16        5       16
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6   7
   3        3          3           5   9  11
   4        3          2          12  13
   5        3          2           8  14
   6        3          2          10  15
   7        3          3           8  10  12
   8        3          1          15
   9        3          2          10  12
  10        3          1          13
  11        3          3          13  14  15
  12        3          1          14
  13        3          1          16
  14        3          1          16
  15        3          1          16
  16        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2
------------------------------------------------------------------------
  1      1     0       0    0
  2      1     2       8    7
         2     3       7    5
         3    10       6    2
  3      1     4       7    5
         2     7       5    5
         3     8       4    5
  4      1     4       9    7
         2     7       8    7
         3     8       6    4
  5      1     3       6   10
         2     5       5    6
         3     5       6    4
  6      1     4       7    4
         2     7       5    2
         3    10       2    2
  7      1     5       6    7
         2     6       3    7
         3     6       4    6
  8      1     3       8    8
         2     4       4    8
         3     4       5    7
  9      1     3       2    8
         2     8       1    7
         3     9       1    6
 10      1     5       8    5
         2     9       6    5
         3    10       3    2
 11      1     6       9    5
         2     6       8    7
         3     7       2    3
 12      1     2       5    8
         2     5       3    5
         3     5       2    6
 13      1     3       9    2
         2     7       8    2
         3    10       6    1
 14      1     6       8    4
         2     7       1    4
         3     7       6    3
 15      1     4       5    8
         2     6       5    3
         3     6       4    5
 16      1     0       0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2
   24   29
************************************************************************
