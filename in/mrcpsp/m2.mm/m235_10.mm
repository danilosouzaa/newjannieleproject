************************************************************************
file with basedata            : cm235_.bas
initial value random generator: 1310693071
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  114
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       33        3       33
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        2          3           9  12  16
   3        2          2          13  14
   4        2          3           5   8  10
   5        2          3           6  11  14
   6        2          2           7  12
   7        2          2           9  17
   8        2          3           9  13  17
   9        2          1          15
  10        2          3          13  14  16
  11        2          2          12  16
  12        2          1          17
  13        2          1          15
  14        2          1          15
  15        2          1          18
  16        2          1          18
  17        2          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     3       7    0    8    6
         2     6       0    9    5    3
  3      1     7       0    5    4    9
         2     9       4    0    2    5
  4      1     7       0   10    4    8
         2     8       0    8    4    7
  5      1     4       1    0    4    2
         2     4       0    3    4    4
  6      1     6       0    6    6   10
         2    10       0    2    4    8
  7      1     7       5    0    8    9
         2     8       3    0    6    7
  8      1     1       0    2    2    1
         2     7       7    0    2    1
  9      1     6       4    0    2    6
         2     6       0    8    5    7
 10      1     3       0    9    9   10
         2     4       0    1    9    8
 11      1     9       6    0    5    2
         2     9       7    0    3    2
 12      1     4       9    0    6    3
         2     8       0    4    1    3
 13      1     1       0   10    6   10
         2     7       6    0    6   10
 14      1     6       0    5    5    6
         2    10       9    0    5    6
 15      1     3       1    0    7    6
         2     5       0    3    6    4
 16      1     4       5    0    4    9
         2     5       5    0    1    7
 17      1     7       6    0    9    8
         2     8       0    2    4    8
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   27   20   71   92
************************************************************************
