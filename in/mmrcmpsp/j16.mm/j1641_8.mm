************************************************************************
file with basedata            : md233_.bas
initial value random generator: 131745160
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  128
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       26        1       26
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6   7
   3        3          2           7   9
   4        3          3           7  10  12
   5        3          3           8  11  12
   6        3          2           9  10
   7        3          1          17
   8        3          3          10  14  17
   9        3          3          15  16  17
  10        3          1          13
  11        3          1          13
  12        3          2          13  16
  13        3          1          15
  14        3          2          15  16
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     2       6    0    4    3
         2     7       5    0    3    3
         3     8       5    0    3    2
  3      1     1       0    7    8    9
         2     3       9    0    8    9
         3     5       8    0    8    6
  4      1     8       2    0    3    3
         2     8       0    6    5    3
         3    10       0    4    1    2
  5      1     3       9    0    5    5
         2     6       7    0    5    3
         3     9       0    2    3    2
  6      1     2       0    4    5    4
         2     4       8    0    5    3
         3     9       0    2    2    2
  7      1     1       0    4    5    6
         2     4       0    3    5    6
         3     9       6    0    3    4
  8      1     7       3    0    2    2
         2     7       0    7    3    2
         3    10       0    6    2    2
  9      1     4       5    0   10    7
         2     5       0    7    7    7
         3     7       4    0    2    6
 10      1     6       0    5    6    6
         2     7       0    5    5    5
         3     7       4    0    5    4
 11      1     1       8    0    2    8
         2     3       5    0    1    8
         3     9       0    5    1    6
 12      1     3       9    0    5    5
         2     3       0    5    5    5
         3     5      10    0    5    4
 13      1     7       9    0    7    6
         2     9       4    0    7    6
         3     9       0    2    7    6
 14      1     1       0    8    7    8
         2     5       0    6    6    8
         3     8       7    0    2    6
 15      1     1      10    0    8    6
         2     4       9    0    7    5
         3     5       0    4    7    5
 16      1     4       4    0    3   10
         2     9       0    8    2    9
         3    10       0    8    2    8
 17      1     8       0    4    5    4
         2     8       7    0    4    7
         3     8       9    0    5    3
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   11    6   73   82
************************************************************************
