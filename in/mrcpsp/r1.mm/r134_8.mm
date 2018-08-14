************************************************************************
file with basedata            : cr134_.bas
initial value random generator: 1485037104
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  116
RESOURCES
  - renewable                 :  1   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       17        7       17
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6  13
   3        3          3           7  10  15
   4        3          3           6  13  17
   5        3          3           9  10  14
   6        3          2           8  15
   7        3          3           8  11  12
   8        3          1          14
   9        3          2          12  15
  10        3          2          11  12
  11        3          2          16  17
  12        3          1          17
  13        3          1          16
  14        3          1          16
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0
  2      1     2       6    3    2
         2     3       0    3    1
         3     4       6    2    1
  3      1     2       0    3    4
         2     9      10    2    4
         3    10       0    2    3
  4      1     1       8    7   10
         2     2       0    5    9
         3     7       0    4    8
  5      1     5       5    7    5
         2     9       0    6    4
         3     9       0    5    5
  6      1     6       8    9    3
         2     6       0    9    5
         3     9       7    9    1
  7      1     1       0    4    8
         2     2       5    3    4
         3     8       4    3    4
  8      1     6       0   10    6
         2     7       0   10    4
         3    10       0    9    4
  9      1     5       5    8    9
         2     7       0    7    7
         3     8       2    7    6
 10      1     3       9    6    3
         2     6       9    2    3
         3     6       0    4    3
 11      1     2       5    4    5
         2     5       0    2    2
         3     5       4    3    1
 12      1     1       3    5    5
         2     2       2    2    3
         3     3       0    2    2
 13      1     4       0    7    6
         2     8       0    6    3
         3    10       7    1    2
 14      1     2       8    6    3
         2     6       6    5    2
         3     8       5    4    2
 15      1     4       0    6    9
         2     6       7    5    5
         3     6       0    4    7
 16      1     1       0    7    2
         2     8       0    6    2
         3     8       0    7    1
 17      1     1       8    9    9
         2     1       0    6   10
         3     5       9    4    5
 18      1     0       0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  N 1  N 2
   18   75   62
************************************************************************
