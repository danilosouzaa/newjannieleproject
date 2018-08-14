************************************************************************
file with basedata            : cr132_.bas
initial value random generator: 1218673837
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  119
RESOURCES
  - renewable                 :  1   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       10        4       10
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2          10  13
   3        3          3           6   8  14
   4        3          2           5  12
   5        3          3           6   8  14
   6        3          3           7   9  10
   7        3          1          15
   8        3          3           9  13  16
   9        3          2          11  17
  10        3          2          16  17
  11        3          1          15
  12        3          2          13  14
  13        3          2          15  17
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
  2      1     4       2    7    0
         2     4       3    5    0
         3    10       2    3    0
  3      1     2       8    0    8
         2     6       4    0    8
         3     8       4    2    0
  4      1     1       7    0    6
         2     1       6    4    0
         3     2       3    0    6
  5      1     1       3    5    0
         2     4       3    2    0
         3     6       2    0    8
  6      1     2       7    0    6
         2     2       6    8    0
         3     3       3    7    0
  7      1     3       5    6    0
         2     5       5    4    0
         3     9       5    0    2
  8      1     1       9    0    2
         2     6       7   10    0
         3     9       4    6    0
  9      1     3       7    7    0
         2     8       5    5    0
         3    10       5    0    8
 10      1     4       8    6    0
         2     7       5    0    4
         3     7       5    5    0
 11      1     1       8    9    0
         2     4       7    3    0
         3     9       5    0    1
 12      1     1       7    0    9
         2     3       6    2    0
         3     6       5    0    6
 13      1     2       8    6    0
         2     5       7    2    0
         3     9       6    0    5
 14      1     1       4    7    0
         2     9       3    0    7
         3    10       3    5    0
 15      1     2       7    8    0
         2     5       6    8    0
         3    10       3    0    6
 16      1     1       5    0    5
         2     4       5    7    0
         3     6       4    6    0
 17      1     1       8    2    0
         2     2       8    0    6
         3     5       8    0    2
 18      1     0       0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  N 1  N 2
   28   96   83
************************************************************************
