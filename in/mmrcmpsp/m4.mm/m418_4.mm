************************************************************************
file with basedata            : cm418_.bas
initial value random generator: 472947272
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  141
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       15       13       15
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        4          3           5  10  11
   3        4          3           6   8  15
   4        4          3           7   8   9
   5        4          3          12  13  17
   6        4          2          11  13
   7        4          2          14  15
   8        4          1          17
   9        4          3          11  12  17
  10        4          2          12  14
  11        4          1          14
  12        4          2          15  16
  13        4          1          16
  14        4          1          16
  15        4          1          18
  16        4          1          18
  17        4          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       4    0    4    0
         2     3       0    8    4    0
         3     4       4    0    0    6
         4     9       0    7    0    2
  3      1     4       0    4    0    4
         2     5       3    0    9    0
         3     5       0    3    0    1
         4     6       0    1   10    0
  4      1     1      10    0    0    4
         2     2       7    0    6    0
         3     8       6    0    0    3
         4     9       4    0    0    2
  5      1     6       9    0    0    7
         2     6       0    9    8    0
         3     9       9    0    6    0
         4    10       0    6    0    7
  6      1     4       0    8    5    0
         2     5       0    7    4    0
         3     8       8    0    0    9
         4    10       0    6    3    0
  7      1     3       0    4    0    9
         2     3       4    0    4    0
         3     4       3    0    0    9
         4     9       0    4    4    0
  8      1     1       9    0    6    0
         2     7       7    0    6    0
         3     9       0    3    0    4
         4    10       4    0    0    2
  9      1     1       9    0    8    0
         2     1       8    0    0    5
         3     2       6    0    0    2
         4     6       6    0    8    0
 10      1     2       0    8    8    0
         2     8       4    0    0    7
         3     9       0    8    6    0
         4    10       0    8    2    0
 11      1     2       0    5    0    3
         2     3       0    4    0    3
         3     8       0    3    0    2
         4     9       8    0    7    0
 12      1     2       0    5   10    0
         2     4       3    0    0    9
         3     5       2    0    0    8
         4    10       2    0   10    0
 13      1     3       4    0    0    7
         2     4       0    6    0    3
         3     6       3    0    8    0
         4     8       2    0    8    0
 14      1     2       0    7    0    7
         2     6       5    0    4    0
         3     6       5    0    0    5
         4    10       4    0    5    0
 15      1     3       0    2    8    0
         2     6       9    0    0    8
         3     7       6    0    8    0
         4    10       5    0    7    0
 16      1     3       0    8    5    0
         2     5       0    7    0    9
         3     9       0    6    0    5
         4     9       1    0    0    7
 17      1     1       9    0    4    0
         2     1       0    5    0    4
         3     4       0    3    0    3
         4     6       7    0    0    3
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   18   14   80   77
************************************************************************
