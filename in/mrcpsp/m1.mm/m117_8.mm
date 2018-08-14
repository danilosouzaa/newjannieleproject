************************************************************************
file with basedata            : cm117_.bas
initial value random generator: 719005379
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  66
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       31        0       31
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        1          2          10  13
   3        1          3           5   7  14
   4        1          2          11  12
   5        1          3           6   8  11
   6        1          1          15
   7        1          3           9  10  12
   8        1          3          10  12  13
   9        1          2          11  13
  10        1          2          15  17
  11        1          1          16
  12        1          2          16  17
  13        1          1          16
  14        1          2          15  17
  15        1          1          18
  16        1          1          18
  17        1          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       0    9    0    6
  3      1     9       0    3    3    0
  4      1     5       7    0    0    2
  5      1     5       6    0    0    2
  6      1     3       0    4    0    2
  7      1     3       5    0    0    7
  8      1     4       5    0    5    0
  9      1     4       0    7    0    4
 10      1     6       1    0    7    0
 11      1     5       4    0    0    1
 12      1     4       0    7    7    0
 13      1     2       4    0    0    8
 14      1     4       2    0    2    0
 15      1     7       4    0    0    6
 16      1     2       0    8    0    3
 17      1     2       0    3    4    0
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
    9   10   28   41
************************************************************************
