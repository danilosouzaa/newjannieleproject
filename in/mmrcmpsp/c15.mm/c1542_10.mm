************************************************************************
file with basedata            : c1542_.bas
initial value random generator: 1572333661
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  127
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       21        8       21
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2           6  10
   3        3          2           5   8
   4        3          1          16
   5        3          3           9  10  13
   6        3          2           7  17
   7        3          1           9
   8        3          1          11
   9        3          2          12  16
  10        3          1          11
  11        3          1          17
  12        3          1          15
  13        3          2          14  17
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
  2      1     6       0    6   10    8
         2     7       9    0   10    8
         3    10       0    4   10    8
  3      1     2       0    8    7    8
         2     9       7    0    5    6
         3     9       0    7    6    6
  4      1     2       9    0    8   10
         2     3       0    8    7    8
         3     4       9    0    7    7
  5      1     1       0    5    9   10
         2     2       9    0    7   10
         3     7       4    0    3    9
  6      1     1       9    0    4    2
         2     8       9    0    3    1
         3    10       0    1    2    1
  7      1     5       0    5    6    7
         2     8       0    4    5    5
         3    10       8    0    5    4
  8      1     2       0    6    2    8
         2     4       0    5    2    6
         3     7       1    0    1    5
  9      1     5       0    8    7    9
         2     6       0    5    4    9
         3     6       4    0    6    8
 10      1     5       0    7    3    1
         2     7       1    0    3    1
         3     7       0    5    3    1
 11      1     1       0    3    8    2
         2     5       0    3    8    1
         3     9       3    0    5    1
 12      1     2       5    0    7    7
         2     4       0   10    7    4
         3     5       0    3    6    2
 13      1     1       0    4    8    9
         2     7       0    2    7    9
         3     8       0    2    7    8
 14      1     1       7    0    3    7
         2     4       3    0    3    7
         3    10       0    4    3    5
 15      1     2       4    0    8    8
         2     7       2    0    8    7
         3     9       0    7    3    4
 16      1     1       7    0   10    6
         2     4       7    0    9    4
         3     6       5    0    9    3
 17      1     5       0    8    4    9
         2     8       7    0    4    7
         3    10       0    8    3    7
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   15   14   90   95
************************************************************************
