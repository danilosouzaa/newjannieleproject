************************************************************************
file with basedata            : cn142_.bas
initial value random generator: 1030208970
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  132
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  1   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       28       11       28
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   9  12
   3        3          3           5   6  10
   4        3          3           6  10  13
   5        3          3          11  15  17
   6        3          3           7  15  17
   7        3          2           8   9
   8        3          3          11  12  14
   9        3          2          11  14
  10        3          1          12
  11        3          1          16
  12        3          1          16
  13        3          1          14
  14        3          1          16
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1
------------------------------------------------------------------------
  1      1     0       0    0    0
  2      1     6       0    2    7
         2     6       1    0    6
         3     8       0    2    3
  3      1     4       0    7    9
         2     4       0    5   10
         3     7       0    4    7
  4      1     5       0    3    7
         2     6       6    0    5
         3     7       3    0    4
  5      1     7       8    0    6
         2     7       0    1    6
         3     8       7    0    6
  6      1     3       0    1    7
         2     8       0    1    5
         3     9       0    1    4
  7      1     1       1    0    6
         2     7       0    5    4
         3    10       0    2    3
  8      1     3       8    0    5
         2     5       7    0    4
         3     8       0    5    3
  9      1     1       5    0    7
         2     5       0    2    7
         3     7       2    0    2
 10      1     4       4    0    4
         2     7       4    0    3
         3     9       3    0    1
 11      1     7       5    0    9
         2     9       0    8    5
         3     9       5    0    3
 12      1     7       0    6    8
         2     8       5    0    7
         3    10       4    0    7
 13      1     5       6    0    4
         2     6       0    4    4
         3     8       0    1    4
 14      1     1       0    5    9
         2     4       0    3    9
         3     5       0    1    8
 15      1     5      10    0    4
         2     6       0    3    4
         3     8       0    2    1
 16      1     8       0    6    4
         2     9       0    6    3
         3    10       0    5    2
 17      1     2       7    0    9
         2     8       0    6    7
         3     9       0    5    5
 18      1     0       0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1
   15   12   85
************************************************************************
