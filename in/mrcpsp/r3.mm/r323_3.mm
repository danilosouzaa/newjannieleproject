************************************************************************
file with basedata            : cr323_.bas
initial value random generator: 748814162
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  106
RESOURCES
  - renewable                 :  3   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       16        2       16
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6  11
   3        3          2           5   8
   4        3          3           8  10  17
   5        3          3           7  12  15
   6        3          2           9  10
   7        3          2          13  17
   8        3          2          15  16
   9        3          3          13  14  17
  10        3          2          12  14
  11        3          2          14  16
  12        3          1          13
  13        3          1          16
  14        3          1          15
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0
  2      1     1       4    3    8    0    8
         2     7       3    2    8    8    0
         3     9       2    2    7    6    0
  3      1     4       5    3    9    0    5
         2     4       5    4    9    8    0
         3    10       5    3    6    0    5
  4      1     4       7    7    4    7    0
         2     7       7    7    3    0    6
         3     9       4    3    2    3    0
  5      1     1       6    7    5    0    7
         2     1       7    6    5    4    0
         3     2       5    5    4    0    7
  6      1     1      10    2    5    0    3
         2     1       8    4    5    0    5
         3     1       7    2    4    2    0
  7      1     4       8    6    7    0    7
         2     8       8    4    4    7    0
         3     8       8    5    1    0    7
  8      1     1       7    6    8    5    0
         2     4       7    6    7    0    1
         3     6       6    6    7    3    0
  9      1     1       4    9    6   10    0
         2     2       3    8    6    9    0
         3    10       2    5    4    6    0
 10      1     2       7    8    9    3    0
         2     2       8    8    9    0    3
         3     6       5    7    9    0    3
 11      1     3       5    7    4    5    0
         2     6       5    4    3    5    0
         3     6       4    5    1    0    8
 12      1     1       5    6    7    0    1
         2     1       6    8    7    4    0
         3     3       4    3    6    2    0
 13      1     3       8    3    6    0    5
         2     3       9    3    7    0    4
         3     4       8    2    5    6    0
 14      1     5       8    8    6    0   10
         2     8       5    7    5    9    0
         3    10       2    7    4    0    9
 15      1     2       6    8    7    5    0
         2     3       6    7    6    2    0
         3     7       2    4    6    2    0
 16      1     4       4    5    5    6    0
         2     5       2    4    2    6    0
         3     5       1    3    5    0    1
 17      1     7       6    7    6    0    4
         2     7       5    5    7    0    3
         3    10       2    3    4    0    3
 18      1     0       0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  N 1  N 2
   22   22   20   69   54
************************************************************************
