************************************************************************
file with basedata            : md206_.bas
initial value random generator: 31122
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  134
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       21       15       21
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3          12  13  16
   3        3          3           5   9  16
   4        3          3           5   8   9
   5        3          3           6   7  12
   6        3          1          13
   7        3          1          10
   8        3          3          11  16  17
   9        3          2          11  14
  10        3          3          11  14  17
  11        3          1          13
  12        3          2          14  17
  13        3          1          15
  14        3          1          15
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       4    6    0    6
         2     6       3    5    9    0
         3     7       3    5    2    0
  3      1     3       6    6   10    0
         2     4       3    4   10    0
         3    10       2    1    0    4
  4      1     4       4    7    0    5
         2     6       3    5    0    4
         3     6       3    6    9    0
  5      1     3       8    7    0    7
         2     8       6    6    8    0
         3    10       4    1    0    5
  6      1     1       9    6    0   10
         2     1       8    7    4    0
         3     5       6    5    0   10
  7      1     3       8    7    0    9
         2     5       5    7    0    3
         3     8       3    3    8    0
  8      1     1      10    8    0    8
         2     2       8    5    8    0
         3     2       9    6    6    0
  9      1     1       9    9    7    0
         2     4       6    8    0    3
         3    10       2    6    6    0
 10      1     1       4    6    0    7
         2     4       4    5    1    0
         3    10       2    2    0    5
 11      1     2       9    9    0    6
         2     6       9    9    0    5
         3    10       8    8    6    0
 12      1     3       9    8    3    0
         2     9       9    5    2    0
         3    10       8    3    1    0
 13      1     1      10    6    0    5
         2     7       8    3    0    3
         3    10       3    2    0    1
 14      1     1      10    7    0    6
         2     4      10    6   10    0
         3     9       9    2   10    0
 15      1     7      10    8    6    0
         2     8      10    6    5    0
         3     9      10    2    4    0
 16      1     3       4    4    9    0
         2     6       4    4    0    4
         3    10       2    2    8    0
 17      1     3       2    6    8    0
         2     4       1    5    6    0
         3     8       1    5    0    9
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   20   17   56   45
************************************************************************
