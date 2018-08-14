************************************************************************
file with basedata            : md187_.bas
initial value random generator: 1495739892
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  16
horizon                       :  115
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     14      0       18        0       18
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   8  13
   3        3          2           6   9
   4        3          3           6   7   8
   5        3          3           6   9  10
   6        3          2          14  15
   7        3          2          12  13
   8        3          2           9  10
   9        3          1          11
  10        3          2          12  15
  11        3          2          12  15
  12        3          1          14
  13        3          1          16
  14        3          1          16
  15        3          1          16
  16        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     3       0   10    5    7
         2     3       4    0    5    6
         3     8       2    0    4    2
  3      1     1       8    0    9    6
         2     3       3    0    8    6
         3     9       0    6    5    5
  4      1     3       6    0    5    9
         2     3       0    7    6    9
         3    10       0    5    2    8
  5      1     2       5    0    6    2
         2     4       4    0    5    1
         3     7       0    9    4    1
  6      1     2       5    0    8    4
         2     5       0    9    7    2
         3     9       0    9    2    1
  7      1     1       3    0    7    5
         2     1       5    0    6    5
         3     1       6    0    6    3
  8      1     4       1    0    8    8
         2     6       0    4    7    5
         3     7       0    1    6    3
  9      1     3       4    0    6    8
         2     6       4    0    5    7
         3    10       2    0    4    7
 10      1     1       9    0    4    3
         2     9       8    0    3    3
         3    10       4    0    2    2
 11      1     3       6    0    6    9
         2     7       0    9    4    8
         3     8       0    9    3    8
 12      1     3       6    0    7    4
         2     8       5    0    7    4
         3    10       0    1    4    4
 13      1     6       0    9    7    5
         2     7       5    0    7    4
         3     8       5    0    5    4
 14      1     2       4    0    5    8
         2    10       3    0    5    3
         3    10       0    3    5    5
 15      1     1       0    5    8    5
         2     7       0    4    5    5
         3     8       4    0    5    3
 16      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   15   21   92   83
************************************************************************
