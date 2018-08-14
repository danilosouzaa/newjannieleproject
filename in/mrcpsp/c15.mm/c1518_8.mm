************************************************************************
file with basedata            : c1518_.bas
initial value random generator: 1048813522
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
    1     16      0       19        3       19
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          1           9
   3        3          2          14  16
   4        3          3           5   8  10
   5        3          3           6   7   9
   6        3          3          11  12  13
   7        3          1          17
   8        3          1          16
   9        3          2          13  17
  10        3          1          12
  11        3          1          15
  12        3          1          17
  13        3          1          16
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
  2      1     2       7    0    3    0
         2     7       0    7    3    0
         3     9       0    5    0    8
  3      1     1       0    8    0    7
         2     4       0    4    9    0
         3     6       6    0    1    0
  4      1     4       0   10    6    0
         2     4      10    0    6    0
         3     8       5    0    5    0
  5      1     3       0    2    0    1
         2     3       5    0    9    0
         3     6       0    3    8    0
  6      1     1       0    7    0    7
         2     1      10    0    9    0
         3     7       0    8    8    0
  7      1     7       0    2    0   10
         2    10       3    0    0    8
         3    10       3    0    2    0
  8      1     5       0    6    0    8
         2     7       0    5    0    8
         3     8       6    0    7    0
  9      1     6       0    5    8    0
         2     8       4    0    6    0
         3     9       0    5    0    8
 10      1     3       4    0    0    3
         2     7       3    0   10    0
         3    10       1    0    9    0
 11      1     1       0    9    3    0
         2     7       4    0    2    0
         3     9       1    0    0    7
 12      1     5       0    8    0    6
         2     5       6    0    6    0
         3     8       5    0    1    0
 13      1     1       0    8    2    0
         2     3       7    0    1    0
         3     5       0    3    0    4
 14      1     4       4    0    6    0
         2     5       2    0    2    0
         3    10       2    0    0    6
 15      1     2       6    0    6    0
         2     7       5    0    0    6
         3    10       0    4    6    0
 16      1     5       4    0    0    6
         2     7       0    8    9    0
         3     8       0    8    0    3
 17      1     1       0    6    9    0
         2     1       3    0    0    9
         3     4       0    6    0    8
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   15   13   79   72
************************************************************************
