************************************************************************
file with basedata            : cm250_.bas
initial value random generator: 351414427
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  112
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       22        5       22
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        2          2           7   9
   3        2          3           5   8  11
   4        2          3           6   8  12
   5        2          3           6  10  12
   6        2          1          16
   7        2          3           8  11  13
   8        2          1          16
   9        2          3          10  11  15
  10        2          2          14  16
  11        2          1          17
  12        2          2          15  17
  13        2          2          14  15
  14        2          1          17
  15        2          1          18
  16        2          1          18
  17        2          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     3       0    3    4    7
         2     5       4    0    3    1
  3      1     8       0    5    7    9
         2     9       0    1    5    5
  4      1     5       6    0    7    4
         2     9       6    0    7    1
  5      1     3       0    7    6    7
         2     7       7    0    4    7
  6      1     5       0    3    7    6
         2     6       5    0    5    6
  7      1     5       0    5    6    9
         2     7       2    0    1    6
  8      1     4       0    5    6    6
         2     5       4    0    5    3
  9      1     7      10    0    5    9
         2    10       4    0    4    7
 10      1     4       6    0    9    4
         2     4       0    9    9    4
 11      1     5       2    0    9    3
         2     7       2    0    7    3
 12      1     5       5    0    9   10
         2    10       0    7    9    6
 13      1     5       0    9    8    4
         2     9       0    1    6    4
 14      1     1       0    2    5    9
         2     7       2    0    5    8
 15      1     6       2    0    7    8
         2     7       0    2    4    7
 16      1     1       7    0    8    7
         2     2       4    0    7    4
 17      1     4       9    0    6    8
         2     8       0    7    5    4
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   14   17  103  102
************************************************************************
