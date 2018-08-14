************************************************************************
file with basedata            : cr525_.bas
initial value random generator: 24763
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  124
RESOURCES
  - renewable                 :  5   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       18        0       18
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           6   7   8
   3        3          2           5  10
   4        3          3          13  16  17
   5        3          3           6   9  11
   6        3          3          12  13  14
   7        3          3           9  14  17
   8        3          2          12  14
   9        3          1          12
  10        3          2          11  13
  11        3          1          15
  12        3          1          16
  13        3          1          15
  14        3          2          15  16
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  R 4  R 5  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0    0    0
  2      1     5       0    0    8    0    0    0    5
         2     8       0   10    0    0    0    7    0
         3     9       8    6    7    0    0    0    5
  3      1     2       0    0    0    0    7    0   10
         2     2       0    2    4    4    0    0    9
         3     4       1    0    0    0    0    0    7
  4      1     7       2    0    0    0    4    0    6
         2     7       0    9    0    0    7    8    0
         3    10       2    6    0    0    0    0    8
  5      1     3       6    1    0    0    0    0    9
         2     3       4    4    7    7    6    0    9
         3     6       0    0    0    4    4    0    8
  6      1     4       3    0   10    3    4    0    5
         2     4       6    0    0    3    0    7    0
         3     7       0    0    0    3    4    0    5
  7      1     4       0    0    0    6    5    0    3
         2     6       6    0   10    0    0    6    0
         3    10       0    2   10    0    0    5    0
  8      1     2       0    0    7    4    3    0    8
         2     4       0    0    6    0    0    0    5
         3     4       0    2    7    3    2    7    0
  9      1     1       3    0    0    0    0    7    0
         2     1       2    8    8    5    6    0    6
         3    10       0    0    8    5    0    8    0
 10      1     5       9    0    0    0    9    0    6
         2     9       0    0    1    1    9    0    6
         3     9       0    7    1    0    9    4    0
 11      1     4       8    9    9    6    0    7    0
         2     6       0    0    8    4    0    1    0
         3    10       0    0    0    0    6    0    5
 12      1     1       2    0    7    0    0    0    5
         2     2       2    0    6    4    3    7    0
         3     9       0    5    6    0    2    7    0
 13      1     1       0    8    0    4    7    8    0
         2     1       0    7    3    0    7    2    0
         3     6       3    6    3    0    7    0    7
 14      1     3       0    0    5    9    0    0    8
         2     7       8    0    5    6    8    5    0
         3     7       0    2    4    0    0    0    7
 15      1     6       0    0   10    0    7    9    0
         2     8       0    0   10    5    0    8    0
         3    10       0   10    9    5    4    0    4
 16      1     3       0    8    0    3    0    0    2
         2     3       8    5    0    0    9    7    0
         3     4       0    4    9    0    9    3    0
 17      1     3       0    4    8    7    0    3    0
         2     5       0    0    0    6    7    2    0
         3     9       0    3    7    0    6    0    2
 18      1     0       0    0    0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  R 4  R 5  N 1  N 2
    6   12   15   10   15   93   93
************************************************************************
