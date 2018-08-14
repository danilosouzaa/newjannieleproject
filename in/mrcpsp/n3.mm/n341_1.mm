************************************************************************
file with basedata            : cn341_.bas
initial value random generator: 17310
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  124
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  3   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       18        8       18
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6  14
   3        3          2          10  13
   4        3          3           9  12  16
   5        3          2           8  12
   6        3          3           7  10  12
   7        3          3          11  13  17
   8        3          2           9  16
   9        3          3          10  11  13
  10        3          1          17
  11        3          1          15
  12        3          1          15
  13        3          1          15
  14        3          2          16  17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2  N 3
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0
  2      1     2       0    9    8   10    7
         2     6       0    7    8    9    7
         3     7       7    0    7    8    7
  3      1     1       5    0    5    7    9
         2     1       0    5    6    7    9
         3     5       6    0    4    7    8
  4      1     2       4    0    6    9    6
         2     3       0   10    5    9    5
         3     6       3    0    4    9    4
  5      1     1       0    6    5    3    7
         2     1       6    0    5    4    7
         3     7       0    6    5    3    3
  6      1     2       0    7    6    6    4
         2     2       6    0    9    6    5
         3     9       6    0    3    5    4
  7      1     4       5    0    4    4    9
         2     7       3    0    3    4    9
         3    10       0    8    2    3    9
  8      1     1       0    4    6    9   10
         2     7       0    3    5    7    7
         3     7       9    0    4    9    7
  9      1     5       0    4    6    6    8
         2     8       0    3    6    5    7
         3    10       5    0    5    4    7
 10      1     1       0    7    5    6    8
         2     3       0    6    4    5    8
         3    10       3    0    4    5    7
 11      1     2       0    8    6   10    3
         2     6       0    7    5    8    3
         3     8       6    0    3    8    2
 12      1     3       8    0    9    7    3
         2     3       9    0    8    8    3
         3     7       8    0    8    7    3
 13      1     6       0    8    5    8   10
         2     7       0    7    2    5   10
         3     9       0    4    1    5    9
 14      1     3       4    0    5    8    5
         2     5       0   10    5    8    5
         3     7       3    0    2    7    2
 15      1     3       5    0    9    6    4
         2     4       3    0    7    6    4
         3     6       3    0    3    5    3
 16      1     5       7    0    4    5    6
         2     8       0    8    2    3    4
         3     8       6    0    2    1    4
 17      1     6       1    0    3   10    5
         2     6       2    0    3    7    5
         3     8       0    7    3    4    3
 18      1     0       0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2  N 3
   13   11   78  102   94
************************************************************************
