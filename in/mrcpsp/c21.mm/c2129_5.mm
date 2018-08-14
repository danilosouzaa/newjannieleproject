************************************************************************
file with basedata            : c2129_.bas
initial value random generator: 1618958250
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  115
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       31        1       31
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           7   8   9
   3        3          3           5  11  13
   4        3          3           6  10  13
   5        3          2           7   8
   6        3          2           8  11
   7        3          3          10  12  17
   8        3          3          12  14  17
   9        3          3          10  11  13
  10        3          1          14
  11        3          2          12  14
  12        3          2          15  16
  13        3          3          15  16  17
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
  2      1     1       7    6    0    5
         2     1       7    9    8    0
         3     5       4    6    0    5
  3      1     6       6    9    0    3
         2     6       7    9    9    0
         3     8       5    5    5    0
  4      1     1       4    3    0    6
         2     1       4    3    5    0
         3     4       1    2    0    5
  5      1     6       8    6    0    8
         2     7       6    5    0    8
         3     9       3    4    1    0
  6      1     3      10    6    0    8
         2     3       9    6    6    0
         3     7       9    6    0    7
  7      1     7       4    7    7    0
         2     8       2    7    0    4
         3     9       1    6    5    0
  8      1     1       7   10    9    0
         2     4       6    8    0    8
         3     5       6    7    9    0
  9      1     4       3   10    0    2
         2     7       3    9    7    0
         3     8       2    9    0    2
 10      1     1       7    5    2    0
         2    10       4    3    1    0
         3    10       5    2    0    2
 11      1     3       6    4    0    4
         2     7       5    4    0    4
         3     9       1    3    0    3
 12      1     6       9    8    7    0
         2     6       9   10    0    8
         3     7       8    8    7    0
 13      1     1       6    4    3    0
         2     1       9    4    0    6
         3     3       4    2    0    5
 14      1     8       4    9    4    0
         2     8       6    9    0    5
         3    10       4    6    0    5
 15      1     3       7    2    3    0
         2     8       3    2    0    8
         3     8       3    2    2    0
 16      1     1       9    8   10    0
         2     2       6    7    7    0
         3     5       4    5    2    0
 17      1     6       5    4    0    4
         2     8       4    3    8    0
         3     8       4    3    0    4
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   13   13   89   81
************************************************************************
