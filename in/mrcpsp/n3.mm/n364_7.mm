************************************************************************
file with basedata            : cn364_.bas
initial value random generator: 682313188
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  133
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  3   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       27        4       27
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           6   7   9
   3        3          3           5  14  16
   4        3          3           8  11  12
   5        3          2           8  13
   6        3          1          13
   7        3          2          10  16
   8        3          2          15  17
   9        3          2          11  12
  10        3          2          11  12
  11        3          1          14
  12        3          2          13  14
  13        3          2          15  17
  14        3          2          15  17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2  N 3
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0
  2      1     2       5   10    4    8    7
         2     3       3   10    4    6    7
         3     7       3    9    3    6    5
  3      1     5       7    6    7    9   10
         2     8       7    5    4    8   10
         3    10       5    3    4    8    9
  4      1     6       6    7    3    9   10
         2     9       6    5    3    9    9
         3    10       5    3    3    7    8
  5      1     2       6    4    6    2    4
         2     2       5    5    5    2    3
         3     9       3    2    3    2    2
  6      1     1       6    9    8    7    6
         2     7       5    8    8    4    5
         3    10       3    5    7    2    3
  7      1     4       7    9    5    7    8
         2     9       5    7    5    5    7
         3     9       6    9    4    6    7
  8      1     4      10    5   10   10    5
         2     7      10    5    9   10    4
         3     9       9    4    8   10    2
  9      1     2       8   10    9    8    8
         2     2       5    9    8    9   10
         3     3       2    9    5    8    4
 10      1     2       4   10    5    9    3
         2     2       3   10    7    8    3
         3     6       2    8    2    6    1
 11      1     5       7    8    4    7    9
         2     6       6    6    3    7    8
         3     8       6    4    2    7    6
 12      1     4       3    5    4    4    7
         2     5       3    4    4    3    7
         3     9       3    2    4    1    4
 13      1     6       8    3    5    8    9
         2     7       6    2    5    8    8
         3     9       5    1    5    8    8
 14      1     7       6    2    8    7   10
         2     7       5    2    9    8   10
         3     8       2    2    7    4    9
 15      1     5       5    5    5    6    4
         2     5       3    6    5    5    4
         3     8       1    5    3    5    3
 16      1     2       6    9    9    6    5
         2     8       6    9    8    3    5
         3     9       4    8    6    3    4
 17      1     7       9    3    4    9    8
         2     8       5    2    4    7    7
         3     9       3    2    2    6    5
 18      1     0       0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2  N 3
   34   41   99  118  115
************************************************************************
