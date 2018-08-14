************************************************************************
file with basedata            : cr446_.bas
initial value random generator: 1428652377
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  130
RESOURCES
  - renewable                 :  4   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       32       10       32
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2          15  17
   3        3          3           5  11  13
   4        3          3           6   8  13
   5        3          2           7  12
   6        3          3          10  11  12
   7        3          2           9  15
   8        3          3          10  11  12
   9        3          1          10
  10        3          2          14  17
  11        3          1          16
  12        3          3          14  15  17
  13        3          1          14
  14        3          1          16
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  R 4  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0    0
  2      1     5       9   10    7    5    9   10
         2     6       8   10    6    5    9    9
         3     9       7    9    6    4    8    9
  3      1     7       8    7    7    3    8    5
         2     9       6    7    4    3    7    3
         3    10       6    6    1    2    5    3
  4      1     1       8    5    9    2    4    8
         2     4       6    5    9    1    3    7
         3     9       4    4    8    1    3    3
  5      1     6       3    6    7    5    3    6
         2     8       2    5    7    5    3    4
         3    10       2    4    7    3    2    3
  6      1     5       9    3   10    7    6    7
         2     7       6    3    9    6    4    7
         3     9       6    3    9    2    4    6
  7      1     5       4    7    8    8    9    9
         2     9       3    6    7    7    6    7
         3    10       2    4    2    7    6    7
  8      1     2       6   10    5    3    9    4
         2     3       4    9    3    2    4    3
         3     4       2    9    3    2    4    2
  9      1     4       2    3    5    8    4    8
         2     8       2    2    3    3    4    6
         3     8       2    3    3    4    3    6
 10      1     2       3   10    8    5   10    4
         2     4       3    5    7    5    9    3
         3    10       3    3    6    3    7    3
 11      1     4       6    7   10    8   10    5
         2     5       3    7    8    7    8    1
         3     5       4    7    8    6    7    3
 12      1     1       9    7    4    8    6    3
         2     2       9    7    2    5    5    3
         3     3       8    6    2    3    5    3
 13      1     7       5    3    3   10    2    8
         2     9       4    2    3    6    2    6
         3     9       4    3    2    7    1    8
 14      1     6       2    5    8    6    8    4
         2     7       1    5    4    3    8    1
         3     7       2    3    5    4    8    4
 15      1     6       6    4   10    6    6    2
         2     8       6    4    7    3    4    2
         3     9       3    1    5    2    4    1
 16      1     2       9   10   10    7    4    1
         2     6       9   10    9    7    3    1
         3     8       8    9    9    6    2    1
 17      1     6       7    7    7    7    7    2
         2     7       7    6    6    7    5    2
         3    10       6    4    6    7    5    1
 18      1     0       0    0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  R 4  N 1  N 2
   20   20   19   17   90   71
************************************************************************
