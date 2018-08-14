************************************************************************
file with basedata            : cr347_.bas
initial value random generator: 1049637680
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  124
RESOURCES
  - renewable                 :  3   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       28        6       28
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2          14  15
   3        3          2           5   9
   4        3          3           6   8   9
   5        3          2           8  10
   6        3          3           7  11  13
   7        3          2          10  14
   8        3          3          12  14  16
   9        3          2          11  16
  10        3          1          12
  11        3          2          12  15
  12        3          1          17
  13        3          3          15  16  17
  14        3          1          17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0
  2      1     5       4    3    4   10    5
         2     5       4    4    3   10    5
         3     9       4    1    1    9    4
  3      1     1       9    9    5    9    6
         2     8       9    9    4    8    5
         3     9       8    6    2    8    5
  4      1     1       5    5    5    2    9
         2     4       4    5    4    1    9
         3     5       4    5    3    1    9
  5      1     2       2    7    6    9    8
         2     5       2    7    5    6    3
         3     7       1    4    4    2    3
  6      1     7       5    4    8    9    8
         2     8       2    3    5    8    7
         3     9       2    3    1    8    5
  7      1     5       5    8    7    6    8
         2     8       5    7    6    4    7
         3     9       5    7    5    2    5
  8      1     1       5    5    3    5    8
         2     1       8    4    3    6    9
         3     4       3    2    2    5    8
  9      1     7       6    5    3    7    2
         2     7       6    5    3    1    7
         3     7       6    7    2    5    3
 10      1     8      10    7    5    4    9
         2     8       8    7    5    5    9
         3    10       8    6    4    3    3
 11      1     3      10   10    7    5    4
         2     6       7   10    6    5    4
         3     9       3    9    4    5    3
 12      1     4       6    7    9    9   10
         2     9       4    7    8    6   10
         3    10       3    7    7    3    9
 13      1     1       8    6    4    1    7
         2     5       8    6    3    1    6
         3    10       8    5    1    1    6
 14      1     1       9    2    2    6    7
         2     2       7    2    2    6    6
         3     6       5    1    1    4    6
 15      1     1       7    8    2    3    2
         2     3       5    6    2    3    2
         3     8       4    6    2    2    2
 16      1     3       7    5    5    6    8
         2     5       5    5    3    4    8
         3     6       2    3    2    3    4
 17      1     3       5    9    3    9    5
         2     3       4    8    3    7    6
         3     6       3    7    1    3    4
 18      1     0       0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  N 1  N 2
   19   20   18   81   96
************************************************************************
