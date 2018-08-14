************************************************************************
file with basedata            : c2143_.bas
initial value random generator: 1306929768
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  138
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       23        6       23
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6  14
   3        3          3          11  13  14
   4        3          3           6   7   8
   5        3          3           9  10  12
   6        3          3           9  10  11
   7        3          3           9  12  14
   8        3          3          10  11  12
   9        3          1          17
  10        3          2          13  15
  11        3          2          15  16
  12        3          2          13  17
  13        3          1          16
  14        3          3          15  16  17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       6    0    7    9
         2     7       5    0    6    6
         3    10       0    5    6    5
  3      1     1       7    0    4    4
         2     4       0    7    3    3
         3     7       5    0    1    3
  4      1     8       8    0    7    5
         2     8       9    0    7    4
         3     9       4    0    6    4
  5      1     7       0    4    5    4
         2     8       2    0    4    2
         3    10       0    3    4    1
  6      1     1       0    4    9    5
         2     2      10    0    9    4
         3    10       0    3    6    3
  7      1     7       5    0   10    7
         2    10       4    0    8    5
         3    10       2    0    9    7
  8      1     2       0    4    9   10
         2     2       0    3   10    8
         3     9       0    2    7    7
  9      1     1       8    0    8    1
         2     8       0    8    6    1
         3    10       3    0    6    1
 10      1     2       7    0    8    7
         2     5       0    6    7    7
         3     7       4    0    7    6
 11      1     2       5    0    8    9
         2     7       0    4    8    7
         3     8       5    0    7    5
 12      1     2       0    5    7    7
         2     3       0    1    7    7
         3     6       4    0    5    3
 13      1     1       0   10    5    8
         2     5       0   10    5    7
         3     9       5    0    3    7
 14      1     4       5    0   10    3
         2     6       4    0    9    3
         3     9       4    0    9    2
 15      1     1       0    8    8    9
         2     6       0    7    5    8
         3     8       0    7    2    8
 16      1     4       8    0    6    4
         2     6       5    0    5    4
         3     8       3    0    3    3
 17      1     3       0    2    9    8
         2     4       8    0    6    6
         3     8       8    0    3    5
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   18   21  102   84
************************************************************************
