************************************************************************
file with basedata            : cn361_.bas
initial value random generator: 108105199
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  126
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  3   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       23       11       23
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5  10  15
   3        3          3           6   8   9
   4        3          2           6   7
   5        3          3          11  13  14
   6        3          2          11  14
   7        3          3          12  15  16
   8        3          2          13  16
   9        3          2          11  15
  10        3          1          12
  11        3          2          16  17
  12        3          2          13  14
  13        3          1          17
  14        3          1          17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2  N 3
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0
  2      1     2       9    6    6    7    8
         2     4       7    5    3    7    7
         3     6       5    4    2    7    6
  3      1     2       5   10    7    2    5
         2     3       4    9    6    2    4
         3     9       3    9    6    2    4
  4      1     8       4   10    9    3    2
         2     9       4    9    7    3    1
         3    10       3    9    6    3    1
  5      1     4       8    7    8    6    7
         2     4       8    7    9    6    6
         3     8       8    6    7    4    4
  6      1     1       6    6   10    3    8
         2     4       5    5    8    3    5
         3     5       3    4    7    2    5
  7      1     3       5    3    5    6   10
         2     6       5    3    4    4    8
         3     8       5    2    4    3    7
  8      1     2       8    5    9    6    7
         2     7       7    4    9    6    5
         3     8       5    3    9    3    2
  9      1     2       6    8    7    8    7
         2     3       6    7    7    4    7
         3     4       6    7    7    2    5
 10      1     3       6    9    9    8    3
         2     8       4    8    5    7    3
         3     9       2    6    2    6    2
 11      1     2       7    4    7   10    2
         2     5       6    4    4    9    1
         3     9       6    3    3    9    1
 12      1     2       7   10    6    8    8
         2     3       3    8    4    7    5
         3     9       1    7    3    7    3
 13      1     5       7    6    8    4    6
         2     6       6    3    7    4    5
         3     7       5    3    6    4    4
 14      1     1      10    8    9    2    8
         2     5       7    8    8    1    8
         3     6       6    7    6    1    5
 15      1     6       6    9    8    6    7
         2     7       6    7    8    5    7
         3    10       5    6    8    1    4
 16      1     6       5    4    8    8    9
         2     9       3    3    8    7    6
         3     9       4    2    8    8    5
 17      1     5       5    2    6    7    4
         2     7       5    2    5    7    3
         3     9       2    2    5    4    3
 18      1     0       0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2  N 3
   14   17  123   94  101
************************************************************************