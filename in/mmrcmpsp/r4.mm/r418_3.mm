************************************************************************
file with basedata            : cr418_.bas
initial value random generator: 1042619514
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  121
RESOURCES
  - renewable                 :  4   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       22       12       22
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           6   7  12
   3        3          3           5   8   9
   4        3          3           5  12  13
   5        3          1          10
   6        3          2          11  17
   7        3          2           9  10
   8        3          3          10  12  13
   9        3          2          16  17
  10        3          1          11
  11        3          2          14  16
  12        3          2          15  16
  13        3          2          14  17
  14        3          1          15
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  R 4  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0    0
  2      1     2       8    9    2    0    5    0
         2     6       5    0    0    0    5    0
         3    10       0    5    0    0    0    3
  3      1     1       0    4    0    7    0    6
         2     4       0    2    0    0    7    0
         3     5       0    0    0    6    0    3
  4      1     2       6    5    7    0    0    9
         2     2       0    0    6    0    6    0
         3     7       0    6    5    9    0   10
  5      1     2       9    4    0    0    0    4
         2     6       9    3    0    0    8    0
         3     8       8    3    0    4    7    0
  6      1     2       0    5    6    0    0    5
         2    10       5    4    5    2    0    2
         3    10       8    2    0    2    8    0
  7      1     3       0    0    5    0    7    0
         2     7       0   10    0    0    0    5
         3     9       8    0    3    0    0    2
  8      1     3       0    0    7    0    0   10
         2     8       4    0    5    0    0    3
         3     8       4    0    0    0    0    4
  9      1     2       0   10    0    0    8    0
         2     3       2    9    7    0    0    6
         3     4       0    0    3    0    0    5
 10      1     4       0    8    0    0    0    7
         2     9       0    8    0    0    0    6
         3     9       3    0    0    0    0    6
 11      1     6       0    8    0    5    0    6
         2     8       0    8    0    0    8    0
         3    10       3    0    0    0    5    0
 12      1     1       0    7    7    0    0    6
         2     5       9    6    3    7    0    5
         3     7       8    0    3    4    0    5
 13      1     3       8    0    0    6    0    3
         2     5       0    3    0    6    8    0
         3     7       3    0    7    0    6    0
 14      1     3       0    5   10    9    3    0
         2     4       8    4    8    8    2    0
         3     5       0    0    0    6    0    6
 15      1     4       6    0    8    0    5    0
         2     5       3    8    8    0    0    9
         3     8       0    0    0    5    5    0
 16      1     1       8    9    7    0    6    0
         2     3       8    8    7    9    2    0
         3     4       6    6    0    5    0   10
 17      1     7       0    2    0    0    6    0
         2     8       1    0    6    0    0    7
         3    10       0    1    0    4    5    0
 18      1     0       0    0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  R 4  N 1  N 2
   21   15   14   13   64   81
************************************************************************
