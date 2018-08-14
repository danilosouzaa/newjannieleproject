************************************************************************
file with basedata            : cr524_.bas
initial value random generator: 161819895
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  132
RESOURCES
  - renewable                 :  5   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       26       14       26
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3          14  15  16
   3        3          2          12  17
   4        3          3           5   6   9
   5        3          2           7   8
   6        3          2          12  13
   7        3          2          10  16
   8        3          2          10  13
   9        3          3          10  13  16
  10        3          3          11  12  14
  11        3          1          17
  12        3          1          15
  13        3          2          14  15
  14        3          1          17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  R 3  R 4  R 5  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0    0    0    0
  2      1     1       7    8    8    7    7    0    7
         2     3       4    8    5    2    6    8    0
         3     3       5    8    4    3    7    0    7
  3      1     1       9    6    2    5    5    7    0
         2     3       6    5    2    3    5    4    0
         3     8       3    2    1    2    5    0   10
  4      1     7       2    7    9    9    7    7    0
         2     7       2    7    7   10    7    7    0
         3     9       2    7    5    9    7    6    0
  5      1     3       8    2    6   10    5    6    0
         2     6       8    2    6    9    5    4    0
         3     9       8    2    3    9    4    0    2
  6      1     4       7    4    1    8    8    6    0
         2     5       7    4    1    7    7    0   10
         3     9       7    3    1    4    3    0    9
  7      1     3      10    9    4    9    9    4    0
         2     4       9    7    3    9    9    0    6
         3     9       9    5    1    8    8    0    3
  8      1     7       9    4    6    3    6    8    0
         2     8       8    4    5    2    3    5    0
         3     9       8    3    5    2    1    5    0
  9      1     1       5   10    4    7    6    7    0
         2     5       3   10    4    5    4    6    0
         3     7       2    9    3    4    4    6    0
 10      1     2       9    5    5    5    5    0    7
         2     4       6    3    4    2    4    0    6
         3     6       5    1    3    2    3    6    0
 11      1     2       8   10    9    4    3    8    0
         2     4       8    9    7    4    3    7    0
         3    10       7    9    2    2    2    5    0
 12      1     2       6   10    8    6    6    0    5
         2     7       6    8    8    6    3    0    4
         3     9       6    7    8    5    3    0    4
 13      1     1       9    8    7    5    9    0    5
         2     1       9    8    6    7    8    7    0
         3     9       7    7    6    1    1    2    0
 14      1     1       5   10   10    7    9    0    8
         2     7       3    7    5    7    9    0    6
         3     9       2    3    3    6    8    8    0
 15      1     2      10    8    8    6   10    9    0
         2     5       8    6    5    5   10    5    0
         3     9       7    6    3    5    9    0    4
 16      1     1       8    7   10    3   10    0    1
         2     3       7    5   10    3    9    5    0
         3    10       7    3    9    3    9    4    0
 17      1     5       5    7    8    2    8    0    9
         2     6       4    7    7    2    8    0    7
         3     7       3    6    7    1    7    0    4
 18      1     0       0    0    0    0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  R 3  R 4  R 5  N 1  N 2
   26   30   27   25   23   78   58
************************************************************************
