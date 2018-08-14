************************************************************************
file with basedata            : mf38_.bas
initial value random generator: 11595300
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  32
horizon                       :  254
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     30      0       24       11       24
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   6   8
   3        3          3           9  15  29
   4        3          3           7  11  13
   5        3          1          23
   6        3          3           7  10  27
   7        3          3          14  17  19
   8        3          3          10  12  16
   9        3          1          20
  10        3          2          13  22
  11        3          2          22  25
  12        3          3          21  22  27
  13        3          1          14
  14        3          2          21  24
  15        3          3          24  25  28
  16        3          2          19  21
  17        3          3          18  20  24
  18        3          1          26
  19        3          2          23  26
  20        3          1          28
  21        3          2          26  30
  22        3          1          28
  23        3          1          31
  24        3          1          30
  25        3          1          30
  26        3          2          29  31
  27        3          1          29
  28        3          1          31
  29        3          1          32
  30        3          1          32
  31        3          1          32
  32        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     1       8    6    5    2
         2     6       6    5    5    2
         3     7       6    5    4    2
  3      1     5       8    8    8    7
         2     7       5    7    6    6
         3    10       3    7    5    5
  4      1     5       5    6    7    8
         2     7       4    3    6    7
         3     8       1    1    6    6
  5      1     2       4    6    5    9
         2     7       4    5    3    7
         3     8       3    3    2    6
  6      1     3       8    7    3    8
         2     3       7    8    4    9
         3    10       7    3    3    6
  7      1     3       7    8    5    4
         2     4       6    5    4    4
         3     6       5    1    1    3
  8      1     2       5    5    7    1
         2     3       4    4    6    1
         3    10       3    2    6    1
  9      1     4       6    8    8    8
         2     6       6    6    6    7
         3    10       4    6    4    7
 10      1     4       3   10    7    8
         2     5       3    9    5    7
         3    10       2    8    4    7
 11      1     4       7    3    9   10
         2     7       7    3    8   10
         3    10       4    2    6   10
 12      1     2       3    3    4    1
         2     6       2    3    3    1
         3     9       1    2    3    1
 13      1     2       6    4    8    2
         2     3       5    3    6    2
         3     9       5    3    5    2
 14      1     1       7    5    3    6
         2    10       5    5    3    5
         3    10       6    4    3    5
 15      1     1       8    9    8    8
         2     4       6    9    7    8
         3     6       6    8    7    7
 16      1     3       8    4    6    9
         2     7       6    3    6    9
         3     9       4    2    6    8
 17      1     2       5    9    3    7
         2     7       4    8    3    6
         3     9       3    5    2    5
 18      1     4       8    8    9    7
         2     5       8    7    8    4
         3     6       7    7    7    3
 19      1     7       5    8    5    2
         2     9       3    7    3    1
         3     9       3    4    5    2
 20      1     2       9    7    9   10
         2     8       8    7    7    8
         3    10       7    3    6    7
 21      1     8       5    8    5    9
         2     8       8    7    5    9
         3    10       2    6    3    8
 22      1     4       5    5    4   10
         2     8       5    4    4    9
         3    10       2    4    3    9
 23      1     2       8   10    6   10
         2     6       7    6    2    9
         3     6       7    4    6   10
 24      1     2       9    8    6    6
         2     4       8    7    3    5
         3    10       5    6    2    4
 25      1     2       7    8    2    4
         2     3       7    7    2    3
         3     8       7    6    1    3
 26      1     1       6    6    6    7
         2     3       4    6    3    4
         3     3       3    6    5    4
 27      1     2       8    7    6    7
         2     3       8    5    6    6
         3     6       7    3    6    4
 28      1     1       7    9    4    6
         2     4       5    8    4    6
         3     8       5    6    3    6
 29      1     4       7    5   10    8
         2     7       6    4    7    7
         3    10       6    3    6    4
 30      1     3       5    4    6    5
         2     3       7    6    5    6
         3    10       4    1    4    5
 31      1     2       1    8    6   10
         2     7       1    6    4    4
         3     7       1    7    2    3
 32      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   27   29  134  164
************************************************************************