#!/usr/bin/python
import os, sys, re

setInstance = ['A-', 'B-', 'X-']
setNI = [1,2,3,4,5,6,7,8,9,10]
for inst in setInstance:
   for i in setNI:
     instance = inst+str(i)  
     os.system('java -jar MRCMPSP-validator.jar ' +instance+'.txt ' +instance+'.sol' )
