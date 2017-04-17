from sys import argv
from itertools import islice
import os
#4-17-2017
#used to read the last CELL_PARAMETERS in an PWSCF output file
script, pwscf_out = argv # accept input file

pwscf_out_file = open(pwscf_out) # open output file
pwscf_out_txt = pwscf_out_file.read() # read output file to string
l = pwscf_out_txt.rfind('CELL_PARAMETERS') # search for last occurence of string
# print l # check of index
pwscf_out_file.seek(l) # seek the postion in file
str1 = os.path.splitext(pwscf_out)[0]
print str1
print 'CELL_PARAMETERS'
for line in pwscf_out_file:
    if 'CELL_PARAMETERS' in line:
        print '\n'.join(islice(pwscf_out_file, 3))  # read and print out next three lines
pwscf_out_file.close
