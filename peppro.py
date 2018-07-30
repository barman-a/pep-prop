# This program is to calculate physio-chemical properties of peptide
# Given a sequence this program will calculate the isoelectric point (PI)
# Charge, Molecular weight, Molar extiction coeff,  Aliphatic index , and average hydropathy 
#
# Usage: python peppro.py sequence     (One liner sequence file sequence.txt)

from __future__ import division
import os
import string
import sys
import math
from sys import argv

try:
   input = sys.argv[1]
except:
       print ("no Input")

seq = open(input+'.txt', "r").readline()
#seq = list(seq)

residues = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
mol_wt = ['71.07', '156.18', '114.10', '115.08', '103.14', '128.13', '129.11', '57.05', '137.14', '113.15', '113.15', '128.17', '131.19', '147.17', '97.11', '87.07', '101.10', '186.21', '163.17', '99.13'] ##FASG760101 
termini = 1.00 + 17.00

hydropathy = ['1.800', '-4.500', '-3.500', '-3.500', '2.500', '-3.500', '-3.500', '-0.400', '-3.200', '4.500', '3.800', '-3.900', '1.900', '2.800', '-1.600', '-0.800', '-0.700', '-0.900', '-1.300', '4.200']  ##kyte-doolittle

seq_old = list(seq)

molwt = 0.0
hydro = 0.0
asp = 0
glu = 0
arg = 0
lys = 0
tyr = 0
cys = 0
his = 0
ntr = 0
ctr = 0

ala = 0
val = 0
leu = 0
ile = 0
trp = 0

for i in range(0,len(seq)):
    res = seq[i]
    if res == 'D':
       asp = asp + 1
    if res == 'E':
       glu = glu + 1
    if res == 'R':
       arg = arg + 1
    if res == 'K':
       lys = lys + 1
    if res == 'Y':
       tyr = tyr + 1
    if res == 'C':
       cys = cys + 1
    if res == 'H':
       his = his + 1
    if res == 'A':
       ala = ala + 1
    if res == 'V':
       val = val + 1
    if res == 'L':
       leu = leu + 1
    if res == 'I':
       ile = ile + 1
    if res == 'W':
       trp = trp + 1

    for j in range(0,len(residues)):
        if res == residues[j]:
           molwt = molwt + float(mol_wt[j])
           hydro = hydro + float(hydropathy[j])

totmolwt = molwt+termini
avg_hydro = hydro / i

qasp = 0.0
qglu = 0.0
qarg = 0.0
qlys = 0.0
qtyr = 0.0
qcys = 0.0
qhis = 0.0
qntr = 0.0
qctr = 0.0

qtot = 0.0

pH = 0.0

while True:
    qntr = -1/(1+pow(10,(3.65-pH)))   
    qarg =  arg/(1+pow(10,(pH-12.48)))
    qlys =  lys/(1+pow(10,(pH-10.54)))
    qhis =  his/(1+pow(10,(pH-6.04)))                                      
    qasp = -asp/(1+pow(10,(3.9-pH)))           
    qglu = -glu/(1+pow(10,(4.07-pH)))           
    qcys = -cys/(1+pow(10,(8.18-pH)))           
    qtyr = -tyr/(1+pow(10,(10.46-pH)))        
    qctr = 1/(1+pow(10,(pH-8.2)))                

    qtot = qntr+qarg+qlys+qhis+qasp+qglu+qcys+qtyr+qctr

    if pH >= 14.0:
       print "pH range more than 14.0"
       break
    if qtot < 0.0:
       break

    pH = pH + 0.01

charge = arg + lys - asp - glu

extinction_coeff = (tyr * 1490) + (trp * 5500) 

ala_ndx = (ala/i) * 100
val_ndx = (val/i) * 100
ile_ndx = (ile/i) * 100
leu_ndx = (leu/i) * 100
aliphatic_ndx = ala_ndx + (val_ndx * 2.9) + (leu_ndx * 3.9) + (ile_ndx * 3.9)


print "  PI  Charge  MolWt   ExtCoeff   AliIndex   Hydropathy"
print "===================================================="
print pH, charge, totmolwt, extinction_coeff, round(aliphatic_ndx,2), round(avg_hydro,2)  
 
