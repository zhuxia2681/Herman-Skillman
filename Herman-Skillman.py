#!/usr/bin/env python3

"""

Jun INAGAKI
2021.10.7

"""
 
import numpy as np
import sys,os,re

# start
print("Herman-Skillman code")


#initialize
x = [0]*521
r = [0]*521
nnlz = [0]*24
wwnl = [0]*24
ee = [0.0]*24
ru  = [0.0]*521
ru2 = [0.0]*521
ru3 = [0.0]*521
v = [0.0]*521

nfiles = 0

# read heading card
f = open("in","r")
lines = f.readlines()

print(lines[0])


# read control cards and input potentials. calculate trial potential
tmp = lines[2].split(",")

key = int(tmp[0])
tol = float(tmp[1])
thresh = float(tmp[2])
mesh = int(tmp[3])
ipratt = int(tmp[4])
maxit = int(tmp[5])
nocopy = int(tmp[6])
kut = int(tmp[7])

#
LimMesh = int(lines[3])
print(LimMesh)

ncards = 90
print(" key = ",key)

if maxit <= 0:
    maxit = 20

nblock = int(mesh / 40)
print(" nblock = ",nblock)


# construct x mesh and r mesh

i = 0
x[i] = 0.0
r[i] = 0.0
deltax = 0.0025
for j in range(nblock):
    for k in range(40):
        i += 1
        x[i] = x[i-1] + deltax
    deltax += deltax

ru2 = [0.0]*521
if key == 1:
    print("read ru2.")
    print("read ru3.")
    ze2 = -1 * ru2[0] * 0.5
    ze3 = -1 * ru3[0] * 0.5
# read in atomic potential
elif key == -1:
    print("read in atomic potential.")
#    m = 0
#    for i in range(1,437,4):
#        ru2[m] = 
elif key == 0:
    n = 4
    m = 0
    for j in range(11):
        tmp = lines[n].split()
        for i in range(10):
            ru2[m] = float(tmp[i])
            m += 4
        n += 1
        
# read
tmp = lines[n].split(",")
z = float(tmp[0])
ncores = int(tmp[1])
nvales = int(tmp[2])
ion = int(tmp[3])

# z
if z <= 0:
    print("Z must over zero.")
    sys.exit()

nfiles += 1
iz = z
ncspvs = ncores + nvales
c = 0.88534138 /z**(1. /3. )
twoion = ion + ion
zzz = ion + 1

#

twozzz = zzz + zzz
m = 1
for i in range(mesh):
    r[m] = c * x[m]
    m += 1

n += 1
m = 0
for i in range(ncspvs):
    tmp = lines[n].split(",")
    nnlz[m] = int(tmp[0])
    wwnl[m] = float(tmp[1])
    ee[m] = float(tmp[2])
    m += 1
    n += 1

www = 0.0

n = 0
for i in range(ncspvs):
    www = www + wwnl[n]
    n += 1

chk = np.abs(z + 1.0 - www - zzz) - 0.001
if chk >= 0:
    sys.exit()

# assume key = 0.

twoz = z + z
for i in range(0,436,4):
    ru[i] = -1 * ru2[i] * twoz
ru[440] = ru[436]
ru[444] = ru[436]
m = 9
for i in range(0,436,4):
    m -= 1
    if m < 0:
        ru[i]   = (22.0*ru[i-1] + 11.0*ru[i+3] - ru[i+7]) / 32.0
        ru[i+1] = (10.0*ru[i-1] + 15.0*ru[i+3] - ru[i+7]) / 24.0
        ru[i+2] = ( 6.0*ru[i-1] + 27.0*ru[i+3] - ru[i+7]) / 32.0
    if m >= 0:
        ru[i]   = (21.0*ru[i-1] + 14.0*ru[i+3] - 3.0*ru[i+7]) / 32.0
        ru[i+1] = ( 3.0*ru[i-1] +  6.0*ru[i+3] -     ru[i+7]) /  8.0
        ru[i+2] = ( 5.0*ru[i-1] + 30.0*ru[i+3] - 3.0*ru[i+7]) / 32.0

# so skip to 31.

v[0] = -9.9E35

m = min(440,mesh)
print(m)




