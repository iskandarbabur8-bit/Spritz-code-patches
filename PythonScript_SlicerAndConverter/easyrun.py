#!/usr/bin/env python
import linecache
import sys
import os
import re
import subprocess
import shutil
from H5_entropyslice     import entropyslicer
from H5_temperatureslice import temperatureslicer

tablepath  = os.path.abspath(sys.argv[1])
Lorenepath = os.path.abspath(sys.argv[2])

print("You have selected EoS in",tablepath)
print("Temperature slice or Entropy slice (T/S)?",end='')
choice = input()
if choice not in ["T","S"]:
    print(choice,"is non acceptable. Aborting.")
    exit()

if (choice == "T"):
    print("Chosen: Tempslice")
    slicedfilepath = temperatureslicer()
elif (choice == "S"):
    print("Chosen: Entropyslice")
    slicedfilepath = entropyslicer()

print("\n********* Now generating files for executing Lorene/nrotstar... **********\n")
nrotstarpath = os.path.join(Lorenepath,"Codes/Nrotstar")

def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)

slicedfilepath = os.path.abspath(slicedfilepath)

os.chdir(nrotstarpath)
foldername = os.path.abspath(input("Write name of folder you want in Lorene/Nrotstar: "))
if not os.path.exists(foldername):
    copyDirectory("LS_LS220_slice0",foldername)
else:
    print("WARNING: folder already exists.")
input("Press ENTER for running nrotstar (overwriting if folder exists)...")

#change EoS in par_eos.d
os.chdir(foldername)
with open('par_eos.d', 'r') as f:
    lines = f.readlines()
    lines[-1] = os.path.join(slicedfilepath,"eos")
with open('par_eos.d', 'w') as f:
    f.writelines(lines)

subprocess.call(["../nrotstar"])

print("\n\n******** Parameters that you have to add in you parfile are: *********\n\n")
print("ID_NrotStar::filename  =  ", os.path.join(foldername,"resu.d"))
print("Spritz_SetBeta::beta_datafile  =  ", os.path.join(slicedfilepath,"eos"))
print("EOS_Omni::nuceos_table_name  =  ", tablepath)
print("********************************************************************")

# os.chdir(directory giusta)
#moifica par.eos
#esegui ../nrotstar
#stampa a schermo ciò che se deve mettere nel parfile.
