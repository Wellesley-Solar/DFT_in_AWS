import os
import subprocess
import sys
import time

p=8
chem_path="Cs125_gpaw"
home_path=os.path.abspath('.')

os.chdir(chem_path)
dirs=[]
for x in filter(os.path.isdir, os.listdir(os.getcwd())):
    print(x)
    dirs.append(x)

for dir in dirs:
    print(home_path+"/"+chem_path+"/"+dir)
    os.chdir(home_path+"/"+chem_path+"/"+dir)
    child = subprocess.Popen(["mpiexec","-np", str(p), "python3", "gpaw_relax.py"])
    time.sleep(3)
    while child.poll() is None:
        time.sleep(10)