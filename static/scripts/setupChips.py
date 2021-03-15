#!/usr/bin/env python
"""This script is part of the initial chips setup: 
1. it generates a {CONDA_ROOT}/envs/chips/etc/conda/activate.d/perl5lib.sh
"""
import os
import sys
from optparse import OptionParser
import subprocess

def main():
    usage = "USAGE: %prog"
    optparser = OptionParser(usage=usage)
    (options, args) = optparser.parse_args(sys.argv)

    #INFER conda root!
    conda_root = ""
    try:
        #check to see if perl5 is installed, i.e. homer is installed 
        #(through chips conda env)
        conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
        #print(conda_root)
    except:
        print("ERROR: conda not installed!  Please install a conda package system\nAnd then create the chips and chips_py2 environments!")
        sys.exit()

    perl5lib=os.path.join(conda_root, "envs", "chips", "lib", "perl5",\
                              "5.18.2")

    if not os.path.isdir(perl5lib):
        print("ERROR: please install the chips environment!\nAND configure the homer pkg")
        sys.exit()
    
    outdir=os.path.join(conda_root,"envs","chips","etc","conda", "activate.d")
    #make outdir:
    subprocess.call(["mkdir","-p",outdir])
    
    #write file:
    f = open(os.path.join(outdir, "perl5lib.sh"), "w")
    f.write("export PERL5LIB=%s\n" % perl5lib)
    f.close()
    
    print("Chips Setup Complete!")

if __name__=='__main__':
    main()
