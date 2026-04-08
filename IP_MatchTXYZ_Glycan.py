''' 
MIT License

Copyright (c) 2021 Chengwen Liu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import shutil
import subprocess
import numpy as np


def _run(cmd, check=True):
  ''' Run a command using subprocess; cmd is a list of arguments. '''
  result = subprocess.run(cmd, capture_output=True, text=True)
  if check and result.returncode != 0:
    sys.exit(f"Command failed: {' '.join(str(c) for c in cmd)}\n{result.stderr}")
  return result


if __name__ == "__main__":
  templates = sys.argv[1:-1]
  fname = sys.argv[-1]
  rootdir = os.path.join(os.path.split(__file__)[0])

  if shutil.which("obabel") is None:
    sys.exit("Error: 'obabel' not found on PATH. Please install Open Babel first.")

  if not os.path.isfile(f"{fname}.pdb"):
    sys.exit(f"Error: PDB file not found: {fname}.pdb")

  _run(["obabel", "-ipdb", f"{fname}.pdb", "-oxyz", "-O", f"{fname}_H.xyz", "-h"])
  _run(["obabel", "-ixyz", f"{fname}.xyz", "-otxyz", "-O", f"{fname}.txyz"])
  _run([sys.executable, os.path.join(rootdir, "IP_MatchTXYZ.py"), "-t"] + templates + ["-d", f"{fname}_H.xyz"], check=False)

  # here the idea is to add H back and run IP_MatchTXYZ
  # added Hs are always appended in the end
  # assign types for the needed fragment
  txyz_file = f"{fname}.txyz"
  matched_file = f"{fname}_H.txyz_2"

  if not os.path.isfile(txyz_file):
    sys.exit(f"Error: expected file not found: {txyz_file}")
  if not os.path.isfile(matched_file):
    sys.exit(f"Error: expected file not found: {matched_file}")

  lines = open(txyz_file).readlines()
  types = np.loadtxt(matched_file, usecols=(5,), dtype='str', unpack=True, skiprows=1)[:len(lines)]

  with open(fname + ".txyz_2", 'w') as f:
    f.write(lines[0])
    for i in range(1, len(lines)):
      d = lines[i].split()
      line = ' '.join(d[:5] + [types[i-1]] + d[6:]) + "\n"
      f.write(line)
