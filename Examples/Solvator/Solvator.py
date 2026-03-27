import os
import sys
import numpy as np
import argparse

# what is this program doing
# 1. Calculate the dimention of tinker xyz file
# 2. Create a water box that solvates the solute
# 3. Remove the overlapped water molecules
# 4. Add ions

def genWaterBox(xyz, abc):
  #figure out the number of molecules needed 
  mass = 0.0
  massdict = {"H": 1.0079, "O":15.9994}
  navogadro = 6.02214086 
  atoms = np.loadtxt(xyz, usecols=(1,), skiprows=1, dtype="str", unpack=True)
  for atom in atoms:
    mass += massdict[atom.upper()]
  volume = abc[0]*abc[1]*abc[2]
  density = 0.997 # water density
  nmolecule = (int(navogadro*volume/(mass/density)*0.1))

  #use xyzedit to generate the waterbox
  with open("tmp.in", "w") as f:
    f.write("23\n%s\n%s %s %s\nN\n"%(nmolecule, abc[0],abc[1],abc[2]))
  os.system('rm -f water.xyz_?')
  cmdstr = "xyzedit.x water.xyz {prmfile} < tmp.in "
  os.system(cmdstr)
  os.system('mv water.xyz_2 waterbox.xyz')
  os.system(f"sed '2 i\{abc[0]:12.6f}{abc[1]:12.6f}{abc[2]:12.6f}   90.000000   90.000000   90.000000' -i waterbox.xyz")
  return

def solvate():
  if not os.path.isfile('waterbox.xyz'):
    sys.exit('Error: I need a waterbox.xyz to solvate')
  
  x,y,z = np.loadtxt('waterbox.xyz', usecols=(2,3,4), skiprows=2, unpack=True)
  center_water = [(min(x)+max(x))/2.0, (min(y)+max(y))/2.0, (min(z)+max(z))/2.0]
  vector = []
  for i in range(3):
    vector.append(center[i] - center_water[i])
  
  for i in range(0, len(x), 3):
    print(f"checking water {i/3}")
    coord1 = np.array([x[i], y[i], z[i]])
    for j in range(len(coords[0])):
      coord2 = np.array([coords[0][j], coords[1][j], coords[2][j]])
      dist = np.linalg.norm(coord1-coord2)
      if (dist < 3.0):
        print(dist, i, j)
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-xyz', dest = 'xyzfile',  help = "Tinker xyz file (solute)", required=True)  
  parser.add_argument('-prm', dest = 'prmfile',  help = "AMOEBA parameter file", required=True)  
  parser.add_argument('-ion', dest = 'ion',  nargs='+', help = "Ions to add", default = ['Na', '0'])  
  args = vars(parser.parse_args())

  global prmfile 
  xyzfile = args["xyzfile"]
  prmfile = args["prmfile"]
  ion = args["ion"]
  
  # get the dimention of xyz file
  x,y,z = np.loadtxt(xyzfile, usecols=(2,3,4), skiprows=1, unpack=True)
  a = max(x) - min(x)
  b = max(y) - min(y)
  c = max(z) - min(z)
  global coords
  coords = [list(x), list(y), list(z)]

  abc = [a + 20, b + 20, c + 20]
  global center
  center = [min(x)+a/2.0, min(y)+b/2.0, min(z)+c/2.0]
  
  # create an appropriate water box
  genWaterBox('water.xyz', abc)
  
  # solvate the solute in water
  solvate()
  return

if __name__ == '__main__':
  main()