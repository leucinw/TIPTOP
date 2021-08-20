''' write template xyz file from amoebabio18.prm for amino acid '''

''' usage: python lWriteTemplate.py '''

import os

def getXYZ(residue):
  cmd = f"sh protein.sh {residue}"
  os.system(cmd)
  tripep = residue + "3.xyz"
  lines = open(tripep).readlines()
  natom = int(lines[0].split()[0])
  n3 = int(natom/3)
  n1 = n3 + 1
  n2 = n3 - 1
  f1 = f"{residue}N.txyz"
  f2 = f"{residue}.txyz"
  f3 = f"{residue}C.txyz"
  with open(f1, "w") as f:
    f.write(f"{n1:>6d} Template\n")
    for i in range(1, n1+1):
      dd = lines[i].split()
      conlist = []
      for n in dd[6:]:
        if int(n) <= n1:
          conlist.append(n)
      linestr = f"{dd[0]:>6s}{dd[1][0]:>3s}{float(dd[2]):12.6f}{float(dd[3]):12.6f}{float(dd[4]):12.6f}{dd[5]:>6s}" + '   ' + '   '.join(conlist)
      f.write(f"{linestr}\n")
  with open(f2, "w") as f:
    f.write(f"{n2:>6d} Template\n")
    for i in range(n1+1, n2+n1+1):
      dd = lines[i].split()
      conlist = []
      for n in dd[6:]:
        if n1 < int(n) <= n1+n2:
          conlist.append(str(int(n)-n1))
      linestr = f"{int(dd[0])-n1:>6d}{dd[1][0]:>3s}{float(dd[2]):12.6f}{float(dd[3]):12.6f}{float(dd[4]):12.6f}{dd[5]:>6s}" + '   ' + '   '.join(conlist)
      f.write(f"{linestr}\n")
  with open(f3, "w") as f:
    f.write(f"{n3:>6d} Template\n")
    for i in range(n1+n2+1, n2+n3+n1+1):
      dd = lines[i].split()
      conlist = []
      for n in dd[6:]:
        if n1+n2 < int(n) <= n1+n2+n3:
          conlist.append(str(int(n)-n1-n2))
      linestr = f"{int(dd[0])-n1-n2:>6d}{dd[1][0]:>3s}{float(dd[2]):12.6f}{float(dd[3]):12.6f}{float(dd[4]):12.6f}{dd[5]:>6s}" + '   ' + '   '.join(conlist)
      f.write(f"{linestr}\n")
  os.system(f"rm {tripep}")     
  return


def main():

  ## amino acids in amoebabio18.prm
  residues = \
  [ 'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'CYS', 'LYD', 'PRO', \
    'PHE', 'TYR', 'TRP', 'HIS', 'ASP', 'ASN', 'GLU', 'GLN', 'MET', 'LYS', \
    'ARG', 'ORN', 'TYD', 'HID', 'HIE', 'HIP', 'ASH', 'GLH', ] # leave CYX, AIB, CYD out for now
  for residue in residues[0:]:
    getXYZ(residue)
  return

if __name__ == "__main__":
  main()
