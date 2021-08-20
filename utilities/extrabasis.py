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

def getBasis(element, basis):
  results = []
  lines = open(basis).readlines()[12:]
  for line in lines:
    d = line.split()
    if (len(d) == 2) and (d[1] == "0") and (d[0].upper() == element.upper()):
      idx1 = lines.index(line)
  for i in range(idx1, len(lines)):
    d = lines[i].split()
    results.append(lines[i])
    if "****" in lines[i]:
      break 
  return results
    
def getECP(element, basis):
  results = []
  lines = open(basis).readlines()[12:]
  idx = len(lines)
  for line in lines:
    d = line.split()
    if (len(d) == 3) and (d[0].upper() == (element.upper() + "-ECP")):
      idx = lines.index(line)
  if idx != len(lines):
    results.append(f"{element}    0\n")
  for i in range(idx, len(lines)):
    d = lines[i].split()
    if (len(d) == 2) and (d[1] == "0") :break 
    results.append(lines[i])
  return results
