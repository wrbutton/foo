#!/usr/bin/python

import gt, sys

def main():

  try:
    path = sys.argv[1]
  except:
    sys.exit('arg1 is code directory to summarize')

  opath = '/Users/WRB/Desktop/codesummary.txt'

  results = []

  flist = gt.get_flist(path, '.py')

  for file in flist:
    fn = file.split('/')[-1]
    if fn[0] == '_':
      continue
    with open(file, 'rU') as df:
      results.append(fn)
      a = 'de'
      b = 'f '
      c = a + b
      for line in df:
        if c in line and 'main' not in line:
          line = line.replace(c, '\t').rstrip().rstrip(':')
          results.append(line)
  with open(opath, 'w') as of:
    for entry in results:
      print(entry, file=of)

if __name__ == '__main__':
  main()