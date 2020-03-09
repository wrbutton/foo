#!/usr/bin/python

import gt, sys

def summary_v0():

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


def summarize_script(path):
    """ reads one script at a time, tallies references of each function by others in the
    script, and then returns list sorted by fewest to most references, high to low level """

    with open(path, 'r') as f:

        text = f.read()

        blocks = text.split('def ') 

        blocks = [x.strip() for x in blocks]

        # make sure just grabbing from commands, get list of cmds
        cmds = []
        for i,b in enumerate(blocks):
            if i == 0:
              continue
            try: 
                cmd = b.split('(')[0]
                cmds.append(cmd)
            except:
                continue

        # do i need to parse docstrings? should
        # for bl in text:
        #     blocks = []
        #     try:
        #       blocks.append(bl.split('"""')[-1])

        counts = []
        for c in cmds:
            counts.append((c, text.count(c)))

        counts = sorted(counts, key=lambda x: x[1])

        for c in counts:
            print(f'{c[1]} uses of {c[0]}')




            # search through all text, count occurances

        # sort counter object as increasing
        # loop through and print out count, command 

if __name__ == '__main__':
  main()