#! /usr/bin/env python3

import sys, os, re
from mycotools.lib.kontools import sysStart, formatPath, fmt_float
from mycotools.lib.biotools import fa2dict, calc_weight


if __name__ == '__main__':

    usage = 'USAGE: Inputs amino acid fasta outputs linear protein weights'
    args = sysStart(sys.argv, usage, 1)

    fa = fa2dict(formatPath(sys.argv[1]))
    print('#protein\tkDa', flush = True)
    for head in fa:
        print(head + '\t' + fmt_float(calc_weight(fa[head]['sequence'])), flush = True)
