#! /usr/bin/env python3

import os
import re
import sys
from mycotools.lib.kontools import sys_start, format_path, fmt_float
from mycotools.lib.biotools import fa2dict, calc_weight


def cli():

    usage = 'USAGE: Inputs amino acid fasta outputs linear protein weights'
    args = sys_start(sys.argv, usage, 1)

    fa = fa2dict(format_path(sys.argv[1]))
    print('#protein\tkDa', flush = True)
    for head in fa:
        print(head + '\t' + fmt_float(calc_weight(fa[head]['sequence'])), flush = True)


if __name__ == '__main__':
    cli()
