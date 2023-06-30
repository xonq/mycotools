#! /usr/bin/env python3

import os
import re
import sys
import argparse
from mycotools.lib.kontools import intro, outro, file2list, format_path, mkOutput

def grab_names(data, query = False):

    if query:
        comp = re.compile(r'Query\:       (.*?)\W+\[M\=\d+\]')
    else:
        comp = re.compile(r'Accession\:   (.*)')

    namesPrep = comp.finditer(data)
#    namesPrep = comp.findall(data)
    names_set = set(x[1] for x in namesPrep)

    return names_set


def sort_hits(hit_list, aln_list):
    hit_dict = {i: v for i, v in enumerate(hit_list)}
    # sort by bitscore because hmmer doesnt have identity :/
    sortd_hits = {i: v for i, v in sorted(hit_dict.items(), key = lambda x: \
                                         float(x[1][2]), reverse = True)}
    hit_list = list(sortd_hits.values())
    aln_list = [aln_list[i] for i in sortd_hits]
    return hit_list, aln_list
                                        


def grab_hits(data, threshold = None, evalue = None, bitscore = None,
             top = None, query = False, accession = False):

    if query:
        comp = re.compile( r'^Query: +' + query + r' +\[M\=(\d+)\]$\n([\s\S]+?\n\n)' \
            + r'([\s\S]*?\n\n\n)', re.M )
    else:
        comp = re.compile( r'Query\:.*?\[M\=(\d+)\]$\n(^Accession\: +' + accession \
            + r'[\s\S]+?\n\n)([\s\S]*?\n\n\n)', re.M )

    hitsAligns = comp.search( data )
    if hitsAligns is not None:
        querySize = int(hitsAligns[1])
        fullHit = hitsAligns[2]
        fullAlign = hitsAligns[3]
        lines = fullHit.split('\n')
        prepLines = []
        for line in lines:
            if re.search(r'^\W+\d', line):
                prepLines.append(line)

        outputLines = []
        for line in prepLines:
            outPrep = line.split(' ')
            outLine = [x for x in outPrep if x != '']
            outputLines.append( outLine )

        hit_str = ''
        for line in outputLines:
            hit_str += f"{' '.join(line[8:])}\t{line[0]}\t{line[1]}\t{line[2]}\t" \
                    + f"{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\n"

        passingSeq, aln_str, count = set(), '', 0
        aligns = re.findall( r'^\>\> (.*?)\n.*?\n[\W\-]+\n(.*?)\n\n', fullAlign, re.DOTALL | re.MULTILINE )
        for align in aligns:
            seq = align[0].rstrip()
            lines = align[1].split('\n')
            for line in lines:
                outPrep = line.split(' ')
                outLine = [x for x in outPrep if x != '' and x != '..' and '[' not in x and ']' not in x]
                if len(outLine) != 13:
                    print('\nINVALID ALIGNMENT HEADERS - CHECK SCRIPT/ALIGNMENT\n', flush = True)
                    sys.exit( 5 )
                if threshold:
                    hmmStart = int(outLine[6])
                    hmmEnd = int(outLine[7])
                    hitSize = hmmEnd - hmmStart + 1
                    if hitSize < threshold * querySize:
                        continue
                    passingSeq.add( str(seq) )
                    

                aln_str += str(seq) + '\t' + outLine[0] + '\t' + outLine[1] + '\t' + outLine[2] + \
                    '\t' + outLine[3] + '\t' + outLine[4] + '\t' + outLine[5] + '\t' + outLine[6] + '\t' + \
                    outLine[7] + '\t' + outLine[8] + '\t' + outLine[9] + '\t' + outLine[10] + '\t' + \
                    outLine[11] + '\t' + outLine[12] + '\n'

        t_hit_list = [hit.split('\t') for hit in hit_str.split('\n') if hit]
        if threshold:
            hit_list = []
            for hit in t_hit_list:
                if hit[0] in passingSeq:
                    hit_list.append(hit)
        else:
            hit_list = t_hit_list
        aln_list = [aln.split('\t') for aln in aln_str.split('\n') if aln]
        hit_list, aln_list = sort_hits(hit_list, aln_list)

        if evalue:
            new_passSeq = set()
            out_hits, out_aligns = [], []
            for hit in hit_list:
                if int(10**200 * float(hit[1])) < int(10**200 * float(evalue)):
                    out_hits.append(hit)
                    new_passSeq.add(hit[0])
            for align in aln_list:
                if align[0] in new_passSeq:
                    out_aligns.append(align)
            hit_list, aln_list = out_hits, out_aligns
        if bitscore:
            new_passSeq = set()
            out_hits, out_aligns = [], []
            for hit in hit_list:
                if float(hit[2]) >= bitscore:
                    out_hits.append(hit)
                    new_passSeq.add(hit[0])
            for align in aln_list:
                if align[0] in new_passSeq:
                    out_aligns.append(align)
            hit_list, aln_list = out_hits, out_aligns

        if top:
            count = 0
            out_hits, out_aligns = [], []
            for hit in hit_list:
                count += 1
                if count > top:
                    break
                out_hits.append( hit )
            count = 0
            for align in aln_list:
                count += 1
                if count > top:
                    break
                out_aligns.append( align )
            hit_list, aln_list = out_hits, out_aligns
        hit_str = '\n'.join(['\t'.join(x) for x in hit_list])
        aln_str = '\n'.join(['\t'.join(x) for x in aln_list])

    else:
        hit_str, aln_str = None, None
        print( 'did not work' , flush = True)

    return hit_str, aln_str


def synthesizeHits( out_dict ):

    check = []
    for hit in out_dict:
        names = set( x.split('\t')[0] for x in out_dict[ hit ][0].split('\n') )
        check.append( names )
    
    passing = set()
    for names in check:
        for name in list(names):
            if all( name in x for x in check ):
                passing.add( name )
    
    new_dict = {}
    for hit in out_dict:
        hit_str, aln_str = '', ''
        for seq in out_dict[ hit ][ 0 ].split('\n'):
            if seq.split('\t')[0] in passing:
                hit_str += seq + '\n'
        for seq in out_dict[ hit ][ 1 ].split( '\n' ):
            if seq.split('\t')[0] in passing:
                aln_str += seq + '\n'
        new_dict[ hit ] = ( hit_str, aln_str )

    return new_dict

    
def main( 
    data, accession, best, threshold, evalue,
    bitscore, query = True, header = True 
    ):

    check, out = None, {}
    if accession:
#        if accession == True:
        check = grab_names(data)
        if not check:
            check = grab_names(data, True)
            accession = False
            query = True
    else:
#        if query == True:
        check = grab_names(data, True)
        if not check:
            accession = True
            query = False
            check = grab_names(data)

    if isinstance(check, set):
        for x in check:
            hits, align = '', ''
            if header: 
                hits = '#seq\tseq_e\tseq_score\tseq_bias\tdom_e\tdom_score\tdom_bias\texp\tN\n'
                align = '#seq\tdomNum\tincl_thresh\tscore\tbias\tc-Evalue\ti-Evalue\thmmFrom\t' + \
                                    'hmmTo\taliFrom\taliTo\tenvFrom\tenvTo\tacc\n'
            if accession:
                t_hits, t_aligns = grab_hits( data, threshold, evalue = evalue, \
                    top = best, accession = x )
            else:
                t_hits, t_aligns = grab_hits( data, threshold, evalue = evalue, \
                    top = best, query = x )
            if t_hits:
                hits += t_hits
                align += t_aligns
                out[x] = ( hits, align )
    else:
        htis, align = '', ''
        if header:
            hits = '#seq\tseq_e\tseq_score\tseq_bias\tdom_e\tdom_score\tdom_bias\texp\tN\n'
            align = '#seq\tdomNum\tincl_thresh\tscore\tbias\tc-Evalue\ti-Evalue\thmmFrom\t' + \
                            'hmmTo\taliFrom\taliTo\tenvFrom\tenvTo\tacc\n'
        t_hits, t_aligns = grab_hits(data, threshold, top = best, \
            query = query, accession = accession)
        hits += t_hits
        align += t_aligns
        if accession:
            out[ accession ] = ( hits, align )
        else:
            out[ query ] = ( hits, align )

    return out
            

def cli():

    parser = argparse.ArgumentParser( description = 'Extracts inputted accession hits and alignments from `hmmsearch` ' + \
        'output. Optionally provide a size-based query threshold to only extract alignments of a certain size.' )
    parser.add_argument( '-i', '--input', required = True, help = '`hmmsearch` output' )
    parser.add_argument( '-a', '--accession', \
        help = 'Accession to extract or new line delimited list of accessions; `-a`|`-q` required' )
    parser.add_argument( '-q', '--query', \
        help = 'Query to extract or new line delimited list of queries; `-a`|`-q` required' )
    parser.add_argument( '-b', '--best', type = int, help = '# top hits to take' )
    parser.add_argument( '-t', '--threshold', type = float, help = 'Percent query coverage to consider a hit (float).' )
    parser.add_argument( '-e', '--evalue', help = 'Negative magnitude of E-value, e.g. ' \
        + '10^(-x) where `x` is the input.', type = float, default = 0 )
    parser.add_argument( '--allq', action = 'store_true', help = 'Extract all queries' )
    parser.add_argument( '--alla', action = 'store_true', help = 'Extract all accessions' )
    parser.add_argument( '-o', '--output', help = 'Output file name/path (extensions automatically applied' )
    args = parser.parse_args()

    # initialize output file structure
    output = format_path(args.output)
    if not args.output:
        output = mkOutput(os.getcwd() + '/', 'extractHmmsearch')
    elif not os.path.isdir(output):
        os.mkdir(args.output)

    args_dict = { 
        'Input': args.input, 'Output': output, 'Accession': args.accession,
        'Query': args.query,  'Best hits': args.best, 'Threshold': args.threshold,
        'E-value': args.evalue
    }
    start_time = intro( 'hmmsearch Extraction', args_dict )


    if args.allq:
        args.query = True
    elif args.alla:
        args.accession = True
    else:
        if args.query:
            if os.path.isfile( args.query ):
                ques = file2list( args.query )
            else:
                ques = [ args.query ]
            args.query = True
        elif args.accession:
            if os.path.isfile( args.accession ):
                ques = file2list( args.accession )
            else:
                ques = [ args.accession ]
            args.accession = True
        else:
            print('\nNeed `-q` or `-a` specified\n', flush = True)
            sys.exit(8)
        

    if not os.path.isfile( args.input ):
        print('\n\tNot a valid input file\n', flush = True)
        sys.eixt(2)
    with open( args.input, 'r' ) as raw:
        data = raw.read()


    if args.allq or args.alla:
        out_dict = main( 
            data, args.accession, args.best, args.threshold, 10**-(args.evalue), query = args.query 
            )
    else:
        for que in ques:
            if args.query:
                args.query = que
            else:
                args.accession = que   
            temp_out_dict = main( 
                data, args.accession, args.best, args.threshold, 
                10**-(args.evalue), header = False, query = args.query 
                )
            out_dict = { **out_dict, **temp_out_dict }
        out_dict = synthesizeHits( out_dict )

    for name in out_dict:
        with open(output + '/' + name + '.hits.tsv', 'w') as out:
            out.write( out_dict[name][0] )
        with open(output + '/' + name + '.aligns.tsv', 'w') as out:
            out.write( out_dict[name][1] )

    outro( start_time )


if __name__ == '__main__':
    cli()
