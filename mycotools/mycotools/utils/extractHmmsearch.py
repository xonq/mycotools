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


def grab_hits(data, threshold = None, evalue = None, 
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
            hit_str += ' '.join(line[8:]) + '\t' + line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + \
                line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\n'

        passingSeq, align_str, count = set(), '', 0
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
                    

                align_str += str(seq) + '\t' + outLine[0] + '\t' + outLine[1] + '\t' + outLine[2] + \
                    '\t' + outLine[3] + '\t' + outLine[4] + '\t' + outLine[5] + '\t' + outLine[6] + '\t' + \
                    outLine[7] + '\t' + outLine[8] + '\t' + outLine[9] + '\t' + outLine[10] + '\t' + \
                    outLine[11] + '\t' + outLine[12] + '\n'

        if threshold:
            hit_list = [ hit for hit in hit_str.split( '\n' ) if hit != '' ]
            out_hits = []
            for hit in hit_list:
                start = re.search(r'^(.*?)\t', hit)[1]
                if start.rstrip() in passingSeq:
                    out_hits.append( hit )
            hit_str = '\n'.join( out_hits )

        if evalue:
            new_passSeq = set()
            hit_list = [x for x in hit_str.split('\n') if x]
            align_list = [x for x in align_str.split('\n') if x]
            out_hits, out_aligns = [], []
            for hit in hit_list:
                if int(10**200 * float( hit.split('\t')[1] )) < int(10**200 * float( '1e-' + str(evalue) )):
                    out_hits.append( hit )
                    new_passSeq.add( hit.split('\t')[0] )
            for align in align_list:
                if align.split('\t')[0] in new_passSeq:
                    out_aligns.append( align )
            hit_str = '\n'.join( out_hits )
            align_Str = '\n'.join( out_aligns )

        if top:
            count = 0
            hit_list, align_list = hit_str.split( '\n' ), align_str.split( '\n' )
            out_hits, out_aligns = [], []
            for hit in hit_list:
                count += 1
                if count > top:
                    break
                out_hits.append( hit )
            count = 0
            for align in align_list:
                count += 1
                if count > top:
                    break
                out_aligns.append( align )
            hit_str = '\n'.join( out_hits )
            align_str = '\n'.join( out_aligns )

    else:
        hit_str, align_str = None, None
        print( 'did not work' , flush = True)

    return hit_str, align_str


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
        hit_str, align_str = '', ''
        for seq in out_dict[ hit ][ 0 ].split('\n'):
            if seq.split('\t')[0] in passing:
                hit_str += seq + '\n'
        for seq in out_dict[ hit ][ 1 ].split( '\n' ):
            if seq.split('\t')[0] in passing:
                align_str += seq + '\n'
        new_dict[ hit ] = ( hit_str, align_str )

    return new_dict

    
def main( 
    data, accession, best, threshold, evalue,
    query = True, header = True 
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
            

if __name__ == '__main__':

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
        + '10^(-x) where `x` is the input.' )
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
            data, args.accession, args.best, args.threshold, args.evalue, query = args.query 
            )
    else:
        for que in ques:
            if args.query:
                args.query = que
            else:
                args.accession = que   
            temp_out_dict = main( 
                data, args.accession, args.best, args.threshold, 
                args.evalue, header = False, query = args.query 
                )
            out_dict = { **out_dict, **temp_out_dict }
        out_dict = synthesizeHits( out_dict )

    for name in out_dict:
        with open(output + '/' + name + '.hits.tsv', 'w') as out:
            out.write( out_dict[name][0] )
        with open(output + '/' + name + '.aligns.tsv', 'w') as out:
            out.write( out_dict[name][1] )

    outro( start_time )
