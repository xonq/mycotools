#! /usr/bin/env python3

from mycotools.lib.kontools import formatPath, sysStart, eprint
from mycotools.lib.fastatools import gff2dict, dict2gff, gff3Comps
import re, sys, os


def hashGff3IDs( gff_list, comps, tag = 'Alias' ):

    id_hash = {}
    for index in range(len(gff_list)):
        entry = gff_list[index]
        hit_id = re.search( comps['id'], entry['attributes'] )[1]
        par_id = re.search( comps['par'], entry['attributes'] )
        alias = re.search( comps[tag], entry['attributes'] )
        if par_id is not None:
            if par_id[1] not in id_hash:
                id_hash[par_id[1]] = [ [], None, [] ]
            if alias is not None:
                id_hash[par_id[1]][1] = alias[1]
            id_hash[par_id[1]][0].append( index )
            id_hash[par_id[1]][2].append( hit_id )

    return id_hash


def acquireFormat( gff_list ):

    prot_comp = re.compile( r'ID=(.*?)[;|$]' )
    gene = False
    for line in gff_list:
        if line['type'] == 'gene':
            gene = True
            break

    if gene:
        if prot_comp.search( line['attributes'] ):
            if 'Name=jgi.p' in line['attributes']:
                return 'jgi_gff3'
            else:
                return 'misc_gff3'
#    else:
 #       prot_comp = re.compile( r'proteinId (.*?);' )
  #      for line in gff_list:
   #         if line['type'] == 'CDS':
    #            break
     #   if prot_comp.search( line['attributes'] ):
      #      return 'jgi_gff'

    return None


#def indices2gff3( gff_list, indices ): 
 #   return [ gff_list[x] for x in indices ]


def checkGff3Alias( alias_dict, cur_list ):

    out_list = []
    for alias in alias_dict:
        temp_list = [cur_list[x] for x in alias_dict[alias]]
        for line in temp_list:
            if not line['attributes'].endswith(';Alias=' + alias):
                line['attributes'] += ';Alias=' + alias
        out_list.extend( temp_list )

    return out_list


def compileGenes( cur_list, par2id, comps = gff3Comps() ):

    genes = [
        (i, re.search(comps['id'], cur_list[i]['attributes'])[1]) \
        for i in range(len(cur_list)) \
        if cur_list[i]['type'] == 'gene'
        ]
    gene_dict = {}
    for info in genes:
        gene = info[1]
        gene_dict[gene] = {
            'children': [info[0]],
            'alias': None,
            'tocheck': { gene }
            }
        while gene_dict[gene]['tocheck']:
            hit = list(gene_dict[gene]['tocheck'])[0]
            if hit in par2id:
                gene_dict[gene]['children'].extend(
                    par2id[hit][0]
                    )
                gene_dict[gene]['tocheck'] = gene_dict[gene]['tocheck'] \
                    | set( par2id[hit][2] )
                if par2id[hit][1]:
                    gene_dict[gene]['alias'] = par2id[hit][1]
            gene_dict[gene]['tocheck'].remove(
                hit
                )

    alias_dict = { 
        gene_dict[x]['alias']: sorted(gene_dict[x]['children']) for x in gene_dict 
        }
    if None in alias_dict:
        del alias_dict[None]
    etc = [ gene_dict[x]['children'] for x in gene_dict if not gene_dict[x]['alias']]
    etc_list = []
    for hit in etc:
        for i in hit:
            etc_list.append( cur_list[i] )

    return alias_dict, etc_list  


def curGff3pass1( gff_list, prot_comp ):

    cur_list = []
    for line in gff_list:
        search = prot_comp.search( line['attributes'] )
        if search:
            prot = search[1]
            if re.search( r'Alias=[^;]$', line['attributes'] ) is None:
                if not line['attributes'].endswith( ';' ):
                    line['attributes'] += ';'
                line['attributes'] += 'Alias=' + ome + '_' + prot
        cur_list.append( line )

    return cur_list
                

def curGff3( gff_list, ome ):

    cur_list, prot_comp = [], None
    for line in gff_list:
    # this order is sketchy, really don't know if it's conserved
        if ';Alias=' in line['attributes']:
            cur_list = gff_list
        if ';proteinId=' in line['attributes'] \
            or line['attributes'].startswith('proteinId='):
            prot_comp = re.compile( r'(?:^|;)proteinId=([^;]+)' )
            break
#        elif re.search( r'[^|;]Name=', line['attributes'] ):
 #           prot_comp = re.compile( r'[^|;]Name=(.*?)[$|;]' )
  #          break
        elif ';protein_id=' in line['attributes'] \
            or line['attributes'].startswith('protein_id='):
            prot_comp = re.compile( r'(?:^|;)protein_id=([^;]+)' )
            break
    if not prot_comp:
        eprint('\tERROR: ' + ome + ' no known protein IDs')
        return gff_list

    if not cur_list:
        cur_list = curGff3pass1( gff_list, prot_comp )

    par2id = hashGff3IDs( cur_list, gff3Comps(), 'Alias' )
    alias_dict, etc_list = compileGenes( cur_list, par2id )
    final_list = checkGff3Alias( alias_dict, cur_list )
    final_list.extend( etc_list )
    

    return final_list

'''
def checkGff2Alias( temp_out, acc ):

    for line in temp_out:
        acc_search = re.search( r'alias "([^"]+)"$', line['attributes'] ) 
        if acc_search is None:
            if not line['attributes'].endswith( ';' ):
                line['attributes'] += ';'
            line['attributes'] += ' alias "' + acc + '"'

    return temp_out


def acc2gff2( out_indices, gff_list ):
    return [ gff_list[x] for x in out_indices ]


def hashGff2IDs( gff_list, comps, tag = 'alias' ):

    id_hash = {}
    for index in range(len(gff_list)):
        entry = gff_list[index]['attributes']
        hit_id = re.search( comps['id'], entry )[1]
        alias = re.search( comps[tag], entry )
        if hit_id not in id_hash:
            id_hash[hit_id] = [ None, [] ]
        if alias is not None:
            id_hash[hit_id][0] = alias[1]
        id_hash[hit_id][1].append( index )

    alias_dict = { 
        id_hash[x][0]: sorted(id_hash[x][1]) for x in id_hash
        }
    if None in alias_dict:
        del alias_dict[None]
    etc = [ id_hash[x][1] for x in id_hash if not id_hash[x][0]]
    etc_list = []
    for hit in etc:
        for i in hit:
            etc_list.append( gff_list[i] )

    return alias_dict, etc_list


def accFromGff2Hash( id_hash, acc ):
    return [id_hash[x][1] for x in id_hash if id_hash[x][0] == acc]


def curGff2( gff_list, ome ):

    cur_list = []
    prot_comp = re.compile( r'proteinId (.*?)[;|$]' )
    for line in gff_list:
        if line['type'] == 'CDS':
            prot = prot_comp.search( line['attributes'] )[1]
            alias = re.search( r' alias ([^;]*)$', line['attributes'] )
            if alias is None:
                if not line['attributes'].endswith( ';' ):
                    line['attributes'] += ';'
                line['attributes'] += ' alias "' + ome + '_' + prot + '"'
            elif not line['attributes'].endswith('"'):
                line['attributes'] = line['attributes'].replace(
                    ' alias ' + ome + '_' + prot, 
                    ' alias "' + ome + '_' + prot + '"'
                    )
        cur_list.append( line )

    alias_dict, etc_list = hashGff2IDs( cur_list, gff2Comps(), 'alias' )
    final_list = []
    for alias in alias_dict:
        prep_out = [ cur_list[x] for x in alias_dict[alias] ]
        final_list.extend(checkGff2Alias( prep_out, alias ))
    final_list.extend( etc_list )

    return final_list
'''

def main( gff_path, ome):

    gff = gff2dict( gff_path )
    typ = acquireFormat( gff )

    if not typ:
        eprint('\tERROR: type unknown ' + gff_path)
        return None

  #  if 'gff3' in typ:
    new_gff = curGff3( gff, ome )
#    else:
 #       new_gff = curGff2( gff, ome )

    return new_gff


if __name__ == '__main__':
    usage = 'Imports gene coordinates file gff3, ome, and curates headers'
    sysStart( sys.argv, usage, 3, files = [sys.argv[1]] )
    cur_gff = main( formatPath(sys.argv[1]), sys.argv[2] )
    print( dict2gff( cur_gff ) )
    sys.exit( 0 )
