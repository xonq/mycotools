#! /usr/bin/env python3

import argparse, sys, os, time, urllib, urllib.request, urllib.parse, json
from pyclustering.cluster import cluster_visualizer
from pyclustering.cluster.optics import optics
from pyclustering.utils import read_sample

# return list of path files from argument string, space delimitted
def filelist(files):
    files = files.split()
    for i in range(len(files)):
        files[i] = os.path.normpath(files[i])
    return files


# remove values outside of tolerances
# tries to import file via pyclustering module
    # if it fails due to presence of headers then import as pandas df, remove header, save, reimport
def data_truncate(data,minrt,maxrt,minint):

    while True:
        try:
            with open(data, 'r') as chrom_data:
                imported_data = []
                for line in chrom_data:
                    if line[0] != '#':
                        if line[-1] == '\n':
                            line = line.rstrip()
                        imported_data.append(line)
            break
        except FileNotFoundError:
            print('One of more file(s) not found. Exit status 3.')
            sys.exit(3)
    curated_list = []
    if maxrt != 0:
        for row in imported_data:
            row = row.split(sep='\t')
            row = [float(i) for i in row]
            if row[0] >= minrt:
                if row[0] <= maxrt:
                    if row[2] >= minint:
                        curated_list.append(row)
    else:
        for row in imported_data:
            row = row.split(sep='\t')
            row = [float(i) for i in row]
            if row[0] >= minrt and row[2] >= minint:
                curated_list.append(row)

    return curated_list


def pyclustering_optics(curated_list,rt_toler,mz_toler,neighbors):

    mz_norm_coeff = rt_toler/mz_toler

    cluster_input = []
    for row in curated_list:
        cluster_input.append([row[0],row[1]*mz_norm_coeff])


    optics_instance = optics(cluster_input, rt_toler, neighbors)
    optics_instance.process()
    cluster_index_list = optics_instance.get_clusters()

    clusters = []
    try:
        for cluster in cluster_index_list:
            clusters.append([])
            for index in cluster:
                clusters[-1].append(curated_list[index])

        for cluster in clusters:
            cluster.sort(key = lambda x: x[0])
    except TypeError:
        print('\n\nNo clusters found')

    return clusters


def cluster_apex(clusters):

    apex_list = []
    for cluster in clusters:
        apex = 0
        for index in range(len(cluster)):
            if cluster[index][2] > apex:
                apex = cluster[index][2]
                apex_index = index
        apex_list.append([cluster[apex_index],apex_index])

    return apex_list

# atom_dict = { minc, maxc, minh, maxh, minn, maxn, mino, maxo, mins, maxs, minp, maxp }, values must be strings
# param_dict = { 'mindus', 'maxdus', 'msrange' }
# returns a list of dictionaries containing mf, unsat, and ppm info for each
def chemcalc( mass, atom_dict, param_dict, ppmax ):

    chemcalcURL = 'http://www.chemcalc.org/chemcalc/em'

    for key in atom_dict:
        atom_dict[key] = str(atom_dict[key])

    param_dict['monoisotopicMass'] = mass
    mfRange = 'C{minc}-{maxc}H{minh}-{maxh}N{minn}-{maxn}O{mino}-{maxo}S{mins}-{maxs}P{minp}-{maxp}'.format(**atom_dict)
    param_dict['mfRange'] = mfRange

    while True:
        try:
            response = urllib.request.urlopen(chemcalcURL, urllib.parse.urlencode(param_dict).encode('utf-8'))
            break
        except ValueError:
            print('Waiting on server')
            time.sleep( 5 )

    jsondata = response.read()
    rawdata = json.loads(jsondata.decode('utf-8'))

    cur_list_dict = []
    for result in rawdata['results']:
        if abs( result['ppm'] ) <= ppmax and result['unsat'] >= param_dict['mindus'] and result['unsat'] <= param_dict['maxdus']:
            cur_list_dict.append({ 'molForm': result['mf'], 'unsat': result['unsat'], 'ppm': result['ppm'] })

    return cur_list_dict
