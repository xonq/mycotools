#! /usr/bin/env python3

from Bio import Entrez
import pandas as pd 
import numpy as np
import time

# generates ome codes from biosample accession (default) and first 3 letters of genus and species
def gen_omes(df,column='#Organism/Name', tag=0):

	newdf = df
	newdf['internal_ome'] = np.nan
	tax_set = set()
	access = '-'

	for index in range(len(df)):
		taxonomy = df[column][index]
		taxonomy = taxonomy.split(' ')
		genus = taxonomy[0]
		species = taxonomy[1]
		name = genus[:3] + species[:3]
		if tag != 0:
			access = df[tag][index]
		number = 1

		while name + str(number) in tax_set:
			number += 1

		tax_set.add(name + str(number))

		if access != '-':
			ome = name + str(number) + '.' + access
		else:
			ome = name + str(number)

		newdf['internal_ome'][index] = ome
	
	return newdf

# collects paths to download proteomes and assemblies
def collect_ftps(ncbi_df,api_key=0,column='Assembly Accession',database="assembly"):

    print('\n\nAssembling ftp directories ...')
    ass_prots = {}
    count = 0

# for each row in the assembly, grab the accession number, form the search term for Entrez, use Entrez,
    for index in range(len(ncbi_df)):
        accession = ncbi_df[column][index]
        ome = ncbi_df['internal_ome'][index]
        print('\t' + ome)
        ass_prots[ome] = {}

        search_term = accession + '[' + column + ']'
        handle = Entrez.esearch(db=database, term=search_term)
        genome_id = Entrez.read(handle)['IdList']

# if more than one genome id is found, this needs to be reported because we're only using the first (shouldn't happen)
        if len(genome_id) > 1:
            print('\t' + accession, 'yields more than one genome ID:')
            for ID in genome_id:
                print('\t\t' + str(ID))
            print('\t\tUsing', str(genome_id[0]))

# obtain the path from a summary of the ftp directory and create the standard paths for proteomes and assemblies
        ID = genome_id[0]
        handle = Entrez.esummary(db=database, id=ID, report="full")
        record = Entrez.read(handle)
        ftp_path = record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        basename = os.path.basename(ftp_path)
        assembly = ftp_path + '/' + basename + '_genomic.fna.gz'
        proteome = ftp_path + '/' + basename + '_protein.faa.gz'

        ass_prots[ome]['assembly'] = assembly
        ass_prots[ome]['proteome'] = proteome

# if no API key is used, we can only generate 3 queries per second, otherwise we can use 10
        count += 1
        if api_key == 0:
            if count == 3:
                time.sleep( 1 )
                count = 0
        else:
            if count == 10:
                time.sleep( 1 )
                count = 0

    return ass_prots


def gather_taxonomy(ncbi_df,api_key=0):

    print('\n\nGathering taxonomy information ...')
    tax_dicts = []
    count = 0

    for i,ncbi_df_row in ncbi_df.iterrows():

        taxid = ncbi_df_row['TaxID']
        handle = Entrez.efetch(db="Taxonomy", id=taxid, remode="xml")
        records = Entrez.read(handle)
        ome = ncbi_df_row['internal_ome']
        print('\t' + ome)

        lineages = records[0]['LineageEx']
        tax_dicts.append({})
        tax_dicts[-1]['internal_ome'] = ome
        for tax in lineages:
             tax_dicts[-1][tax['Rank']] = tax['ScientificName']
    
        count += 1
        if count == 3 or count == 10:
            if api_key == 0:
                time.sleep( 1 )
                count = 0
            elif api_key == 1:
                if count == 10:
                    time.sleep( 1 )
                    count = 0

    return tax_dicts


def assimilate_taxonomy(ncbi_db,tax_dicts,forbid=['no rank', 'superkingdom', 'subkingdom', 'genus', 'species', 'species group', 'varietas', 'forma']):

    tax_df = pd.DataFrame(tax_dicts)
    for forbidden in forbid:
        if forbidden in tax_df.keys():
            del tax_df[forbidden]

    output = {}
    keys = tax_df.keys()
    for i,row in tax_df.iterrows():
        output[row['internal_ome']] = {}
        for key in keys:
            output[row['internal_ome']][key] = row[key]

    for index in range(len(ncbi_db)):
        ncbi_db['taxonomy'][index] = str(output[ncbi_db['internal_ome'][index]])

    return ncbi_db
