#!/usr/bin/env python3

import os
import glob
import argparse
import xml.etree.ElementTree as ET
from subprocess import call as unix
from joblib import Parallel, delayed


# Author: Peter Mulhair
# Date: 04/03/2022
# Usage: python3 get_genes.py --genomes <ncbi xml file> --genes <ensembl csv file>

'''
Script to download available annotation 
data for insect species from DToL. A csv 
file downloaded from ensembl rapid release 
site is required to get species names with 
related data.
'''

parse = argparse.ArgumentParser()

parse.add_argument("--genomes",type=str, help="ncbi xml file with genome info",required=True)
parse.add_argument("--genes",type=str, help="ensembl csv file with gene annotation info",required=True)

args = parse.parse_args()

#Check for species with data already downloaded
sp_complete = []
for sp in glob.glob('proteins/*fa'):
    sp_name = sp.split('/')[-1].split('-')[0]
    sp_complete.append(sp_name)

#Parse ensembl csv file to get species names with annotation data
GCA_list = []
with open(args.genes) as f:
    next(f)
    for line in f:
        GCA = line.split('","')[4]
        GCA_list.append(GCA)
        

#Parse ncbi xml file to create dictionary of species name to GCA assembly ID
genome_dict = {}
with open(args.input) as f:
    xml = f.read()
    
ncbi_file = ET.fromstring("<root>" + '\n' + xml + "</root>")
for elem in ncbi_file:
    
    for subelem in elem.findall('AssemblyName'):
        spID = subelem.text
        
    for subelem in elem.findall('AssemblyAccession'):
        GCA = subelem.text
        GCAid = GCA.split('.')[0]
        if 'GCF' in GCA:
            GCA = 'GCA_' + GCA.split('_')[1]
            GCAid = 'GCA_' + GCAid.split('_')[1]
        
    for subelem in elem.findall('Organism'):
        sp_name = subelem.text
        sp_name = sp_name.split(' (')[0]
        
    if 'alternate' not in spID:
        if spID[:1] == 'i':
            if GCA in GCA_list:
                sp_name = sp_name.replace(' ', '_')
                if sp_name == 'Nymphalis_io':
                    sp_name = 'Inachis_io'
                print(spID[:2])
                if sp_name not in sp_complete:
                    genome_dict[sp_name] = GCA

if len(genome_dict) > 0:
    print('Downloading data for', len(genome_dict), 'species')
else:
    print('No new data to download')
        
num_list = ['04', '05', '06', '07', '08', '09', '10', '11', '12']

def genome_download(species, genome):
    os.makedirs('proteins', exist_ok=True)
    os.makedirs('cds', exist_ok=True)
    os.makedirs('gff', exist_ok=True)
    os.makedirs('gtf', exist_ok=True)
    
    print(species)

    os.chdir('proteins/')
    for num in num_list:
        try:
            unix('sudo wget -q ftp://ftp.ensembl.org/pub/rapid-release/species/' + species + '/' + genome + '/geneset/2021_' + num + '/*pep*', shell=True)
        except:
            continue
    unix('sudo gzip -d *gz', shell=True)
    os.chdir('../cds/')
    for num in num_list:
        try:
            unix('sudo wget -q ftp://ftp.ensembl.org/pub/rapid-release/species/' + species + '/' + genome + '/geneset/2021_' + num + '/*cds*', shell=True)
        except:
            continue
    unix('sudo gzip -d *gz', shell=True)
    os.chdir('../gff/')
    for num in num_list:
        try:
            unix('sudo wget -q ftp://ftp.ensembl.org/pub/rapid-release/species/' + species + '/' + genome + '/geneset/2021_' + num + '/*gff*', shell=True)
        except:
            continue
    unix('sudo gzip -d *gz', shell=True)

    os.chdir('gtf/')
    for num in num_list:
        try:
            unix('sudo wget -q ftp://ftp.ensembl.org/pub/rapid-release/species/' + species + '/' + genome + '/geneset/2021_' + num + '/*gtf*', shell=True)
        except:
            continue
    unix('sudo gzip -d *gz', shell=True)
    os.chdir('../') 

if len(genome_dict) > 0:
    Parallel(n_jobs=1)(delayed(genome_download)(k, v) for k, v in genome_dict.items())  
