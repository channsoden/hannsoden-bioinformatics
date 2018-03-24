#!/usr/bin/env python
""" This Python script contains dictionaries of fasta genomes for use in the WholeGenomePhylogeny.py
pipeline and other bioinformatics tools by Christopher Hann-Soden.

A Python dictionary is a set of key:value pairs. A valid dictionary for these tools should be composed
of two strings (surrounded by " or ' characters). The keys of the dictionary should be unique, short or
abbreviated names for the genomes. The values of the dictionary should be the full path to the fasta
file of the genome.

The subgroup functions and merge_dict functions can be used to easily create subsets or supersets of
datasets.
"""

import re
def subgroup_shortname(regex, d):
    """Returns a dictionary composed of elements whose keys (short names) contain the regex."""
    grp = {k:v for k,v in list(d.items()) if re.search(regex, k)}
    return grp
def subgroup_longname(regex, d):
    """Returns a dictionary composed of elements whose values (long names) contain the regex."""
    grp = {k:v for k,v in list(d.items()) if re.search(regex, v)}
    return grp

def merge_dicts(d1, d2):
    """Returns a dictionary composed of elements from both d1 and d2."""
    new = {k:v for k, v in list(d1.items())}
    for k, v in list(d2.items()):
        new[k] = v
    return new

molds = {'Nc': '/global/home/users/username/genomes/Neurospora_crassa_OR74A_FungiDB-3.1.fasta',
         'Nd': '/global/home/users/username/genomes/Neurospora_discreta_FGSC8579_FungiDB-3.1.fasta',
         'Np': '/global/home/users/username/genomes/Neurospora_pannonica_FGSC7221.fasta',
         'Ns': '/global/home/users/username/genomes/Neurospora_sublineolata_FGSC5508.fasta',
         'Sm': '/global/home/users/username/lustre_scratch/genomes/Sordaria_macrospora_k-hell_FungiDB-3.1.fasta',
         'Mg': '/global/home/users/username/genomes/magnaporthe_grisea_70-15.fasta',
         'Pa': '/global/home/users/username/genomes/Podospora_anserina_NCBI100915.fasta',
         'Ta': '/global/home/users/username/genomes/Trichoderma_atroviride_IMI_206040.17-JUL-2007.nt.fasta',
         'Cg': '/global/home/users/username/genomes/Chaetomium_globosum_fullgenome_NCBI100915.fasta',
         'Fo': '/global/home/users/username/genomes/Fusarium_oxysporum_f_sp_lycopersici_4286.02-JUL-2007.nt.fasta'}

yeasts = {'Sc': '/global/home/users/username/genomes/Saccharomyces_cerevisiae.final.scaffolds.fasta',
          'Sb': '/global/home/users/username/genomes/Saccharomyces_bayanus_2HBK.fasta',
          'Sp': '/global/home/users/username/genomes/Saccharomyces_paradoxus_4MFT.fasta',
          'Ca': '/global/home/users/username/genomes/Candida_albicans_NCBI_072014.fasta',
          'Scp': '/global/home/users/username/genomes/Schizosaccharomyces_pombe_v12.fasta',
          'Mv': '/global/home/users/username/genomes/Microbotryum_violaceum.final.scaffolds.fasta'}

neurospora = subgroup_shortname('^N', molds)
saccharomyces = subgroup_longname('Saccharomyces', yeasts)
fungi = merge_dicts(molds, yeasts)
