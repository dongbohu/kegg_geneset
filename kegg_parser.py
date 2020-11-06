#!/usr/bin/env python

from datetime import date
import glob
import os
import re
import unicodedata
import requests

from biothings.utils.dataload import dict_sweep, unlist
import mygene

ROOT_URL="http://rest.kegg.jp/get/"

# ================================================================
# This function is copied from `slugify()` function in
# `annotation-refinery/slugify.py`
# ================================================================
def slugify(value, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces to hyphens.
    Remove characters that aren't alphanumerics, underscores, or hyphens.
    Convert to lowercase. Also strip leading and trailing whitespace.
    """
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
        value = re.sub('[^\w\s-]', '', value, flags=re.U).strip().lower()

    value = re.sub('[^\w\s-]', '', value).strip().lower()
    value = re.sub('[-\s]+', '-', value)
    return value

# ================================================================
# This function is based on `get_kegg_sets_members()` function in
# `annotation-refinery/process_kegg.py`
# ================================================================
def get_ksets_in_type(kegg_set_filename):
    """
    Reads in file with all the KEGG sets and associated members

    Arguments:
      kegg_sets_filename: path of the file that contains all KEGG sets and
      their members/genes for a specific type of sets (e.g. pathway, module, etc.).

    Returns:
      A dictionary whose keys are KEGG set IDs, and values are sets of the members/genes.

    """

    ksets = dict()
    with open(kegg_set_filename, 'r') as fh:
        for line in fh:
            tokens = line.strip().split()
            kset_id = tokens[0].split(':')[1]    # kset_id is listed first with a prefix

            if tokens[0].split(':')[0] == 'md':
                # Modules are prefixed with species and underscore
                kset_id = kset_id.split('_').pop()

            gene_id = tokens[1].split(':')[1]  # gene listed second, has prefix
            if kset_id not in ksets:
                ksets[kset_id] = list()
            if gene_id not in ksets[kset_id]:
                ksets[kset_id].append(gene_id)

    return ksets


# ================================================================
# This function is based on `get_kegg_set_info()` function in
# `annotation-refinery/process_kegg.py`
# ================================================================
def get_kset_info(kset_id, organism):
    """
    Requst kegg info based on `kset_id`, then parse the response and
    create a dictionary.

    Arguments:
      kset_id: str
      organism: str, organism's scientific name

    Returns:
      a dictionary with 'title' and 'abstract', which correspond to 'NAME'
      and 'DESCRIPTION' lines respectively in kset_info_filename.

    """

    kset_info = dict()

    url = ROOT_URL + kset_id
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Invalid response from {url}")
        return kset_info

    kset_type = None
    kset_title = ""

    # Strip the trailing blank lines in content text, then split and
    # parse it line by line:
    texts = response.text.rstrip('\n').split("\n")
    for line in texts:
        if line.startswith('ENTRY'):
            tokens = line.split()
            kset_type = tokens[-1]
        if line.startswith('NAME'):
            kset_title = ' '.join(line.split()[1:])

    kset_info['name'] = kset_title

    id_str = 'KEGG-' + slugify(organism) + "-"
    if kset_type:
        id_str += kset_type + "-"
    id_str += kset_id
    #kset_info["_id"] = id_str
    kset_info["_id"] = kset_id

    return kset_info


def rm_gene_fields(genes_info):
    """Remove some redundant fields returned by mygene.info API."""
    for g in genes_info:
        if 'query' in g:
            del g['query']
        if '_version' in g:
            del g['_version']

    return genes_info


def load_data(data_folder):
    """A generator that produces KEGG genesets."""

    # Global geneset information
    organism = "Homo sapiens"
    creator = "kegg_parser"
    is_public = True

    mg = mygene.MyGeneInfo()
    for ktype_filename in glob.glob(os.path.join(data_folder, "*")):
        ksets_in_type = get_ksets_in_type(ktype_filename)

        for kset_id, genes in ksets_in_type.items():
            kset_info = get_kset_info(kset_id, organism)

            # Populate other fields in `kset_info`:
            #mg_fields=['symbol', 'name', 'taxid']
            #genes_info = mg.getgenes(genes, fields=mg_fields)
            #kset_info['genes'] = rm_gene_fields(genes_info)
            kset_info['genes'] = genes
            kset_info['creator'] = creator
            kset_info['date'] = date.today().isoformat()
            kset_info['is_public'] = is_public
            kset_info['organism'] = organism

            # Clean up the dict
            kset_info = dict_sweep(kset_info)
            kset_info = unlist(kset_info)

            yield {'_id': kset_id, 'info': kset_info}


# Test harness: ~33 minutes w/o calling mygene API, 90 minutes when calling mygene API.
# 2,426 genesets, each with a unique "_id".
#   - "disease": 1,915
#   - "module":    173
#   - "pathway":   358
#
if __name__ == "__main__":
    import json

    genesets = load_data("./test_data")
    for gs in genesets:
        print(json.dumps(gs, indent=2))
