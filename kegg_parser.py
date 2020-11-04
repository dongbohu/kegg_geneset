#!/usr/bin/env python

from datetime import date
import glob
import os
import requests

from biothings.utils.dataload import dict_sweep, unlist
import mygene

ROOT_URL="http://rest.kegg.jp/get/"


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
    kset_title = None

    # Strip the trailing blank lines in content text, then split and parse
    # it line by line:
    texts = response.text.rstrip('\n').split("\n")
    for line in texts:
        if line.startswith('ENTRY'):
            tokens = line.split()
            kset_type = tokens[-1]
        if line.startswith('NAME'):
            kset_title = ' '.join(line.split()[1:])

    if kset_title:
        kset_info['_id'] = (
            'KEGG-' + kset_type + '-' + kset_id + ': ' + kset_title
        )

    return kset_info


def load_data(data_folder):
    """A generator that produces KEGG genesets."""

    # Global geneset information
    organism = "Homo sapiens"
    creator = "kegg_parser"
    is_public = True

    for ktype_filename in glob.glob(os.path.join(data_folder, "*")):
        ksets_in_type = get_ksets_in_type(ktype_filename)

        for kset_id, genes in ksets_in_type.items():
            kset_info = get_kset_info(kset_id, organism)

            # Populate other fields in `kset_info`:
            kset_info['genes'] = [int(x) for x in genes]
            kset_info['creator'] = creator
            kset_info['date'] = date.today().isoformat()
            kset_info['is_public'] = is_public
            kset_info['organism'] = organism

            # TODO: get genes info from mygene API
            kset_info = dict_sweep(kset_info)
            kset_info = unlist(kset_info)

            yield kset_info


# Test harness
if __name__ == "__main__":
    import json

    genesets = load_data("./test_data")
    for gs in genesets:
        print(json.dumps(gs, indent=2))
