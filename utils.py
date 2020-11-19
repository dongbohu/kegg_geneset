#!/usr/bin/env python3

"""Utility functions"""

import requests

KEGG_INFO_URL = "http://rest.kegg.jp/info/kegg"


def get_release(self):
    release = ""
    resp = requests.get(KEGG_INFO_URL)
    if resp.status_code != 200:
        raise Exception(f"Failed to request {url}")

    text_lines = resp.text.strip('\n').split('\n')
    for line in text_lines:
        tokens = line.strip().split()
        if len(tokens) > 1 and tokens[0] == "kegg" and tokens[1] == "Release":
            release = " ".join(tokens[1:])
            break

    return release


# Test harness
if __name__ == "__main__":
    print(get_release(None))
