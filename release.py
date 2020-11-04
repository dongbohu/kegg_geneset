def get_release(self):
    import requests

    release = ""
    kegg_release_resp = requests.get("http://rest.kegg.jp/info/kegg")
    for line in kegg_release_resp.text.split("\n"):
        tokens = line.strip().split()
        if len(tokens) > 1 and tokens[0] == "kegg" and tokens[1] == "Release":
            release = " ".join(tokens[1:])
            break

    return release
