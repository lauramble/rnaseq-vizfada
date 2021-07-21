#!/usr/bin/python3.8

import requests
import flatten_json as fj
import pandas as pd
import sys
import subprocess

installed = set(sys.modules.keys())
required = {"flatten-json", "requests"}
missing = required - installed

if missing:
    python = sys.executable
    for pkg in missing:
        subprocess.check_call(
            [python, '-m', 'pip', 'install', pkg], stdout=subprocess.DEVNULL)


def json_to_df(json, flatten=True, root_keys_to_ignore=set()):
    ignore = root_keys_to_ignore
    if (flatten):
        for i in range(len(json)):
            ignore = ignore.union({key for key in json[i]["_source"].keys(
            ) if isinstance(json[i]["_source"][key], list)})
        print(ignore)
        temp = [dict({"id": json[i]["_id"]}, **fj.flatten_json(json[i]
                     ["_source"], root_keys_to_ignore=ignore)) for i in range(len(json))]
    else:
        temp = [dict({"id": json[i]["_id"]}, **json[i]["_source"])
                for i in range(len(json))]
    temp = pd.DataFrame(temp)
    temp.index = temp.id
    return temp


def get_from_faang_api(url, post_data, size=20000, verbose=False):
    newURL = url + '/?size=' + str(size)
    response = requests.post(newURL, json=post_data)
    result = response.json()
    if verbose:
        print(f"Number of hits: {result['hits']['total']}")
    if (result['hits']['total'] > size):
        print(
            f"WARNING: more hits than expected ({size}).\nA new request will be made to include every hit.")
        return(get_from_faang_api(url, post_data, size=result['hits']['total'], verbose=verbose))
    else:
        return result['hits']['hits']


def match_exps_to_files(experiments, files, verbose=False):
    exps = {d["_id"]: {} for d in experiments}
    noExp = []

    for f in files:
        x = f["_source"]["experiment"]["accession"]
        r = f["_source"]["run"]["accession"]
        if x in exps.keys():
            exps[x].setdefault(r, [])
            furl = f["_source"]["url"]
            try:
                furl = 'ftp://' + furl[furl.index('ftp.sra.ebi'):]
            except ValueError:
                if verbose: print(f"Unusual file url: {furl}")
            exps[x][r].append(furl)
        else:
            if verbose: print(f"WARNING: {x} not found in RNA-Seq experiments !")
            noExp.append(x)
    if verbose and noExp:
        print("WARNING: The following experiments were not included in the list of experiments from FAANG :")
        print("\n\t".join(noExp))

    exps = {k: exps[k] for k in exps.keys() if len(exps[k]) != 0}
    return(exps)


def split_exp_inputs(expList, n, name="input", all=False, verbose=False):
    tot = len(expList)
    imax = tot // n
    ilist = [str(i).zfill(len(str(imax))) for i in range(0, imax+1)]

    if all:
        with open(f"{name}.txt", 'w') as f:
            f.write('\n'.join(expList))
        if verbose: print(f"Saved full list to {name}.txt")

    for a, i in zip(range(0, tot+n, n), ilist):
        with open(f"{name}_{i}.txt", 'w') as f:
            f.write('\n'.join(expList[a:min(a+n, tot)]))
        if verbose: print(f"{name}_{i}.txt")


# MAIN
if __name__ == "__main__":
    species = sys.argv[1]
    nExp = sys.argv[2]

    post_experiment = {
        "query": {
            "bool": {
                "must": [
                    {"match": {"standardMet": "FAANG"}},
                    {"exists": {"field": "RNA-seq"}}
                ]
            }
        }
    }

    post_files = {
        "query": {
            "bool": {
                "must": [
                    {"match": {"experiment.standardMet": "FAANG"}},
                    {"match": {"species.text": species}}
                ]
            }
        }
    }

    post_specimens = {
        "query": {
            "bool": {
                "must": [
                    {"match": {"organism.organism.text": species}},
                    {"match": {"standardMet": "FAANG"}}
                ]
            }
        }
    }

    # GET DATA
    experiments = get_from_faang_api(
        'http://data.faang.org/api/experiment/_search', post_data=post_experiment, verbose=True)
    files = get_from_faang_api(
        'http://data.faang.org/api/file/_search', post_data=post_files, verbose=True)
    specimens = get_from_faang_api(
        'http://data.faang.org/api/specimen/_search', post_data=post_specimens, verbose=True)

    # MATCH INPUT AND DATA

    exps = match_exps_to_files(experiments, files)

    split_exp_inputs(list(exps.keys()).sort(), nExp)

    # with open(f'rnaseq_{species.replace(' ', '_')}.json', 'w') as f:
    #    f.write(json.dumps(exps, indent=2, sort_keys=True))

    # GET METADATA

    filesDf = json_to_df(files)
    experimentsDf = json_to_df(experiments)
    experimentsDf = experimentsDf[experimentsDf.index.isin(exps.keys())]
    # inputDf = json_to_df(input_dna)
    specDf = json_to_df(specimens)

    # allChipDf = inputDf.append(experimentsDf)
    meta = pd.merge(filesDf,
                    experimentsDf,
                    how='right',
                    left_on="experiment_accession",
                    right_on=experimentsDf.index,
                    suffixes=["_file", "_exp"],
                    sort=True
                    )
    meta = pd.merge(meta, specDf,
                    how="left",
                    left_on="specimen",
                    right_on=specDf.index,
                    suffixes=["", "_spec"])

    meta.dropna(how='all', axis=1, inplace=True)

    with open("metadata.tsv", "w") as f:
        f.write(meta.to_csv(index=False, sep="\t"))


"""
d = [len(exps[id]) for id in exps.keys()]

######### MAKE DESIGN FILES

paired = []
single = []
done = []

for f in files:
    fname = f["_id"]
    #print(fname)
    furl = "ftp://" + f["_source"]["url"]
    fexp = f["_source"]["experiment"]["accession"]
    if fname[-2:] != "_2":
        done.append(fexp)
        design = {"group" : fexp, "replicate": done.count(fexp)}
        design["fastq_1"] = furl
        design["fastq_2"] = ""
        design["antibody"] = f["_source"]["experiment"]["target"]
        design["control"] = all_chip[fexp]
        if fname[-2] == "_":
            design["fastq_2"] = furl[:-10] + "2" + furl[-9:]
            paired.append(design)
        else:
            single.append(design)

# PAIRED

if paired:
    df = pd.DataFrame(paired)
    df.replace("input DNA", "", inplace=True)
    df=df.sort_values(["control", "antibody",  "group", "replicate"])

    print(df)

    with open(f"paired.csv", 'w') as f:
        f.write(df.to_csv(index=False))

    nExpsPerControl = df['control'].value_counts()
    controls = nExpsPerControl.index
    nGroups = df['group'].value_counts()
    allGroups = nGroups.index

    for inputDNA in controls:
        if inputDNA in allGroups:
            csv = df[df.control == inputDNA]
            csv = csv.append(df[df.group == inputDNA])
            with open(f"paired_{inputDNA}.csv", 'w') as f:
                f.write(csv.to_csv(index=False))

# SINGLE

if single:
    df = pd.DataFrame(single)
    df.replace("input DNA", "", inplace=True)
    df=df.sort_values(["control", "antibody",  "group", "replicate"])

    print(df)

    with open(f"single.csv", 'w') as f:
        f.write(df.to_csv(index=False))

    nExpsPerControl = df['control'].value_counts()
    controls = nExpsPerControl.index
    nGroups = df['group'].value_counts()
    allGroups = nGroups.index

    for inputDNA in controls:
        if inputDNA in allGroups:
            csv = df[df.control == inputDNA]
            csv = csv.append(df[df.group == inputDNA])
            with open(f"single_{inputDNA}.csv", 'w') as f:
                f.write(csv.to_csv(index=False))



"""
