import os
import sys
import yaml, json, sys, itertools

# list of the main dependecies used in the workflow
dependencies = {
    "antismash": r"workflow/envs/antismash.yaml",
    "prokka": r"workflow/envs/prokka.yaml",
    "mlst": r"workflow/envs/mlst.yaml",
    "eggnog-mapper": r"workflow/envs/eggnog.yaml",
    "roary": r"workflow/envs/roary.yaml",
    "refseq_masher": r"workflow/envs/refseq_masher.yaml",
    "seqfu": r"workflow/envs/seqfu.yaml",
}


def get_dependency_version(dep, dep_key):
    """
    return dependency version tags given a dictionary (dep) and its key (dep_key)
    """
    with open(dep[dep_key]) as file:
        result = []
        documents = yaml.full_load(file)
        for i in documents["dependencies"]:
            if i.startswith(dep_key):
                result = i.split("=")[-1]
    return str(result)


def write_dependecies_to_json(outfile, dep=dependencies):
    """
    write dependency version to a json file
    """
    with open(outfile, "w") as file:
        dv = {}
        for ky in dep.keys():
            vr = get_dependency_version(dep, ky)
            dv[ky] = vr
        json.dump(
            dv,
            file,
            indent=2,
        )
        file.close()
    return dv


if __name__ == "__main__":
    write_dependecies_to_json(sys.argv[1])
