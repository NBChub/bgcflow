import json
import logging
import sys

import yaml

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

# list of the main dependecies used in the workflow
dependencies = {
    "antismash": r"workflow/envs/antismash.yaml",
    "bigslice": r"workflow/envs/bigslice.yaml",
    "cblaster": r"workflow/envs/cblaster.yaml",
    "prokka": r"workflow/envs/prokka.yaml",
    "eggnog-mapper": r"workflow/envs/eggnog.yaml",
    "roary": r"workflow/envs/roary.yaml",
    "seqfu": r"workflow/envs/seqfu.yaml",
    "checkm": r"workflow/envs/checkm.yaml",
    "gtdbtk": r"workflow/envs/gtdbtk.yaml",
    "gecco": r"workflow/envs/gecco.yaml",
}


def get_dependency_version(dep, dep_key, antismash_version="7"):
    """
    return dependency version tags given a dictionary (dep) and its key (dep_key)
    """
    if dep_key == "antismash":
        logging.info(f"AntiSMASH version is: {antismash_version}")
        if antismash_version == "6":
            dep[dep_key] = "workflow/envs/antismash_v6.yaml"
    logging.info(f"Getting software version for: {dep_key}")
    with open(dep[dep_key]) as file:
        result = []
        documents = yaml.full_load(file)
        for i in documents["dependencies"]:
            if type(i) == str:
                if i.startswith(dep_key):
                    result = i.split("=")[-1]
            elif type(i) == dict:
                assert list(i.keys()) == ["pip"], i.keys()
                for p in i["pip"]:
                    if dep_key in p:
                        if p.startswith("git+"):
                            result = p.split("@")[-1]
                            result = result.replace("-", ".")
                        else:
                            result = p.split("=")[-1]
    return str(result)


def write_dependecies_to_json(outfile, antismash_version, dep=dependencies):
    """
    write dependency version to a json file
    """
    with open(outfile, "w") as file:
        dv = {}
        for ky in dep.keys():
            vr = get_dependency_version(dep, ky, antismash_version=antismash_version)
            dv[ky] = vr
        json.dump(
            dv,
            file,
            indent=2,
        )
        file.close()
    return dv


if __name__ == "__main__":
    write_dependecies_to_json(sys.argv[1], sys.argv[2])
