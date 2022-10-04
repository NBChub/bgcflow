import sys
from pathlib import Path
import json, yaml
import peppy
import logging
log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def refine_bgcflow_project(p_bgcflow, p):
    """
    Refine a pep bgcflow project created from sample table.
    Exist for back compatibility with bgcflow=<0.3.3

    Arguments:
        p_bgcflow:
        p:

    Returns:
        p_bgcflow:
    """
    for k in p.keys():
        if k == 'samples':
            p_bgcflow.config['sample_table'] = p[k]
            p_bgcflow['_config_file'] = str(Path(p[k]).parent / "project_config.yaml")
        elif k in p_bgcflow.keys():
            p_bgcflow[k] = p[k]
        elif k == 'rules':
            with open(p['rules'], "r") as file:
                rules = yaml.safe_load(file)
                p_bgcflow.config[k] = rules['rules']
        else:
            p_bgcflow.config[k] = Path(p[k]).name

    return p_bgcflow

def get_bgcflow_metadata(bgcflow_path="."):
    # BGCFlow config
    BGCFlow_path = Path(bgcflow_path).resolve()

    logging.info("Getting config metadata...")
    config_path = BGCFlow_path / "config/config.yaml"
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    logging.info("Getting rules information...")
    # rules_dict
    rules_dict_path = BGCFlow_path / "workflow/rules/rules_main.json"
    with open(rules_dict_path, "r") as f:
        rules_dict = json.load(f)
    
    return BGCFlow_path, config, rules_dict

def get_all_metadata(config, rules_dict, BGCFlow_path, bgcflow_version):
    
    # Get Metadata
    logging.info("Getting metadata from projects...")
    project_metadata = {}

    for p in config['projects']:
        # process only pep projets
        if p['name'].endswith(".yaml"):
            project = peppy.Project(str(BGCFlow_path / p['name']), sample_table_index="genome_id")
            if 'description' in project.config.keys():
                    name = project['name']
                    description = project.config['description']
                    project_metadata[name] = {'description' : description}
        else:
            # non peppy input
            project = peppy.Project(str(BGCFlow_path / p['samples']), sample_table_index="genome_id")
            p = refine_bgcflow_project(project, p)
            project_metadata[p['name']] = "No description provided."

        # get what rules are being used
        rule_used = {}
        if 'rules' in project.config.keys():
            rules = project.config["rules"]
        else:
            rules = config["rules"]
        for r in rules.keys():
            if rules[r]:
                bgcflow_rules = rules_dict[r]
                rule_used[r] = bgcflow_rules
        project_metadata[name].update({'rule_used' : rule_used})
        
        # get sample size
        project_metadata[name]['sample_size'] = len(project.sample_table)
        
        # get citations
        citation_all = []
        for r in rule_used:
            citations = rules_dict[r]['references']
            citation_all.extend(citations)
        citation_all.sort()
        project_metadata[name].update({'references' : citation_all})
        project_metadata[name]['references'] = list(set(project_metadata[name]['references']))

        # get bgcflow_version
        project_metadata[name]['bgcflow_version'] = bgcflow_version
    return project_metadata

def get_project_metadata(project_name, outfile, bgcflow_path=".", bgcflow_version="unknown"):
    BGCFlow_path, config, rules_dict = get_bgcflow_metadata(bgcflow_path)
    all_metadata = get_all_metadata(config, rules_dict, BGCFlow_path, bgcflow_version)
    logging.info(f"Extracting project {project_name} metadata to {outfile}")
    with open(outfile, "w") as f:
        json.dump({project_name : all_metadata[project_name]}, f, indent=2)
    return

if __name__ == "__main__":
    get_project_metadata(sys.argv[1], sys.argv[2], bgcflow_version=sys.argv[3])