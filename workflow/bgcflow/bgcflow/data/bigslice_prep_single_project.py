from pathlib import Path
import pandas as pd
import sys, logging

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

# create intermediate folder


def prepare_bigslice(input_folder, taxonomy_file, project_name, output_folder):
    """
    Preparing bigslice input folder for a single project
    """
    logging.info(f"Preparing bigslice input folder for {project_name}")
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    assert input_folder.is_dir()
    output_folder.mkdir(parents=True, exist_ok=True)

    # create symlink to antismash
    logging.info("Creating symlink to AntiSMASH BGCs...")
    antismash_dir = output_folder / project_name
    try:
        antismash_dir.symlink_to(input_folder.resolve(), target_is_directory=True)
    except FileExistsError:
        pass

    # create tax info
    logging.info("Getting taxonomic information...")
    taxonomy_dir = output_folder / "taxonomy"
    taxonomy_dir.mkdir(parents=True, exist_ok=True)
    df_tax = pd.read_csv(taxonomy_file, sep="\t")
    df_tax = df_tax.loc[
        :,
        [
            "#Genome Folder",
            "Domain",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
            "Organism",
        ],
    ]
    df_tax.to_csv(taxonomy_dir / f"{project_name}.tsv", sep="\t", index=False)

    # build metadata
    logging.info("Preparing metadata...")
    metadata = {
        "# Dataset name": project_name,
        "Path to folder": f"{output_folder.name}/",
        "Path to taxonomy": f"taxonomy/{project_name}.tsv",
        "Description": project_name,
    }
    df_metadata = pd.DataFrame.from_dict(metadata, orient="index").T
    df_metadata.to_csv(output_folder / "datasets.tsv", sep="\t", index=False)
    logging.info("Job done!")
    return


if __name__ == "__main__":
    prepare_bigslice(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
