import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    format="%(levelname)-8s %(asctime)s   %(message)s",
    datefmt="%d/%m %H:%M:%S",
    level=logging.DEBUG,
)


def annotate_mmseqs2_cog(mmseqs2_cog: str, gbk_path: str, outfile: str) -> None:
    """
    Annotate an MMseqs2 COG file with information from a GenBank file.

    Args:
        mmseqs2_cog: Path to the MMseqs2 COG file.
        gbk_path: Path to the GenBank file.
        outfile: Path to the output file.

    Returns:
        None
    """
    gbk_path = Path(gbk_path)
    mmseqs2_cog = Path(mmseqs2_cog)

    logging.info("Reading MMseqs2 COG file...")
    df_mmseqs2 = pd.read_csv(mmseqs2_cog)

    logging.info("Parsing GenBank file...")
    with open(gbk_path, "r") as handle:
        output = {}
        ctr = 0
        for record in SeqIO.parse(handle, "genbank"):
            file_id = gbk_path.stem
            seq_id = record.id
            for f in record.features:
                if f.type == "CDS":
                    value = {k: ",".join(v) for k, v in f.qualifiers.items()}
                    if "locus_tag" not in value.keys():
                        logging.warning(
                            f"Could not find locus_tag in {file_id} {seq_id}. Available keys: {value.keys()}"
                        )
                        for locus_tag in ["protein_id", "gene"]:
                            if locus_tag in value:
                                logging.info(
                                    f"Using {locus_tag} as locus_tag in {file_id} {seq_id}."
                                )
                                value["locus_tag"] = value[locus_tag]
                                break
                    value["file_id"] = file_id
                    value["seq_id"] = seq_id
                    value["start"] = int(f.location.start)
                    value["stop"] = int(f.location.end)
                    value["strand"] = f.location.strand
                    value["type"] = f.type
                    output[ctr] = value
                    ctr = ctr + 1
    df_annotation = pd.DataFrame.from_dict(output).T

    logging.info("Merging dataframes...")
    logging.info(f"Length of MMseqs2 COG dataframe: {df_mmseqs2.shape}")
    logging.info(f"Length of annotation dataframe: {df_annotation.shape}")
    df_mmseqs2.set_index("locus_tag").to_csv("mmseqs2_cog2.csv")
    df_annotation.set_index("locus_tag").to_csv("annotation2.csv")
    df = df_mmseqs2.merge(
        df_annotation, right_on="locus_tag", left_on="locus_tag", how="outer"
    )
    logging.info(f"Length of merged dataframe: {df.shape}")
    df.set_index("locus_tag").to_csv("merged2s.csv")
    assert len(df_mmseqs2) == len(
        df
    ), "Error: Merged dataframe has different length than original MMseqs2 COG dataframe."

    logging.info("Writing output file...")
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outfile, index=False)

    logging.info("Job done!")


if __name__ == "__main__":
    annotate_mmseqs2_cog(sys.argv[1], sys.argv[2], sys.argv[3])
