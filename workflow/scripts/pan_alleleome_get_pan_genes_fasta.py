import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_data(roary_path):
    """
    Load data from Roary output.

    Parameters:
    roary_path (Path): Path to Roary output

    Returns:
    tuple: A tuple containing two pandas DataFrames
    """
    df_gene_presence_binary = pd.read_csv(
        roary_path / "df_gene_presence_binary.csv", index_col="Gene", low_memory=False
    )
    df_gene_presence_locustag = pd.read_csv(
        roary_path / "df_gene_presence_locustag.csv", index_col="Gene", low_memory=False
    )
    df_gene_presence_locustag.index = [
        str(i).replace("/", "_") for i in list(df_gene_presence_locustag.index)
    ]
    return df_gene_presence_binary, df_gene_presence_locustag


def parse_genbank_files(df_gene_presence_locustag, gbk_folder):
    """
    Parse GenBank files.

    Parameters:
    df_gene_presence_locustag (DataFrame): DataFrame with gene presence information
    gbk_folder (Path): Path to folder containing GenBank files

    Returns:
    DataFrame: DataFrame with parsed GenBank data
    """
    all_locustag_list = []
    for genome_id in list(df_gene_presence_locustag.columns):
        genbank_file_path = gbk_folder / f"{genome_id}.gbk"
        for record in SeqIO.parse(genbank_file_path, "genbank"):
            for feature in record.features:
                genome_data_list = []
                tag = feature.qualifiers.get("locus_tag")
                if tag:
                    genome_data_list.append(tag[0])  # Locus tag
                    genome_data_list.append(genome_id)  # Genome ID
                    if "product" in feature.qualifiers.keys():
                        genome_data_list.append(
                            feature.qualifiers["product"][0]
                        )  # Prokka annotation
                    else:
                        genome_data_list.append("")
                    genome_data_list.append(
                        int(feature.location.start)
                    )  # Start position
                    genome_data_list.append(int(feature.location.end))  # End position
                    genome_data_list.append(
                        str(feature.extract(record.seq))
                    )  # Nucleotide Seq
                    if "translation" in feature.qualifiers.keys():
                        genome_data_list.append(
                            feature.qualifiers["translation"][0]
                        )  # Amino Acid Seq
                    else:
                        genome_data_list.append("")
                    all_locustag_list.append(genome_data_list)
    all_locustag_df = pd.DataFrame(
        all_locustag_list,
        columns=[
            "Locus_Tag",
            "Genome_ID",
            "Prokka_Annotation",
            "Start_Position",
            "End_Position",
            "Nucleotide_Seq",
            "Amino_Acid_Seq",
        ],
    )
    all_locustag_df.index = all_locustag_df["Locus_Tag"]
    return all_locustag_df


def get_pan_genes(pangene_summary_path):
    """
    Get list of pan genes.

    Parameters:
    pangene_summary_path (Path): Path to pangenome summary file

    Returns:
    list: List of pan genes
    """
    gene_class_table = pd.read_csv(pangene_summary_path, index_col=0)
    pan_gene_list = list(
        gene_class_table.loc[gene_class_table["pangenome_class_2"], :].index
    )
    return pan_gene_list


def process_genes(
    pan_gene_list, df_gene_presence_locustag, all_locustag_df, output_folder
):
    """
    Process genes and save output files.

    Parameters:
    pan_gene_list (list): List of pan genes
    df_gene_presence_locustag (DataFrame): DataFrame with gene presence information
    all_locustag_df (DataFrame): DataFrame with all locus tags
    output_folder (Path): Path to output folder
    """
    for gene_id in [str(i).replace("/", "") for i in pan_gene_list]:
        logging.info(f"   Processing gene: {gene_id}")
        alignment_dir_path = output_folder / "pangenome_alignments" / gene_id / "input"
        alignment_dir_path.mkdir(exist_ok=True, parents=True)
        Neu_fasta_filename = alignment_dir_path / "pangenes.fna"
        AA_fasta_filename = alignment_dir_path / "pangenes.faa"
        gene_locustag = []
        for locus_tag_str in df_gene_presence_locustag.loc[gene_id, :].dropna():
            gene_locustag.extend(
                locus_tag_str.split("\t") if "\t" in locus_tag_str else [locus_tag_str]
            )
        nucleotide_records = []
        amino_acid_records = []
        gene_locustag_set = set(gene_locustag)
        filtered_df = all_locustag_df[all_locustag_df.index.isin(gene_locustag_set)]
        nucleotide_records = [
            SeqRecord(
                Seq(row["Nucleotide_Seq"]),
                id=locustag,
                description=f"{row.get('Prokka_Annotation', 'Unknown')} | {row['Genome_ID']}",
            )
            for locustag, row in filtered_df.iterrows()
        ]
        amino_acid_records = [
            SeqRecord(
                Seq(row["Amino_Acid_Seq"]),
                id=locustag,
                description=f"{row.get('Prokka_Annotation', 'Unknown')} | {row['Genome_ID']}",
            )
            for locustag, row in filtered_df.iterrows()
        ]
        with open(Neu_fasta_filename, "w") as fasta_n_file:
            SeqIO.write(nucleotide_records, fasta_n_file, "fasta")
        with open(AA_fasta_filename, "w") as fasta_aa_file:
            SeqIO.write(amino_acid_records, fasta_aa_file, "fasta")
        logging.info(f"   Finished gene: {gene_id}")


def main():
    parser = argparse.ArgumentParser(description="Process some files.")
    parser.add_argument(
        "--roary_path", type=str, required=True, help="Path to Roary output"
    )
    parser.add_argument(
        "--gbk_folder", type=str, required=True, help="Folder containing GenBank files"
    )
    parser.add_argument(
        "--pangene_summary_path",
        type=str,
        required=True,
        help="Path to pangenome summary file",
    )
    parser.add_argument(
        "--output_folder", type=str, required=True, help="Folder to save output files"
    )

    args = parser.parse_args()

    roary_path = Path(args.roary_path)
    gbk_folder = Path(args.gbk_folder)
    pangene_summary_path = Path(args.pangene_summary_path)
    output_folder = Path(args.output_folder)

    df_gene_presence_binary, df_gene_presence_locustag = load_data(roary_path)
    all_locustag_df = parse_genbank_files(df_gene_presence_locustag, gbk_folder)
    pan_gene_list = get_pan_genes(pangene_summary_path)
    process_genes(
        pan_gene_list, df_gene_presence_locustag, all_locustag_df, output_folder
    )


if __name__ == "__main__":
    main()
