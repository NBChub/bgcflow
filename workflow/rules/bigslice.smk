rule install_bigslice:
    output:
        "resources/bigslice.txt"
    conda:
        "../envs/bigslice.yaml"
    shell:
        """
        download_bigslice_hmmdb
        bigslice --version . > {output}
        """

rule big_slice_prep:
    input:
        #antismash_dir = "data/interim/antismash/",
        antismash = expand("data/interim/antismash/{strains}/", strains = STRAINS)
    output:
        table = "data/interim/antismash/datasets.tsv",
        taxonomy_placement = "data/interim/antismash/placement_list.txt",
    run:
        dataset = {"# Dataset name":["all"], 
               "Path to folder":["."], 
               "Path to taxonomy":["all_taxonomy.tsv"], 
               "Description":["all_samples"]}
        df = pd.DataFrame.from_dict(dataset)
        df.to_csv(output.table, sep="\t", index=False)

        # ncbi placement for bigslice
        custom_placement = df_samples[df_samples.source.eq("custom")].closest_placement_reference.tolist()
        ncbi_placement = df_samples[df_samples.source.eq("ncbi")].genome_id.tolist()
        placement_tax = custom_placement+ncbi_placement

        with open(output.taxonomy_placement, "w") as out:
            for acc in placement_tax:
                out.write(acc + "\n")
            out.close()    

rule fetch_gtdb_taxonomy:
    input: 
        taxonomy_placement = "data/interim/antismash/placement_list.txt"
    output:
        taxonomy = "data/interim/antismash/all_taxonomy_raw.tsv",
    shell:
        """
        wget https://raw.githubusercontent.com/medema-group/bigslice/master/misc/assign_gtdb_taxonomy/fetch_taxonomy_from_api.py -P workflow/scripts/ -nc
        python workflow/scripts/fetch_taxonomy_from_api.py {input.taxonomy_placement} {output.taxonomy}
        """

rule fix_gtdb_taxonomy:
    input: 
        taxonomy = "data/interim/antismash/all_taxonomy_raw.tsv",
    output:
        taxonomy = "data/interim/antismash/all_taxonomy.tsv",
    run:
        df_tax = pd.read_csv("data/interim/antismash/all_taxonomy_raw.tsv", sep="\t")
        df_tax = df_tax.set_index("#Genome folder", drop=False)

        dict_tax = {i : str(df_samples.loc[i, "closest_placement_reference"])+"/" for i in df_samples[df_samples.source.eq("custom")].genome_id}
        dict_tax.update({i : str(i)+"/" for i in df_samples[df_samples.source.eq("ncbi")].genome_id})

        container = []
        for i in df_samples.genome_id:
            line = df_tax.loc[dict_tax[i], :].to_dict()
            line.update({'#Genome folder': f"{i}/"})
            container.append(line)
            
        df = pd.DataFrame(container)
        df.to_csv(output.taxonomy, sep="\t", index=False)

rule bigslice:
    input: 
        ref = "resources/bigslice.txt",
        antismash = expand("data/interim/antismash/{strains}/", strains = STRAINS),
        taxonomy = "data/interim/antismash/all_taxonomy.tsv",
    output:
        folder = directory("data/interim/bigslice/all/")
    conda:
        "../envs/bigslice.yaml"
    threads: 12
    shell:
        """
        bigslice -i data/interim/antismash/ {output.folder}
        """