rule gtdb_prep:
    input:
        fasta = expand("data/interim/fasta/{strains}.fna", strains = STRAINS)
    output:
        taxonomy_placement = "data/interim/gtdb/placement_list.txt",
    run:
        # NCBI accession or closest NCBI
        custom_placement = df_samples[df_samples.source.eq("custom")].closest_placement_reference.tolist()
        ncbi_placement = df_samples[df_samples.source.eq("ncbi")].genome_id.tolist()
        patric_placement = df_samples[df_samples.source.eq("patric")].closest_placement_reference.tolist()
        placement_tax = custom_placement + ncbi_placement + patric_placement

        with open(output.taxonomy_placement, "w") as out:
            for acc in placement_tax:
                out.write(acc + "\n")
            out.close()    

rule fetch_gtdb_taxonomy:
    input: 
        taxonomy_placement = "data/interim/gtdb/placement_list.txt"
    output:
        taxonomy = "data/interim/gtdb/all_taxonomy_raw.tsv",
    params:
        version = 'R202'    
    conda:
        "../envs/gtdb.yaml"
    shell:
        """
        wget https://raw.githubusercontent.com/medema-group/bigslice/master/misc/assign_gtdb_taxonomy/fetch_taxonomy_from_api.py -P workflow/scripts/ -nc
        python workflow/scripts/fetch_taxonomy_from_api.py {input.taxonomy_placement} {output.taxonomy} --gtdb {params.version}
        """

rule fix_gtdb_taxonomy:
    input: 
        taxonomy_raw = "data/interim/gtdb/all_taxonomy_raw.tsv",
    output:
        taxonomy = "data/interim/gtdb/all_taxonomy.tsv",
        meta = "data/processed/tables/df_gtdb_meta.csv"
    run:
        df_tax = pd.read_csv(input.taxonomy_raw, sep="\t")
        df_tax = df_tax.set_index("#Genome folder", drop=False)

        dict_tax = {i : str(df_samples.loc[i, "closest_placement_reference"])+"/" for i in df_samples[df_samples.source.eq("custom")].genome_id}
        dict_tax.update({i : str(df_samples.loc[i, "closest_placement_reference"])+"/" for i in df_samples[df_samples.source.eq("patric")].genome_id})
        dict_tax.update({i : str(i)+"/" for i in df_samples[df_samples.source.eq("ncbi")].genome_id})

        container = []
        df_meta = pd.DataFrame(index=df_samples.index, columns=df_tax.columns)
        for i in df_samples.genome_id:
            line = df_tax.loc[dict_tax[i], :].to_dict()
            line.update({'#Genome folder': f"{i}/"})
            container.append(line)
            df_meta.loc[i, :] = df_tax.loc[dict_tax[i], :]
            
        df = pd.DataFrame(container)
        df.to_csv(output.taxonomy, sep="\t", index=False)
        df_meta.drop(columns='#Genome folder', inplace=True)
        df_meta.to_csv(output.meta)
