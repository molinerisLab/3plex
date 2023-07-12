
if config['species']:
    config['genome_fasta'] = "{gencode_dir_prefix}/{species}/{gencode_version}/{genome_version}.primary_assembly.genome.fa".format(
        gencode_dir_prefix=config['gencode_dir_prefix'],
        species=config['species'],
        gencode_version=config['species_gencode_version_map'][config['species']],
        genome_version=config['species_genome_version_map'][config['species']]
        )

    config['shuffle_blacklist'] = "{gencode_dir_prefix}/{species}/{gencode_version}/shuffle_blacklist.bed".format(
        gencode_dir_prefix=config['gencode_dir_prefix'],
        species=config['species'],
        gencode_version=config['species_gencode_version_map'][config['species']],
        )

    config['chrom_info'] = "{gencode_dir_prefix}/{species}/{gencode_version}/chrom.info.no_alt".format(
        gencode_dir_prefix=config['gencode_dir_prefix'],
        species=config['species'],
        gencode_version=config['species_gencode_version_map'][config['species']],
        )
else:
    config['genome_fasta'] = ""
    config['shuffle_blacklist'] = ""
    config['chrom_info'] = ""


if config['dsDNA_predefined']:
    config['dsDNA_predefined_fa']="{gencode_dir_prefix}/{species}/{gencode_version}/{dsDNA_predefined}.fa".format(
        gencode_dir_prefix=config['gencode_dir_prefix'],
        species=config['species'],
        gencode_version=config['species_gencode_version_map'][config['species']],
        dsDNA_predefined=config['dsDNA_predefined']
    )
    config['dsDNA_predefined_bed']="{gencode_dir_prefix}/{species}/{gencode_version}/{dsDNA_predefined}.bed".format(
        gencode_dir_prefix=config['gencode_dir_prefix'],
        species=config['species'],
        gencode_version=config['species_gencode_version_map'][config['species']],
        dsDNA_predefined=config['dsDNA_predefined']
    )
else:
    config['dsDNA_predefined_fa'] = ""
    config['dsDNA_predefined_bed'] = ""
