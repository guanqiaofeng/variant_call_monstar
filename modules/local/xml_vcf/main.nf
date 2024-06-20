process XML_VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/pandas' :
        'biocontainers/pandas:2.2.1' }"


    input:
    tuple val(meta), path(xml)

    output:
    // tuple val(meta), path "*.vcf", emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    main.py \\
        -i ${xml} \\
        -o ${prefix}.short_variant.vcf
    """
}
