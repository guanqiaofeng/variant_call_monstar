process XML_VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/pandas' :
        'biocontainers/pandas:2.2.1' }"
        // 'docker.io/gfeng2023/pandas-pyfaidx-image:latest' }"

    input:
    tuple val(meta), path(xml)
    path (hg19_fa)
    path (hg19_fai)

    output:
    tuple val(meta), path ("*.short_variant.vcf"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    shortvariant.py \
        -i ${xml} \
        -r ${hg19_fai} \
        -o ${prefix}.short_variant.vcf
    """
}
//
    // rearrangement.py \\
    //     -i ${xml} \\
    //     -r ${hg19_fa} \\
    //     -r2 ${hg19_fai} \\
        // -o ${prefix}.rearrangement.vcf
