/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_REF } from '../modules/local/download_ref/main'
include { XML_VCF } from '../modules/local/xml_vcf/main'
include { PICARD_LIFTOVERVCF } from '../modules/nf-core/picard/liftovervcf/main'
// include { PAYLOAD_VARIANT_CALL } from '../modules/local/payload/variantcall/main'
// include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCALL {

    // take:
    // xml

    // main:

    ch_versions = Channel.empty()

    // Download Reference fasta, fai, liftover_chain files
    // DOWNLOAD_REF (
    //     params.fasta_url,
    //     params.fasta_checksum,
    //     params.fai_url,
    //     params.fai_checksum,
    //     params.chain_url,
    //     params.chain_checksum
    // )

    // XML to VCF conversion

    xml_ch = Channel.fromPath(params.xml)
                .map { path -> [ [id: '2001205343'], path ] }

    XML_VCF (
        xml_ch,
        Channel.fromPath(params.hg19_ref_fa),
        Channel.fromPath(params.hg19_ref_fai)
    )
    // ch_versions = ch_versions.mix(XML_VCF.out.versions)

    // VCF lift over
    hg38_ref_ch = Channel.fromPath(params.hg38_ref_fa)
                            .map{ path -> [ [id: 'fasta'], path ] }

    hg38_ref_dict = Channel.fromPath(params.hg38_ref_dict)
                            .map{ path -> [ [id: 'dict'], path ] }

    hg19_to_hg38_chain_ch = Channel.fromPath(params.hg19_to_hg38_chain)
                            .map{ path -> [ [id: 'chain'], path ] }

    PICARD_LIFTOVERVCF (
        XML_VCF.out.vcf,
        hg38_ref_dict,
        hg38_ref_ch,
        hg19_to_hg38_chain_ch
    )
    ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF.out.versions)

    // // Payload generation
    // PAYLOAD_VARIANT_CALL (
    //     LIFT_OVER.out.updatedvcf
    // )
    // ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL.out.versions)

    // // Upload
    // SONG_SCORE_UPLOAD (
    //     PAYLOAD_ALIGNMENT_H.out.payload_files
    // ) // [val(meta), path("*.payload.json"), path VCF]
    // ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD.out.versions)

    emit:

     versions = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
