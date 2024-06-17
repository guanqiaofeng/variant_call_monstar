/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XML_VCF } from '../modules/local/xml_vcf/main'
include { PICARD_LIFTOVERVCF } from '../modules/nf-core/picard/liftovervcf/main'
include { PAYLOAD_VARIANT_CALL } from '../modules/local/payload/variantcall/main'
include { SONG_SCORE_UPLOAD } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCALL {

    take:
    xml

    main:

    ch_versions = Channel.empty()

    // XML to VCF conversion
    XML_VCF (
        xml
    )
    ch_versions = ch_versions.mix(XML_VCF.out.versions)

    // VCF lift over
    PICARD_LIFTOVERVCF (
        XML_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(LIFT_OVER.out.versions)

    // Payload generation
    PAYLOAD_VARIANT_CALL (
        LIFT_OVER.out.updatedvcf
    )
    ch_versions = ch_versions.mix(PAYLOAD_VARIANT_CALL.out.versions)

    // Upload
    SONG_SCORE_UPLOAD (
        PAYLOAD_ALIGNMENT_H.out.payload_files
    ) // [val(meta), path("*.payload.json"), path VCF]
    ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD.out.versions)

    emit:

     versions = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
