/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_REF } from '../modules/local/download_ref/main'
include { XML_VCF } from '../modules/local/xml_vcf/main'
include { PICARD_LIFTOVERVCF as  PICARD_LIFTOVERVCF_SV } from '../modules/nf-core/picard/liftovervcf/main'
include { PICARD_LIFTOVERVCF as PICARD_LIFTOVERVCF_RA } from '../modules/nf-core/picard/liftovervcf/main'
// include { sanityCheck } from '../modules/local/payload/main'
include { sanityCheck } from '../modules/local/prep/metadata/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_SV } from '../modules/local/payload/main'
include { PAYLOAD_VARIANT_CALL as  PAYLOAD_VARIANT_CALL_RA } from '../modules/local/payload/main'
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

    // SanityCheck
    // sanityCheck(
    //     file(params.experiment_info_tsv),
    //     params.api_token,
    //     params.song_url,
    //     params.clinical_url,
    //     params.skip_duplicate_check
    // )
    // ch_versions = ch_versions.mix(sanityCheck.out.versions)

    // XML to VCF conversion

    xml_ch = Channel.fromPath(params.xml)
                .map { path -> [ [id: '2001205343'], path ] }

    XML_VCF (
        xml_ch,
        Channel.fromPath(params.hg19_ref_fa),
        Channel.fromPath(params.hg19_ref_fai)
    )
    ch_versions = ch_versions.mix(XML_VCF.out.versions)

    // VCF lift over
    hg38_ref_ch = Channel.fromPath(params.hg38_ref_fa)
                            .map{ path -> [ [id: 'fasta'], path ] }

    hg38_ref_dict = Channel.fromPath(params.hg38_ref_dict)
                            .map{ path -> [ [id: 'dict'], path ] }

    hg19_to_hg38_chain_ch = Channel.fromPath(params.hg19_to_hg38_chain)
                            .map{ path -> [ [id: 'chain'], path ] }

    // Short Vraint
    // lift over
    PICARD_LIFTOVERVCF_SV (
        XML_VCF.out.short_variant_vcf,
        hg38_ref_dict,
        hg38_ref_ch,
        hg19_to_hg38_chain_ch
    )
    ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_SV.out.versions)

    //Payload Generation
    sanity_ch = Channel.fromPath(params.sanity)

    PICARD_LIFTOVERVCF_SV.out.vcf_lifted
    .combine(PICARD_LIFTOVERVCF_SV.out.vcf_lifted_index)
    .combine(sanity_ch)
    .map{
        metaA, vcf, metaB, index, metadata_analysis ->
        [
            metaA, [vcf, index], metadata_analysis
        ]
    }.set{vcf_and_index_sv}

    // PICARD_LIFTOVERVCF_SV.out.vcf_lifted
    // .combine(PICARD_LIFTOVERVCF_SV.out.vcf_lifted_index)
    // .combine(sanityCheck.out.updated_experiment_info_tsv)
    // .map{
    //     metaA, vcf, metaB, index, metadata_analysis ->
    //     [
    //         metaA, [vcf, index], metadata_analysis
    //     ]
    // }.set{vcf_and_index}

    PAYLOAD_VARIANT_CALL_SV (
        vcf_and_index_sv,
        Channel.empty()
        .mix(XML_VCF.out.versions)
        .mix(PICARD_LIFTOVERVCF_SV.out.versions)
        .collectFile(name: 'collated_versions.yml')
    )

    // Upload


    // Rearrangement
    // lift over
    PICARD_LIFTOVERVCF_RA (
        XML_VCF.out.rearrangement_vcf,
        hg38_ref_dict,
        hg38_ref_ch,
        hg19_to_hg38_chain_ch
    )
    ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF_RA.out.versions)

    //Payload Generation
    sanity_ch = Channel.fromPath(params.sanity)

    PICARD_LIFTOVERVCF_RA.out.vcf_lifted
    .combine(PICARD_LIFTOVERVCF_RA.out.vcf_lifted_index)
    .combine(sanity_ch)
    .map{
        metaA, vcf, metaB, index, metadata_analysis ->
        [
            metaA, [vcf, index], metadata_analysis
        ]
    }.set{vcf_and_index_ra}

    // PICARD_LIFTOVERVCF_SV.out.vcf_lifted
    // .combine(PICARD_LIFTOVERVCF_SV.out.vcf_lifted_index)
    // .combine(sanityCheck.out.updated_experiment_info_tsv)
    // .map{
    //     metaA, vcf, metaB, index, metadata_analysis ->
    //     [
    //         metaA, [vcf, index], metadata_analysis
    //     ]
    // }.set{vcf_and_index}

    PAYLOAD_VARIANT_CALL_RA (
        vcf_and_index_ra,
        Channel.empty()
        .mix(XML_VCF.out.versions)
        .mix(PICARD_LIFTOVERVCF_RA.out.versions)
        .collectFile(name: 'collated_versions.yml')
    )

    // Payload generation
    // PAYLOAD_VARIANT_CALL (
    //     LIFT_OVER.out.vcf_lifted
    //     molecular_mapping_tsv
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
