/*
    Run BGC prediction tools (antiSMASH, DeepBGC, GECCO).
    Prediction only: no modes, no QUAST, no bgc-quast. Those live in the comparison subworkflow.
*/

include { UNTAR as UNTAR_ANTISMASHDB           } from '../../modules/nf-core/untar/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES } from '../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { ANTISMASH_ANTISMASH                  } from '../../modules/nf-core/antismash/antismash/main'
include { DEEPBGC_DOWNLOAD                     } from '../../modules/nf-core/deepbgc/download/main'
include { DEEPBGC_PIPELINE                     } from '../../modules/nf-core/deepbgc/pipeline/main'
include { GECCO_RUN                            } from '../../modules/nf-core/gecco/run/main'

workflow BGC_PREDICTION {
    take:
    fastas // tuple val(meta), path(fasta) -- not consumed here; kept to match the caller (#5)
    faas   // tuple val(meta), path(faa)   -- not consumed here; kept to match the caller (#5)
    gbks   // tuple val(meta), path(gbk)   -- the input all three tools read

    main:
    ch_versions       = Channel.empty()
    ch_antismash_json = Channel.empty()
    ch_deepbgc_tsv    = Channel.empty()
    ch_gecco_clusters = Channel.empty()

    // ANTISMASH
    if (!params.bgc_skip_antismash) {
        // User-supplied DB (gz tarball or directory) if given, else download it.
        if (params.bgc_antismash_db && file(params.bgc_antismash_db, checkIfExists: true).extension == 'gz') {
            UNTAR_ANTISMASHDB([[id: 'antismashdb'], file(params.bgc_antismash_db, checkIfExists: true)])
            ch_antismash_databases = UNTAR_ANTISMASHDB.out.untar.map { _meta, dir -> [dir] }
            ch_versions = ch_versions.mix(UNTAR_ANTISMASHDB.out.versions)
        }
        else if (params.bgc_antismash_db && file(params.bgc_antismash_db, checkIfExists: true).isDirectory()) {
            ch_antismash_databases = Channel.fromPath(params.bgc_antismash_db, checkIfExists: true).first()
        }
        else {
            ANTISMASH_ANTISMASHDOWNLOADDATABASES()
            ch_versions = ch_versions.mix(ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.versions)
            ch_antismash_databases = ANTISMASH_ANTISMASHDOWNLOADDATABASES.out.database
        }

        ANTISMASH_ANTISMASH(gbks, ch_antismash_databases, [])
        ch_versions       = ch_versions.mix(ANTISMASH_ANTISMASH.out.versions)
        ch_antismash_json = ANTISMASH_ANTISMASH.out.json_results
    }

    // DEEPBGC
    if (!params.bgc_skip_deepbgc) {
        if (params.bgc_deepbgc_db) {
            ch_deepbgc_database = Channel.fromPath(params.bgc_deepbgc_db, checkIfExists: true).first()
        }
        else {
            DEEPBGC_DOWNLOAD()
            ch_deepbgc_database = DEEPBGC_DOWNLOAD.out.db
            ch_versions = ch_versions.mix(DEEPBGC_DOWNLOAD.out.versions)
        }

        DEEPBGC_PIPELINE(gbks, ch_deepbgc_database)
        ch_versions    = ch_versions.mix(DEEPBGC_PIPELINE.out.versions)
        ch_deepbgc_tsv = DEEPBGC_PIPELINE.out.bgc_tsv
    }

    // GECCO
    if (!params.bgc_skip_gecco) {
        // GECCO reads the GBK too (same as antiSMASH/DeepBGC). hmm = [], model_dir = [] (bundled model).
        ch_gecco_input = gbks
            .groupTuple()
            .map { meta, gbk -> [meta, gbk, []] }

        GECCO_RUN(ch_gecco_input, [])
        ch_versions       = ch_versions.mix(GECCO_RUN.out.versions)
        ch_gecco_clusters = GECCO_RUN.out.clusters
    }

    emit:
    versions       = ch_versions          // channel: [ path(versions.yml) ]
    antismash_json = ch_antismash_json    // channel: [ val(meta), path(*.json) ]
    deepbgc_tsv    = ch_deepbgc_tsv        // channel: [ val(meta), path(*.bgc.tsv) ]   (optional per sample)
    gecco_clusters = ch_gecco_clusters    // channel: [ val(meta), path(*.clusters.tsv) ] (optional per sample)
}
