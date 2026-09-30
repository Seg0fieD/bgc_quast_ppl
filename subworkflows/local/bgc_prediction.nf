/*
    BGC prediction with antiSMASH, DeepBGC and GECCO, each unless skipped;
    the comparison modes, QUAST and bgc-quast run in a separate subworkflow.
*/

include { UNTAR as UNTAR_ANTISMASHDB           } from '../../modules/nf-core/untar/main'
include { ANTISMASH_ANTISMASHDOWNLOADDATABASES } from '../../modules/nf-core/antismash/antismashdownloaddatabases/main'
include { ANTISMASH_ANTISMASH                  } from '../../modules/nf-core/antismash/antismash/main'
include { DEEPBGC_DOWNLOAD                     } from '../../modules/nf-core/deepbgc/download/main'
include { DEEPBGC_PIPELINE                     } from '../../modules/nf-core/deepbgc/pipeline/main'
include { GECCO_RUN                            } from '../../modules/nf-core/gecco/run/main'

workflow BGC_PREDICTION {
    take:
    fastas // tuple val(meta), path(fasta): unused, kept for argument order
    faas   // tuple val(meta), path(faa): unused, kept for argument order
    gbks   // tuple val(meta), path(gbk): input to all three tools

    main:
    ch_versions       = Channel.empty()
    ch_antismash_json = Channel.empty()
    ch_antismash_gbk  = Channel.empty()
    ch_deepbgc_tsv    = Channel.empty()
    ch_deepbgc_gbk    = Channel.empty()
    ch_gecco_clusters = Channel.empty()
    ch_gecco_gbk      = Channel.empty()

    if (!params.bgc_skip_antismash) {
        // User-provided DB path, else download it.
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
        ch_antismash_gbk  = ANTISMASH_ANTISMASH.out.gbk_results
    }

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
        // One multi-record GBK per sample; split per BGC before BiG-SCAPE.
        ch_deepbgc_gbk = DEEPBGC_PIPELINE.out.bgc_gbk
    }

    if (!params.bgc_skip_gecco) {
        ch_gecco_input = gbks
            .groupTuple()
            .map { meta, gbk -> [meta, gbk, []] }

        GECCO_RUN(ch_gecco_input, [])
        ch_versions       = ch_versions.mix(GECCO_RUN.out.versions)
        ch_gecco_clusters = GECCO_RUN.out.clusters
        // Already one file per cluster, so no split step is needed.
        ch_gecco_gbk      = GECCO_RUN.out.gbk
    }

    emit:
    // A sample with no BGCs is absent here, except from antismash_json.
    versions       = ch_versions       // [ versions.yml ]
    antismash_json = ch_antismash_json // [ meta, *.json ]
    antismash_gbk  = ch_antismash_gbk  // [ meta, [ *region*.gbk ] ]
    deepbgc_tsv    = ch_deepbgc_tsv    // [ meta, *.bgc.tsv ]
    deepbgc_gbk    = ch_deepbgc_gbk    // [ meta, *.bgc.gbk ]
    gecco_clusters = ch_gecco_clusters // [ meta, *.clusters.tsv ]
    gecco_gbk      = ch_gecco_gbk      // [ meta, [ *_cluster_*.gbk ] ]
}
