/*
    Assemble inputs per mode and run bgc-quast.
    Modes: compare-tools, compare-samples, compare-to-reference (+ auto = TODO).
    QUAST runs ONLY in compare-to-reference, inside this subworkflow, UNLESS the user
    supplies their own QUAST output dir via --bgc_quast_quastdir (then QUAST is skipped).

    Contract with the caller (#5):
      - Query prediction channels are keyed by the query sample id (meta.id).
      - Reference prediction/genome channels are ALSO keyed by the matching query
        sample id (#5 fans the single reference out to each query). Empty in the other modes.

    Output folder naming:
      - Every BGCQUAST run carries meta.leaf, the final folder name used by publishDir:
          compare-samples      -> tool      (antiSMASH | DeepBGC | GECCO)
          compare-tools        -> sample id (e.g. sample_1)
          compare-to-reference -> sample id / tool  (3 tool reports per sample, kept apart)
      - publishDir then writes to: <outdir>/bgc_quast/<mode>/<meta.leaf>/
*/

include { QUAST    } from '../../modules/nf-core/quast/main'
include { BGCQUAST } from '../../modules/local/bgcquast'

workflow BGCQUAST_COMPARISON {
    take:
    antismash_json     // [ meta, json ]        query
    deepbgc_tsv        // [ meta, tsv  ]        query (optional per sample)
    gecco_clusters     // [ meta, tsv  ]        query (optional per sample)
    genomes            // [ meta, fasta ]       query contigs (for --genome and QUAST consensus)
    ref_antismash_json // [ meta, json ]        reference, keyed by query id (empty unless compare-to-reference)
    ref_deepbgc_tsv    // [ meta, tsv  ]        reference
    ref_gecco_clusters // [ meta, tsv  ]        reference
    ref_genome         // [ meta, fasta ]       reference genome, keyed by query id
    ref_name           // val: reference display name (bgc-quast --ref-name); empty unless compare-to-reference

    main:
    ch_versions    = Channel.empty()
    ch_bgcquast_in = Channel.empty()
    def mode       = params.bgc_quast_mode

    // Proper-case folder names for the three tools (used in leaf paths).
    def proper = [antismash: 'antiSMASH', deepbgc: 'DeepBGC', gecco: 'GECCO']

    if (mode == 'compare-tools') {
        // Per sample: gather the present tool outputs into one ordered list -> one run per sample.
        def tool_order = ['antismash', 'deepbgc', 'gecco']

        ch_bgcquast_in = antismash_json.map { meta, f -> [meta, 'antismash', f] }
            .mix(deepbgc_tsv.map    { meta, f -> [meta, 'deepbgc', f] })
            .mix(gecco_clusters.map { meta, f -> [meta, 'gecco', f] })
            .groupTuple(by: 0)
            .map { meta, tools, files ->
                // Keep a stable tool order (antismash, deepbgc, gecco) for the report columns.
                def idx           = (0..<tools.size()).toList().sort { tool_order.indexOf(tools[it]) }
                def ordered_files = idx.collect { files[it] }
                [meta, ordered_files]
            }
            .join(genomes)
            .map { meta, files, genome ->
                // mining_results, genome, quast_dir, reference_mining, reference_genome
                // No --names here: bgc-quast auto-labels each column by detected tool.
                // leaf = sample id
                [meta + [leaf: "${meta.id}"], files, genome, [], [], []]
            }
    }
    else if (mode == 'compare-samples') {
        // Per tool: collect that tool's output across all samples -> one run per tool.
        def by_tool = { ch, tool ->
            ch.join(genomes).map { meta, f, g -> [tool, meta.id, f, g] }
        }

        ch_bgcquast_in = by_tool(antismash_json, 'antismash')
            .mix(by_tool(deepbgc_tsv, 'deepbgc'))
            .mix(by_tool(gecco_clusters, 'gecco'))
            .groupTuple(by: 0)
            .map { tool, ids, files, gens ->
                // ids[i] <-> files[i] <-> gens[i] stay aligned through groupTuple
                // leaf = proper-case tool name
                [[id: "compare_samples_${tool}", bgcquast_names: ids.join(','), leaf: proper[tool]], files, gens, [], [], []]
            }
    }
    else if (mode == 'compare-to-reference') {
        // QUAST per query (its contigs vs the reference) -- unless the user supplies their own dir.
        if (params.bgc_quast_quastdir) {
            def user_quast = file(params.bgc_quast_quastdir, checkIfExists: true)
            ch_quast_dir = genomes.map { meta, g -> [meta, user_quast] } // same dir attached to each query
        }
        else {
            ch_quast_in = genomes.join(ref_genome) // [meta, query_contigs, ref_genome]

            QUAST(
                ch_quast_in.map { meta, q, _r -> [meta, q] },
                ch_quast_in.map { meta, _q, r -> [meta, r] },
                ch_quast_in.map { meta, _q, _r -> [meta, []] },
            )
            // QUAST version flows via its `topic: versions` channel; not mixed here.
            ch_quast_dir = QUAST.out.results // [meta, dir]
        }

        // Per (sample, tool): join query pred + reference pred + QUAST dir + genomes, all on sample id.
        // leaf = "<sample id>/<tool>" so the 3 tool reports per sample go to separate folders.
        def per_tool = { qch, rch, tool ->
            qch.join(rch)              // [meta, qfile, rfile]
                .join(ch_quast_dir)    // + qdir
                .join(genomes)         // + genome
                .join(ref_genome)      // + ref_genome
                .combine(ref_name)     // + rid (single value)
                .map { meta, qfile, rfile, qdir, genome, refgenome, rid ->
                    [meta + [bgcquast_names: meta.id, ref_name: rid, leaf: "${meta.id}/${proper[tool]}"], [qfile], genome, qdir, rfile, refgenome]
                }
        }

        ch_bgcquast_in = per_tool(antismash_json, ref_antismash_json, 'antismash')
            .mix(per_tool(deepbgc_tsv, ref_deepbgc_tsv, 'deepbgc'))
            .mix(per_tool(gecco_clusters, ref_gecco_clusters, 'gecco'))
    }
    else {
        // TODO: auto mode undecided. Decide how to infer the mode from inputs, then implement.
        error("[bgc_quast_ppl] bgc_quast_mode='${mode}' is not supported yet. Use compare-tools, compare-samples, or compare-to-reference.")
    }

    BGCQUAST(ch_bgcquast_in)
    ch_versions = ch_versions.mix(BGCQUAST.out.versions)

    emit:
    results  = BGCQUAST.out.results // [ meta, files ]
    tsv      = BGCQUAST.out.tsv     // [ meta, report.tsv ]
    versions = ch_versions          // [ path(versions.yml) ]
}
