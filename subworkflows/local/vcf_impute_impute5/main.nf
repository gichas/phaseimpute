include { IMPUTE5_CONVERTREF } from '../../../modules/nf-core/impute5/convertref/main'
include { IMPUTE5_CHUNK      } from '../../../modules/nf-core/impute5/chunk/main'
include { IMPUTE5_IMPUTE     } from '../../../modules/nf-core/impute5/impute/main'

workflow VCF_IMPUTE_IMPUTE5 {

    take:
    ch_input  // channel: [ [id, chr], vcf, tbi ]
    ch_panel  // channel: [ [id, chr], vcf, tbi ]  
    ch_map    // channel: [ [chr], map ]

    main:
    ch_versions = Channel.empty()

    // Convert reference panel to XCF format
    ch_panel_with_region = ch_panel.map { meta, vcf, tbi ->
        tuple(meta, vcf, tbi, meta.chr)
    }
    
    IMPUTE5_CONVERTREF(ch_panel_with_region)
    ch_versions = ch_versions.mix(IMPUTE5_CONVERTREF.out.versions)

    // Generate chunks
    ch_chunk_input = ch_input
        .map { meta, vcf, tbi -> [meta.chr, meta, vcf, tbi] }
        .combine(
            IMPUTE5_CONVERTREF.out.xcf_file.map { meta, xcf, idx, bin, fam -> 
                [meta.chr, xcf, idx, bin, fam] 
            }, 
            by: 0
        )
        .map { chr, target_meta, vcf, tbi, ref_xcf, ref_idx, ref_bin, ref_fam ->
            tuple(target_meta, ref_xcf, ref_idx, vcf, tbi, target_meta.chr)
        }

    IMPUTE5_CHUNK(ch_chunk_input)
    ch_versions = ch_versions.mix(IMPUTE5_CHUNK.out.versions)

    // Parse chunks file to create individual chunk jobs (nf-core style)
    ch_parsed_chunks = IMPUTE5_CHUNK.out.chunks
        .flatMap { meta, chunks_file ->
            def chunks = []
            chunks_file.readLines().eachWithIndex { line, index ->
                if (index > 0 && line.trim() && !line.startsWith('#')) {
                    def cols = line.split('\t')
                    if (cols.size() >= 4) {
                        def chunk_meta = meta + [chunk_id: cols[0]]
                        chunks.add([chunk_meta, cols[3]]) // [meta_with_chunk, impute_region]
                    }
                }
            }
            return chunks
        }

    // Prepare imputation input: combine chunks with all required data
    ch_impute_input = ch_parsed_chunks
        .map { chunk_meta, region -> [chunk_meta.chr, chunk_meta, region] }
        .combine(
            ch_input.map { meta, vcf, tbi -> [meta.chr, vcf, tbi] }, 
            by: 0
        )
        .combine(
            IMPUTE5_CONVERTREF.out.xcf_file.map { meta, xcf, idx, bin, fam -> 
                [meta.chr, xcf, idx, bin, fam] 
            }, 
            by: 0
        )
        .combine(
            ch_map.map { meta, map -> [meta.chr, map] }, 
            by: 0
        )
        .map { chr, chunk_meta, region, vcf, tbi, ref_xcf, ref_idx, ref_bin, ref_fam, map ->
            tuple(chunk_meta, vcf, tbi, ref_xcf, ref_idx, ref_bin, ref_fam, map, region)
        }

    // Run imputation on each chunk individually (simple module call)
    IMPUTE5_IMPUTE(ch_impute_input)
    ch_versions = ch_versions.mix(IMPUTE5_IMPUTE.out.versions)

    // Group chunks back by sample/chromosome for output
    ch_imputed_vcf_tbi = IMPUTE5_IMPUTE.out.vcf
        .map { chunk_meta, vcf -> 
            def sample_meta = [id: chunk_meta.id, chr: chunk_meta.chr, tools: "impute5"]
            [sample_meta, vcf]
        }

    emit:
    vcf_tbi  = ch_imputed_vcf_tbi // channel: [ [id, chr, tools], vcf ]
    versions = ch_versions        // channel: [ versions.yml ]
}