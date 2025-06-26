include { IMPUTE5_CONVERTREF } from '../../../modules/nf-core/impute5/convertref/main'
include { IMPUTE5_CHUNK      } from '../../../modules/nf-core/impute5/chunk/main'
include { IMPUTE5_IMPUTE     } from '../../../modules/nf-core/impute5/impute/main'
include { BCFTOOLS_CONCAT    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CHUNKS } from '../../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FINAL  } from '../../../modules/nf-core/bcftools/index/main'
workflow VCF_IMPUTE_IMPUTE5 {

    take:
    ch_target_data    // tuple val(meta), path(target_bcf), path(target_index), val(region)
    ch_ref_data       // tuple val(meta), path(ref_bcf), path(ref_index), val(region)
    ch_genetic_map    // tuple val(meta), path(genetic_map), val(region)

    main:
    ch_versions = Channel.empty()

    //
    // Convert reference panel to xcf format for IMPUTE5
    //
    IMPUTE5_CONVERTREF (
        ch_ref_data
    )
    ch_versions = ch_versions.mix(IMPUTE5_CONVERTREF.out.versions)

    //
    // Prepare input data for chunking by combining target and converted reference
    //
    ch_chunk_input = ch_target_data
        .combine(IMPUTE5_CONVERTREF.out.xcf_file)
        .filter { target_meta, target_bcf, target_index, target_region,
                  ref_meta, ref_xcf, ref_index, ref_bin, ref_fam ->
            target_meta.chr == ref_meta.chr
        }
        .map { target_meta, target_bcf, target_index, target_region,
               ref_meta, ref_xcf, ref_index, ref_bin, ref_fam ->
            tuple(target_meta, ref_xcf, ref_index, target_bcf, target_index, target_region)
        }

    //
    // Generate optimal chunks for imputation
    //
    IMPUTE5_CHUNK (
        ch_chunk_input
    )
    ch_versions = ch_versions.mix(IMPUTE5_CHUNK.out.versions)

    //
    // Parse chunks file to extract optimal imputation regions
    //
    ch_chunks_parsed = IMPUTE5_CHUNK.out.chunks
        .flatMap { meta, chunks_file ->
            def chunks = []
            try {
                chunks_file.readLines().eachWithIndex { line, index ->
                    if (index > 0 && line.trim() && !line.startsWith('#')) { // Skip header and empty lines
                        def cols = line.split('\t')
                        if (cols.size() >= 6) {
                            def chunk_id = cols[0]
                            def chr = cols[1]
                            def full_region = cols[2]      // FULL_REGION (with buffer)
                            def impute_region = cols[3]    // IMPUTE_REGION (optimal imputation region)
                            def size = cols[4]
                            def target_snp = cols[5]

                            def meta_with_chunk = meta + [
                                chunk_id: chunk_id,
                                full_region: full_region,
                                target_snp: target_snp,
                                size: size
                            ]
                            chunks.add([meta_with_chunk, impute_region, full_region])
                        }
                    }
                }
            } catch (Exception e) {
                log.warn "Failed to parse chunks file for ${meta.id}: ${e.message}"
                // Fallback: use original region
                def target_region = meta.region ?: "chr${meta.chr}:1-50000000"
                chunks.add([meta + [chunk_id: "fallback"], target_region, target_region])
            }
            return chunks
        }

    //
    // Prepare input channels for imputation by combining all required data
    //
    ch_impute_input = ch_chunks_parsed
        .map { meta_chunk, impute_region, full_region ->
            [meta_chunk.chr, meta_chunk, impute_region, full_region]
        }
        .combine(ch_target_data.map { meta, bcf, idx, region ->
            [meta.chr, meta, bcf, idx]
        }, by: 0)
        .combine(IMPUTE5_CONVERTREF.out.xcf_file.map { meta, xcf, idx, bin, fam ->
            [meta.chr, meta, xcf, idx, bin, fam]
        }, by: 0)
        .combine(ch_genetic_map.map { meta, gmap, region ->
            [meta.chr, meta, gmap]
        }, by: 0)
        .map { chr, meta_chunk, impute_region, full_region,
               target_meta, target_bcf, target_idx,
               ref_meta, ref_xcf, ref_idx, ref_bin, ref_fam,
               map_meta, gen_map ->
            tuple(meta_chunk, target_bcf, target_idx, ref_xcf, ref_idx, ref_bin, ref_fam, gen_map, impute_region)
        }

    //
    // Perform imputation on each chunk
    //
    IMPUTE5_IMPUTE (
        ch_impute_input
    )
    ch_versions = ch_versions.mix(IMPUTE5_IMPUTE.out.versions)

    //
    // Index imputed BCF files for concatenation
    //
    BCFTOOLS_INDEX_CHUNKS (
        IMPUTE5_IMPUTE.out.bcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_CHUNKS.out.versions)

    //
    // Prepare input for concatenation by grouping chunks per sample/chromosome
    //
    ch_merge_input = IMPUTE5_IMPUTE.out.bcf
        .join(BCFTOOLS_INDEX_CHUNKS.out.csi)
        .map { meta, bcf, index ->
            def group_key = [id: meta.id, chr: meta.chr]
            [group_key, bcf, index]
        }
        .groupTuple()
        .map { group_meta, bcf_list, index_list ->
            // Add tool name and create final prefix for concatenated file
            def final_meta = group_meta + [
                tools: "impute5",
                ext: [prefix: "${group_meta.id}_chr${group_meta.chr}_imputed"]
            ]
            [final_meta, bcf_list.sort(), index_list.sort()]
        }

    //
    // Concatenate imputed chunks into final chromosome-level files
    //
    BCFTOOLS_CONCAT (
        ch_merge_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    //
    // Use concatenated VCF files with their indices directly from CONCAT
    //
    ch_imputed_vcf_tbi = BCFTOOLS_CONCAT.out.vcf
        .join(BCFTOOLS_CONCAT.out.csi)

    emit:
    vcf_tbi     = ch_imputed_vcf_tbi              // channel: [ [id, chr, tools], vcf, tbi ]
    imputed_log = IMPUTE5_IMPUTE.out.log          // channel: [ [meta], log ] - individual chunk logs
    chunks      = IMPUTE5_CHUNK.out.chunks        // channel: [ [meta], chunks_file ]
    xcf_files   = IMPUTE5_CONVERTREF.out.xcf_file // channel: [ [meta], xcf, index, bin, fam ]
    versions    = ch_versions                     // channel: [ versions.yml ]
}
