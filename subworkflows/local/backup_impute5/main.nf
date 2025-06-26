include { IMPUTE5_CHUNK      } from '../../../modules/nf-core/impute5/chunk/main'
include { IMPUTE5_CONVERTREF } from '../../../modules/nf-core/impute5/convertref/main'
include { IMPUTE5_IMPUTE     } from '../../../modules/nf-core/impute5/impute/main'

workflow VCF_IMPUTE_IMPUTE5 {
    take:
    ch_target_data    // tuple val(meta), path(target_bcf), path(target_index), val(region)
    ch_ref_data       // tuple val(meta), path(ref_bcf), path(ref_index), val(region)
    ch_genetic_map    // tuple val(meta), path(genetic_map), val(region)

    main:
    ch_versions = Channel.empty()

    // Étape 1: Convertir la référence au format XCF
    IMPUTE5_CONVERTREF (
        ch_ref_data
    )
    ch_versions = ch_versions.mix(IMPUTE5_CONVERTREF.out.versions)

    // Étape 2: Créer les chunks pour l'imputation
    ch_chunk_input = ch_target_data
        .join(IMPUTE5_CONVERTREF.out.xcf_file, by: 0)
        .map { meta, target_bcf, target_index, target_region, ref_xcf, ref_index, ref_bin, ref_fam ->
            tuple(meta, ref_xcf, ref_index, ref_bin, ref_fam, target_bcf, target_index, target_region)
        }

    IMPUTE5_CHUNK (
        ch_chunk_input
    )
    ch_versions = ch_versions.mix(IMPUTE5_CHUNK.out.versions)

    // Étape 3: Parser les chunks et créer un channel par chunk
    ch_chunks_parsed = IMPUTE5_CHUNK.out.chunks
        .map { meta, chunks_file ->
            // Lire le fichier chunks et créer une entrée par chunk
            def chunks = file(chunks_file).readLines()
                .findAll { !it.startsWith('#') && it.trim() } // Ignorer commentaires et lignes vides
                .collect { line ->
                    def fields = line.split('\t')
                    return [
                        chunk_id: fields[0],
                        chr: fields[1],
                        buffer_region: fields[2],
                        impute_region: fields[3]
                    ]
                }

            // Retourner une liste de tuples [meta_with_chunk, chunk_info]
            return chunks.collect { chunk ->
                def meta_chunk = meta + [chunk_id: chunk.chunk_id]
                return tuple(meta_chunk, chunk.impute_region)
            }
        }
        .flatten() // Aplatir pour avoir un élément par chunk
        .map { it } // Chaque élément est déjà un tuple(meta_chunk, region)

    // Étape 4: Préparer les données pour chaque chunk
    ch_impute_input = ch_chunks_parsed
        .combine(ch_target_data.map { meta, bcf, idx, region -> tuple(meta.id, bcf, idx) }, by: 0)
        .combine(IMPUTE5_CONVERTREF.out.xcf_file.map { meta, xcf, idx, bin, fam -> tuple(meta.id, xcf, idx, bin, fam) }, by: 0)
        .combine(ch_genetic_map.map { meta, gmap, region -> tuple(meta.id, gmap) }, by: 0)
        .map { meta_id, meta_chunk, chunk_region, target_bcf, target_idx, ref_xcf, ref_idx, ref_bin, ref_fam, gen_map ->
            // Reconstruire le tuple pour IMPUTE5_IMPUTE
            tuple(meta_chunk, target_bcf, target_idx, ref_xcf, ref_idx, ref_bin, ref_fam, gen_map, chunk_region)
        }

    // Étape 6: Effectuer l'imputation par chunk (en parallèle)
    IMPUTE5_IMPUTE (
        ch_impute_input
    )
    ch_versions = ch_versions.mix(IMPUTE5_IMPUTE.out.versions)

    emit:
    // Résultats par chunk (pour analyses parallèles ou ligation manuelle)
    imputed_chunks = IMPUTE5_IMPUTE.out.bcf     // tuple val(meta_chunk), path("*_imputed.bcf")
    imputed_logs   = IMPUTE5_IMPUTE.out.log     // tuple val(meta_chunk), path("*_imputed.log")

    // Métadonnées
    chunks         = IMPUTE5_CHUNK.out.chunks   // tuple val(meta), path("chunks_*.txt")
    xcf_files      = IMPUTE5_CONVERTREF.out.xcf_file // référence convertie
    versions       = ch_versions                // versions de tous les outils
}
