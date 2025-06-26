def chunkPrepareChannel(ch_chunks, tool) {
    if(tool == "glimpse"){
        ch_chunks = ch_chunks.map { chr, txt -> [chr, file(txt)]}
                .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
                .map { meta, it -> [meta, it["RegionIn"], it["RegionOut"]]}
    } else if(tool == "quilt") {
        ch_chunks = ch_chunks.map { chr, txt -> [chr, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { meta, it ->
                def startEnd = it["RegionIn"].split(':')[1].split('-')
                [meta, meta.chr, startEnd[0], startEnd[1]]
            }
    } else {
        error "Only 'glimpse' and 'quilt' output format are supported. Got ${tool}"
    }
    return ch_chunks // channel:   [ [meta], regionstart, regionend ]
}
