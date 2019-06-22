library(data.table)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(magrittr)
library(seqinr)
library(diffloop)

# Import data
peaks <- fread("../data/GSE123576_mousebrain_peaks_revision.bed.gz", col.names = c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
tss <- fread("../data/mm10.refGene.TSS.bed", col.names = c("chr", "start", "end", "extra")) %>% makeGRangesFromDataFrame()

# Define sequences for peaks 
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks)
names <- paste0("peak_", as.character(1:length(peaks)))

# Export to fasta
output_file <- "../output/mouse_brain_seqs.fasta"
write.fasta(as.list(as.character(seqs)), names, output_file,
            open = "w", nbchar = 60, as.string = FALSE)
system(paste0("gzip ", output_file))

# Define peaks that overlap tss +/- 1kb
promoters <- padGRanges(tss, 1000)
overlaps <- findOverlaps(peaks, promoters)

# How many? about 27 k which is reasonable
overlaps %>% queryHits() %>% unique() %>% length()

# Export to file
promoter_peak_ids <- names[overlaps %>% queryHits() %>% unique()]
write.table(data.frame(promoter_peak_ids), file = "../output/mouse_brain_promoter_peaks.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)