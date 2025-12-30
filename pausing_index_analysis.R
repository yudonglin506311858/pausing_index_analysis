
setwd("/data6/yudonglin/帮忙分析/jw/251204-POLII-CHIP/bowtie2")
library(rtracklayer)
library(GenomicRanges)
library(GenomicAlignments)
library(dplyr)

# 1. 读取注释文件
annotation <- read.table("/data/yudonglin/reference/mouse/mm10/mm10_ucsc.tsv", 
                         sep="\t", skip=1)
colnames(annotation) <- c("chrom", "start", "end", "name", "score", 
                          "strand", "thickStart", "thickEnd", "itemRgb", 
                          "blockCount", "blockSizes", "blockStarts")
head(annotation)
# 2. 计算TSS和TES
annotation$TSS <- ifelse(annotation$strand == "+", 
                         annotation$start, annotation$end)
annotation$TES <- ifelse(annotation$strand == "+", 
                         annotation$end, annotation$start)
annotation$gene_length <- annotation$end - annotation$start

# 3. 过滤短基因
annotation <- annotation[annotation$gene_length >= 400, ]

# 4. 定义promoter和gene body区域
annotation$promoter_start <- annotation$TSS - 200
annotation$promoter_end <- annotation$TSS + 200
annotation$promoter_start <- ifelse(annotation$promoter_start < 0, 
                                    0, annotation$promoter_start)

# 对于正链基因
pos_strand <- annotation$strand == "+"
annotation$genebody_start[pos_strand] <- annotation$TSS[pos_strand] + 400
annotation$genebody_end[pos_strand] <- annotation$TES[pos_strand]

# 对于负链基因
neg_strand <- annotation$strand == "-"
annotation$genebody_start[neg_strand] <- annotation$TES[neg_strand]
annotation$genebody_end[neg_strand] <- annotation$TSS[neg_strand] - 400
head(neg_strand)

# 5. 创建GRanges对象
promoter_gr <- GRanges(
  seqnames = annotation$chrom,
  ranges = IRanges(start = annotation$promoter_start, 
                   end = annotation$promoter_end),
  strand = annotation$strand,
  gene_id = annotation$name
)

genebody_gr <- GRanges(
  seqnames = annotation$chrom,
  ranges = IRanges(start = annotation$genebody_start, 
                   end = annotation$genebody_end),
  strand = annotation$strand,
  gene_id = annotation$name
)

# 6. 定义BAM文件列表
bam_files <- c(
  "35-Ser2_sorted_rmdup.bam",
  "35-Ser5_sorted_rmdup.bam",
  "36-Ser2_sorted_rmdup.bam",
  "36-Ser5_sorted_rmdup.bam",
  "B16-Ser2_sorted_rmdup.bam",
  "B16-Ser5_sorted_rmdup.bam",
  "P-Ser2_sorted_rmdup.bam",
  "P-Ser5_sorted_rmdup.bam"
)

# 7. 计算每个BAM文件的read counts
results <- list()

for(bam_file in bam_files) {
  sample_name <- gsub(".bam", "", basename(bam_file))
  
  # 读取BAM文件
  aln <- readGAlignments(bam_file)
  
  # 计算promoter counts
  promoter_counts <- countOverlaps(promoter_gr, aln)
  
  # 计算genebody counts
  genebody_counts <- countOverlaps(genebody_gr, aln)
  
  # 计算区域长度
  promoter_length <- width(promoter_gr)
  genebody_length <- width(genebody_gr)
  
  # 计算密度
  promoter_density <- promoter_counts / promoter_length
  genebody_density <- genebody_counts / genebody_length
  
  # 计算pausing index
  pausing_index <- promoter_density / (genebody_density + 0.001)  # 避免除零
  
  # 存储结果
  results[[sample_name]] <- data.frame(
    gene_id = annotation$name,
    promoter_count = promoter_counts,
    genebody_count = genebody_counts,
    promoter_density = promoter_density,
    genebody_density = genebody_density,
    pausing_index = pausing_index
  )
}

# 8. 保存结果
output_dir <- "pausing_index_R"
dir.create(output_dir, showWarnings = FALSE)

for(sample_name in names(results)) {
  write.table(results[[sample_name]],
              file = paste0(output_dir, "/", sample_name, "_pi.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# 9. 创建汇总矩阵
all_genes <- unique(annotation$name)
pi_matrix <- data.frame(gene_id = all_genes)

for(sample_name in names(results)) {
  pi_matrix[[sample_name]] <- results[[sample_name]]$pausing_index[
    match(all_genes, results[[sample_name]]$gene_id)
  ]
}

write.table(pi_matrix,
            file = paste0(output_dir, "/all_samples_pi_matrix.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Analysis completed!\n")
