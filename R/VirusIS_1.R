#' @import magrittr
#' @title FeatureCounts_txt_to_csv
#' @description process the file and reorder it.
#' @param file input txt from Featurecounts
#' @param out_path output filename
#'
#' @return reordered csv file
#' @export
FeatureCounts_txt_to_csv <- function(file = '/home/SharedData/VIS/first_test/04-featureCounts_by_exon.txt',
                                     out_path= '/home/SharedData/VIS/first_test/05-reordered-gene-counts.csv'
){

  dt <- fread(file)
  dt
  len <- ncol(dt) - 7

  names(dt)[7:(7+len)] <- paste('Sample', 1:(1+len), 'Counts', sep = '-')




  dt %>% head

  length <- dt %>% gather(key = 'Sample', value = 'reads', 7:ncol(dt)) %>% as.data.table()
  length <- length[reads > 0,]
  setorder(length, -reads)
  out <- length %>% spread(key = Sample,  value = reads)


  out

  write.csv(out, out_path)
  return(out)
}


#' @title import bed file
#'
#' @param file bed filename
#'
#' @return the contents of the bed file
#' @export
read_bed <- function(file){
  fread(file, header = F,  sep = '\t')
}


#' @title write dt to bed file
#'
#' @param dt data.table to be write
#' @param bed_file bed filename for output
#' @export
write_bed <- function(dt, bed_file){
  write_delim(x = dt, file = bed_file, col_names = F, delim = '\t')
}




#' @title seperate the gtf for reads ext
#'
#' @param input sorted.bed
#'
#' @return bed files were writing into files.
#' @export
seperate_gtf <- function(input = '/home/SharedData/Ref_Data/GFF_Files_for_Gene_Annotation/hsa_GRCh38_latest_genomic.gff-sorted.bed'){
  # sep gtf.bed to exon.bed and gene.bed
  # thus the bed can be used for other condition.
  # Such as intro mapping/ext, non-gene ext
  bed <- read_bed(input)
  # bed
  colnames(bed) <- c('chr', 'start', 'end', 'V4', 'V5', 'strand', 'source', 'type', 'V9', 'info')
  bed$info <- NULL
  bed <- bed[,.(chr, start, end, strand, source, type)]
  gene.bed <- bed[type == 'gene',]
  exon.bed <- bed[type == 'exon',]

  gene.file <- input %>% sub('.bed$', '-gene.bed',.)
  exon.file <- input %>% sub('.bed$', '-exon.bed',.)

  write_bed(dt = gene.bed, bed_file = gene.file)
  write_bed(dt = exon.bed, bed_file = exon.file)
}


#' @title rename chr name to 'chr1' mod.
#'
#' @param chrs chrs
#'
#' @return chrs with new names formatted as 'chr1' to 'chrY'
#' @export
rename_chr <- function(chrs = c('NC_000012.12', 'NC_000001.11')){
  # Note
  # 'NC_000023.11' -> chrX
  # 'NC_000024.11' -> chrY
  # 'NC_012920.1'  -> chrM
  # 'NC_00000X.10' %>% sub('NC_[0]+([1-9XYM].*)[.][0-9]+$', 'chr\\1', .)
  # NT: unlocalized contigs.

  new_chrs <- chrs %>% sub('NC_[0]+([1-9].*)[.][0-9]+$', 'chr\\1', .)
  new_chrs[new_chrs == 'chr23'] <- 'chrX'
  new_chrs[new_chrs == 'chr24'] <- 'chrY'
  new_chrs[new_chrs == 'chr12920'] <- 'chrM'
  return(new_chrs)
}


#' @title Get bedgraph file for IGV from bed file
#'
#' @param bed_file input bed file
#' @param set_end_equal_start T or F; when set as T(default), the start will be set as end too.
#'
#' @return bedgraph for IGV: with 4 columns: ① chr;  ② start; ③ end; ④
#' @export
bedgraph_for_IGV <- function(bed_file, set_end_equal_start = T){
  dt <- read_bed(bed_file)
  dt$V1 <- dt$V1 %>% rename_chr
  dt <- dt[grepl('chr', V1),] ## only subset those can matched to the 'chr'
  # dt$V1 %>% table
  if(set_end_equal_start){
    dt$V3 <- dt$V2 # 对于VIS位点，插入的起始位置最重要。不需要保留插入终止位置。
  }
  dt <- dt[,.(V1, V2, V3)]
  setorder(dt, V1, V2)
  out_bed <- bed_file %>% sub('.bed$', '-renamed+sorted.bed', .)

  write_bed(dt, bed_file = out_bed)

  bedgraph <- dt[, .(rank = 1:.N, sum = .N), by = list(V1, V2, V3)]
  bedgraph <- bedgraph[rank == 1, ]
  bedgraph$rank <- NULL
  # bedgraph[grep('chr', V1),]$V1 %>% table %>% names
  # bedgraph[grep('chr', V1),]$V1 %>% table %>% names
  # bedgraph
  out_graph <-  bed_file %>% sub('.bed$', '-renamed+sorted.bedgraph', .)
  write_bed(bedgraph, bed_file = out_graph)
}


#' Title get tss dist plot from tss-dist.bed file
#'
#' @param tss_dist_bed bed contains tss_dist
#' @param chr_col num for chr column.
#' @param dist_cols list out the position from two files 列出连个文件中的位置所在的列，两个文件各列出1个碱基就够了。
#' The distance will be calculated again.
#'
#'
#' @return pic by ggplot
#' @export
tss_dist_plot <- function(tss_dist_bed = '/home/SharedData/VIS/2022-07-second-sample/demo/Res-11-gene-mapped-reads-dist_to_tss.bed', chr_col = 1, dist_cols = c(7,11)){
  # 可用于测试的参考输入'/home/SharedData/VIS/2022-07-second-sample/demo/Res-09-dist_to_tss.bed'
  # '/home/SharedData/VIS/2022-07-second-sample/demo/Res-11-gene-mapped-reads-dist_to_tss.bed'


  bed <- read_bed(file = tss_dist_bed)
  #文件1 中毒额位置 - 文件2中的位置才是实际的距离。
  head(bed)
  tail(bed)
  bed %>% nrow

  index <- grepl('chr[0-9XY]+$', bed$V1) ## 提取出能够匹配的，且匹配到染色体上面的结果。不带线粒体，如果带线粒体就加个M
  bed <- bed[index,]
  bed %>% unique() %>% nrow

  ## 提取出两个文件中的距离, 使用with = F，直接传入数据形式的变量。
  bed <- bed[,c(chr_col, dist_cols), with = F] %>% unique
  bed %>% nrow
  # bed
  colnames(bed) <- c('chr', 'position.in.A', 'position.in.B')

  bed$dist <- bed$position.in.A - bed$position.in.B

  dt <- bed
  dt

  dt$dist <- dt$dist %>% as.numeric()
  dt$tss_dist <- cut(dt$dist/1000,
                     breaks = c(-Inf, -500,          -50,          -5,         0,       5,         50,          500,    Inf),
                     labels = c('<-500', '-500 to -50', '-50 to -5', '-5 to 0',  '0 to 5', '5 to 50', '50 to 500', '> 500'), ordered_result = T)


  pic_dt <- dt[,.(count=.N), by = tss_dist]
  pic_dt$percent <- 100 * pic_dt$count/sum(pic_dt$count)
  pic <- ggplot(pic_dt, aes(x = tss_dist, y = percent))+
    ylab('Percent')+ xlab("base Distance (K)")+
    geom_bar(fill = '#3db874', color = 'grey20', stat = 'identity') + theme_bw()+
    theme(
      text = element_text(size = 18, face = 'bold'),
      axis.text =   element_text(size = 18, color = 'black'),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
      axis.title  = element_text(size = 20, face = 'bold'),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      plot.subtitle = element_text(color = 'black', size = 12, face = 'plain', hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "cm") ### leave margins
    )
  pic
  return(pic)




  if(F){
    # test
    dt[,.(chr, position.in.B)] %>% unique
    dt[,.(chr, position.in.A)] %>% unique
    min_max <- dt[,.(min=min(`position.in.A`), max = max(position.in.A), counts = .N),by=chr] %>% copy
    min_max$range <- min_max$max - min_max$min
    min_max$range %>% sum / 59621
    min_max

    min_max$range %>% sum / 223865

  }
}




#' Title sort the bed file for tss dist mapping.
#'
#' @param bed_path path to bed file contains dist-to-tss
#'
#' @return dist-files
#' @export
sort_bed_for_tss <- function(bed_path = '/home/SharedData/VIS/first_test/find-VIS-08-b-bed-from-sam.bed'){
  # bed文件是tab分隔的。
  bed <- fread(bed_path, header = F, sep = '\t')

  names(bed) <- c('chr', 'start', 'end', 'info', 'MAPQ', 'strand')
  # MAPQ stands for mapping quality. a zero MAPQ indicates that the seq may be multi mapped.

  # bed %>% head
  bed$chr %>% table
  bed %>% is.data.table()
  dt <- bed[grepl('^NC_', chr),]
  dt[,chr_new := sub('^NC_0+([1-9]+[0-9]?)[.].*$', '\\1', chr) %>% as.numeric]
  dt$chr_new %>% table
  # dt %>% str
  # dt
  dt$chr_new <- dt$chr_new %>% factor(., levels = 1:24, labels = c(1:22, 'X', 'Y'))
  dt$chr_new <- paste('chr', dt$chr_new, sep = '')
  dt$chr_new %>% table
  dt$chr <- dt$chr_new
  dt$chr_new <- NULL
  dt$end <- dt$start
  setorder(dt, chr, start) # reorder the bed file. by chr(in character) and position.
  # dt %>% str
  # dt$chr %>% unique
  out_path <- bed_path %>% sub('.bed$', '-sorted.bed', .)
  write_delim(dt,delim = '\t', col_names = F, file = out_path)
  message("the bed file was sorted and chr was renamed.")
}


#' @title get top 10 mapped gene.
#' @description get 10 mapped .
#' @param path_to_txt input a path to the
#' @return pic for top 10
top10_exon_plot <- function(path_to_txt = '/home/SharedData/VIS/2022-07-second-sample/Rm_contam-06-gene-featurecounts-by-exon-cleaned.txt'){
  cleaned <- fread(path_to_txt)
  sum <- cleaned$`Sample-1-Counts` %>% sum
  cleaned$percent <- 100 * cleaned$`Sample-1-Counts`/sum
  # cleaned

  pic <- cleaned[1:10,] %>% ggplot(data = ., aes(x = 'Top 10',y = percent, fill = Geneid)) +
    geom_bar(stat = "identity", position = "stack", width = .88, color = 'black') +
    # scale_fill_manual(values = c('grey45', '#A0EEE1')) +
    theme_bw()  +
    xlab(NULL)+
    ylab('Percentage %')+
    theme(
      text = element_text(size = 18, face = 'bold'),
      axis.text =   element_text(size = 18, color = 'black'),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
      axis.title  = element_text(size = 20, face = 'bold'),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      plot.subtitle = element_text(color = 'black', size = 12, face = 'plain', hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "cm") ### leave margins
    )
  return(pic)
}



#' Title get the percent plot for the top 10 sites.
#'
#' @param cleaned_bedgraph cleaned bedgraph
#'
#' @return plot with the top10 VIS
#' @export
top10_site_plot <- function(cleaned_bedgraph = '/home/SharedData/VIS/2022-07-second-sample/Rm_contam-10-exon-renamed+sorted.bedgraph'){

  # top10_bed_site() # 获取前十位点对应的比例
  for_pic <- fread(cleaned_bedgraph)

  colnames(for_pic) <- c('chr', 'position', 'position-another', 'count')
  for_pic$matcher <- paste(for_pic$chr, for_pic$position)
  setorder(for_pic, -count, matcher)
  # for_pic[,rank:=1:.N, by = matcher] # 不应该有重复的
  # for_pic[rank> 1,] #
  for_pic[, percent:=100 * count/(sum(count))]
  for_pic

  pic <- for_pic[1:10,] %>% ggplot(data = ., aes(x = 'Top 10',y = percent, fill = matcher)) +
    geom_bar(stat = "identity", position = "stack", width = .88, color = 'black') +
    # scale_fill_manual(values = c('grey45', '#A0EEE1')) +
    theme_bw()  +
    xlab(NULL)+
    ylab('Percentage %')+
    labs(fill = 'Position')+
    theme(
      text = element_text(size = 18, face = 'bold'),
      axis.text =   element_text(size = 18, color = 'black'),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
      axis.title  = element_text(size = 20, face = 'bold'),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      plot.subtitle = element_text(color = 'black', size = 12, face = 'plain', hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "cm") ### leave margins
    )
  # pic
  return(pic)

}







#' Title get logo plot
#'
#' @param cleaned_mapped_fastq  fq with the mapped reads.
#'
#' @return logo-plot
#' @export
get_logo_plot <- function(cleaned_mapped_fastq = '/home/SharedData/VIS/first_test/get-break-sub-03-trim_half_adapter.fastq'){
  fq <- fread(cleaned_mapped_fastq, sep = NULL, header = F)
  # fq %>% head

  len <- nrow(fq)/4
  index <- (1:len) * 4 - 2
  seqs <- fq[index]
  # seqs %>% head

  { # 试图提取reads的id，但是这里好像看不出来 R1 还是R2.
    names <- (1:len) * 4 - 3
    reads_id <- fq[names]

  }

  seqs$V1 <- seqs$V1 %>% substr(., 1, 20)
  index <- seqs$V1 %>% nchar == 20
  (index == F) %>% sum
  seqs <- seqs[index,]
  # nrow(seqs)
  {
    # logo_plot <- '/home/SharedData/VIS/first_test/get-break-sub-07-Logo-plot.pdf'
    pic <- ggseqlogo(seqs$V1 # , seq_type="aa"
    ) +
      theme_logo() +
      theme(text = element_text(size = 12),
            axis.text =  element_text(size = 14),
            axis.title = element_text(size = 18),
            title = element_text(size = 20),
            strip.text = element_text(size = 20))


  }
  pic
  return(pic)
}


