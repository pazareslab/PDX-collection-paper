library(rstatix)
library(tidyverse)
library(progeny) 
library(GenVisR)
library(DESeq2)
library(pheatmap)
library(reshape2)
library(OmnipathR)
library(decoupleR)
library(pheatmap)

# Set the path where de .rds files are located:

path_to_rds <- "~/Desktop/figure_generation/"

# Load clinical data:

Clinical <- readRDS(paste0(path_to_rds,"clinical_data.rds"))
ADCnames <- Clinical %>% filter(Histology=="ADC") %>% pull(Sample)
SCCnames <- Clinical %>% filter(Histology=="SCC") %>% pull(Sample)

# Load CNV data and create PDF plots:

CNVs <- readRDS(paste0(path_to_rds,"CNVs.rds"))

ADC_NOX <- CNVs %>%
  filter(Sample %in% ADCnames & weight>1000 & cn!=2) %>%
  select(sample=Sample, chromosome, start, end, segmean=cn)

SCC_NOX <- CNVs %>%
  filter(Sample %in% SCCnames & weight>1000 & cn!=2) %>%
  select(sample=Sample, chromosome, start, end, segmean=cn)

pdf(paste0(path_to_rds,"ADC_CNV_PDX.pdf"),width = 20,height = 10)
cnSpec(ADC_NOX, genome="hg38",CNscale = "absolute")
dev.off()

pdf(paste0(path_to_rds,"SCC_CNV_PDX.pdf"),width = 20,height = 10)
cnSpec(SCC_NOX, genome="hg38",CNscale = "absolute")
dev.off()


# Get info about the TMB and the top mutated genes pero histology, useful for the waterfall plot and the expression heatmap annotations:

Aberrations <- readRDS(paste0(path_to_rds,"Aberrations.rds"))
TMB_PDL1_ADC <- Clinical %>% filter(Histology=="ADC") %>% select(Sample,"TMB High" = `TMB Consensus Status`,"PDL1 Positive" ="PDL1 Status")  %>% column_to_rownames("Sample")
top_smallvars_ADC <- Aberrations %>% filter(Histology == "ADC" &VariantCategory=="Small Variant" ) %>% select(Sample,Gene) %>% group_by(Gene) %>% 
  mutate(occ=n()) %>% ungroup %>% filter(occ>=3) %>% unique %>% mutate(matrix=1) %>% select(-occ) %>% spread(Sample,matrix,fill=0) %>% column_to_rownames("Gene") %>% t %>% as.data.frame()
colnames(top_smallvars_ADC) <- top_smallvars_ADC %>% colnames %>% paste0(.," Mut")

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

annotation_progeny_ADC <- cbind(top_smallvars_ADC,TMB_PDL1_ADC)

TMB_PDL1_SCC <- Clinical %>% filter(Histology=="SCC") %>% select(Sample,"TMB High" = `TMB Consensus Status`,"PDL1 Positive" ="PDL1 Status") %>% column_to_rownames("Sample")

top_smallvars_SCC <- Aberrations %>% filter(Histology == "SCC" &VariantCategory=="Small Variant" ) %>% select(Sample,Gene) %>% group_by(Gene) %>% 
  mutate(occ=n()) %>% ungroup %>% filter(occ>=4) %>% unique %>% mutate(matrix=1) %>% select(-occ) %>% spread(Sample,matrix,fill=0) %>% column_to_rownames("Gene") %>% t %>% as.data.frame()

colnames(top_smallvars_SCC) <- top_smallvars_SCC %>% colnames %>% paste0(.," Mut")
top_CNVs_SCC <- Aberrations %>% filter(Histology == "SCC" &VariantCategory=="CNV" ) %>% select(Sample,Gene) %>% group_by(Gene) %>% 
  mutate(occ=n()) %>% ungroup %>% filter(occ>=3) %>% unique %>% mutate(matrix=1) %>% select(-occ) %>% spread(Sample,matrix,fill=0) %>% column_to_rownames("Gene") %>% t %>% as.data.frame 
add_to_cnvs <- outersect(SCCnames,top_CNVs_SCC %>% rownames)
if (length(add_to_cnvs)>0){
  for (i in 1:length(add_to_cnvs)){
    top_CNVs_SCC[add_to_cnvs[i],] <- 0
  }
  
}

colnames(top_CNVs_SCC) <- top_CNVs_SCC %>% colnames %>% paste0(.," Amp")
annotation_progeny_SCC <- cbind(top_smallvars_SCC,top_CNVs_SCC,TMB_PDL1_SCC)

common_annotations <- intersect(annotation_progeny_ADC %>% colnames(), annotation_progeny_SCC %>% colnames())
unique_annotations <- outersect(annotation_progeny_ADC %>% colnames(), annotation_progeny_SCC %>% colnames())

annotation_colors_ADC <-  list("TP53 Mut" = c("1"="black", "0" = "white"),
                               "TMB High" = c("1"="darkblue", "0" = "white"),
                               "PDL1 Positive" = c("1" = "darkgreen", "0" = "white"),
                               "EGFR Mut" = c("1" = "#B15928", "0" = "white"),                             
                               "KRAS Mut" = c("1" = "#FFFF33", "0" = "white"),                             
                               "MGA Mut" = c("1" = "#FF7F00", "0" = "white"),                             
                               "PIK3CA Mut" = c("1" = "#984EA3", "0" = "white"),                             
                               "RB1 Mut" = c("1" = "#DECBE4", "0" = "white"),                             
                               "STK11 Mut" = c("1" = "#7FC97F", "0" = "white"),                             
                               "FGFR1 Del" = c("1" = "#36454F", "0" = "white")                          
)

rownames_adc <- annotation_progeny_ADC %>% rownames
annotation_progeny_ADC <- annotation_progeny_ADC %>% apply(.,2,as.character) %>% as.data.frame
rownames(annotation_progeny_ADC) <- rownames_adc

annotation_colors_SCC <-  list("TP53 Mut" = c("1"="black", "0" = "white"),
                               "MYC Amp" = c("1"="darkred", "0" = "white"),
                               "TMB High" = c("1"="darkblue", "0" = "white"),
                               "PDL1 Positive" = c("1" = "darkgreen", "0" = "white"),
                               "PIK3CA Mut" = c("1" = "#984EA3", "0" = "white"),                             
                               "PIK3CA Amp" = c("1" = "#7FC97F", "0" = "white"),
                               "MET Amp" = c("1" = "#DECBE4", "0" = "white"),
                               "NOTCH1 Mut" = c("1" = "#E5D8BD", "0" = "white"),
                               "KMT2D Mut" = c("1" = "#A65628", "0" = "white"),
                               "KEAP1 Mut" = c("1" = "#FFFF33", "0" = "white"),
                               "GNAS Mut" = c("1" = "#B15928", "0" = "white"),
                               "FAT1 Mut" = c("1" = "#984EA3", "0" = "white"),
                               "ESR1 Mut" = c("1" = "#FFD92F", "0" = "white"),
                               "CDKN2A Mut" = c("1" = "#FF7F00", "0" = "white"),
                               "ARID1A Mut" = c("1" = "#F781BF", "0" = "white")
)
rownames_scc <- annotation_progeny_SCC %>% rownames
annotation_progeny_SCC <- annotation_progeny_SCC %>% apply(.,2,as.character) %>% as.data.frame
rownames(annotation_progeny_SCC) <- rownames_scc

# Load RNA-seq gene quantification data:
Exp <- readRDS(paste0(path_to_rds,"Expression.rds"))

# Filtering of low expressed and low variable genes:
vst_ADC <- Exp %>% select(Gene,ADCnames) %>% column_to_rownames("Gene") %>%
  filter_at(.,vars(contains("TP")),any_vars(. !=0)) %>%
  filter_at(.,vars(contains("TP")),any_vars(. > 10)) %>%
  mutate(media=dplyr::select(.,contains("TP")) %>% rowMeans) %>%
  filter(media>10) %>% dplyr::select(-media) %>%  
  as.matrix %>% varianceStabilizingTransformation()

vst_SCC <- Exp  %>% select(Gene,SCCnames) %>% column_to_rownames("Gene") %>%
  filter_at(.,vars(contains("TP")),any_vars(. !=0)) %>%
  filter_at(.,vars(contains("TP")),any_vars(. > 10)) %>%
  mutate(media=dplyr::select(.,contains("TP")) %>% rowMeans) %>%
  filter(media>10) %>% dplyr::select(-media) %>%  as.matrix %>%
  varianceStabilizingTransformation()


# Load progeny:
net <- get_progeny(organism = 'human', top = 100)

# Run wmean:
sample_acts_ADC <- run_wmean(mat=vst_ADC %>% as.matrix, net=net, .source='source', .target='target',
                             .mor='weight', times = 100, minsize = 5) %>% filter(statistic=="norm_wmean") %>% adjust_pvalue(p.col = "p_value",method = "BH") %>% add_significance(p.col = "p_value")

sample_acts_SCC <- run_wmean(mat=vst_SCC %>% as.matrix, net=net, .source='source', .target='target',
                             .mor='weight', times = 100, minsize = 5) %>% filter(statistic=="norm_wmean") %>% adjust_pvalue(p.col = "p_value",method = "BH") %>% add_significance(p.col = "p_value")

sample_acts_ADC %>% dplyr::rename(pathway=source,sample=condition) %>% write_tsv(paste0(path_to_rds,"results_progeny_ADC_8ago22.tsv"))
sample_acts_SCC %>% dplyr::rename(pathway=source,sample=condition) %>% write_tsv(paste0(path_to_rds,"results_progeny_SCC_8ago22.tsv"))



# Transform to wide matrix:
sample_acts_mat_ADC <- sample_acts_ADC %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

sample_acts_mat_SCC <- sample_acts_SCC %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Scale per sample:
sample_acts_mat_ADC_scaled <- scale(sample_acts_mat_ADC)
sample_acts_mat_SCC_scaled <- scale(sample_acts_mat_SCC)

# Choose color palette:
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot progeny heatmaps:
pdf(paste0(path_to_rds,"ADC_pheatmap_progeny.pdf"), width = 10, height = 8)
results_pheamtap_ADC <- pheatmap(sample_acts_mat_ADC_scaled %>% t, border_color = NA, color=my_color, breaks = my_breaks, annotation_col = annotation_progeny_ADC,
                                 annotation_colors = annotation_colors_ADC, annotation_legend = F) 
dev.off()

pdf(paste0(path_to_rds,"SCC_pheatmap_progeny.pdf"), width = 10, height = 8)
results_pheamtap_SCC <- pheatmap(sample_acts_mat_SCC_scaled %>% t, border_color = NA, color=my_color, breaks = my_breaks, annotation_col = annotation_progeny_SCC,
                                 annotation_colors = annotation_colors_SCC, annotation_legend = F) 
dev.off()


# Get Waterfall mutation plots per histology:

Aberrations$VariantEffect[Aberrations$VariantEffect=="frameshift_insertion"] <- "Frameshift"
Aberrations$VariantEffect[Aberrations$VariantEffect=="frameshift_deletion"] <- "Frameshift"
Aberrations$VariantEffect[Aberrations$VariantEffect=="frameshift_substitution"] <- "Frameshift"
Aberrations$VariantEffect[Aberrations$VariantEffect=="nonframeshift_insertion"] <- "Nonframeshift"
Aberrations$VariantEffect[Aberrations$VariantEffect=="nonframeshift_deletion"] <- "Nonframeshift"
Aberrations$VariantEffect[Aberrations$VariantEffect=="nonframeshift_substitution"] <- "Nonframeshift"
Aberrations$VariantEffect[Aberrations$VariantEffect=="missense"] <- "Missense"
Aberrations$VariantEffect[Aberrations$VariantEffect=="splicing"] <- "Splicing"
Aberrations$VariantEffect[Aberrations$VariantEffect=="stopgain"] <- "Stopgain"
Aberrations$VariantEffect[Aberrations$VariantEffect=="startloss"] <- "Startloss" 
Aberrations$VariantEffect[Aberrations$VariantEffect=="synonymous"] <- "Splicing" #These are the "TP53:NM_000546:exon4:c.G375T:p.T125T"

Aberrations$ExtraInfo <- Aberrations$ExtraInfo %>% gsub("p.","",.)

oncoADC <-   Aberrations %>% filter(Histology=="ADC") %>% 
  select(sample=Sample, gene=Gene, variant_type=ExtraInfo, variant_class=VariantEffect) %>% as.data.frame 

most_deleterious <- c("Double","Fusion","Startloss","Stopgain","Frameshift","Splicing","Del","Missense","Amp","Nonframeshift")
colors <- c("black","slategrey","#762A83","magenta","darkorange","dodgerblue","blue","darkgreen","darkred","yellowgreen")
clinicalADC <- Clinical %>% filter(Histology=="ADC") %>% select(sample=Sample, "TMB" = "TMB Consensus Status") #%>% gather(key="variable", value= value )
clinicalADC <- melt(clinicalADC, id.vars = c("sample"))
clinicalADC <- clinicalADC %>% mutate(value = replace(value, value==0, "Low (<10)")) %>% 
  mutate(value = replace(value, value==1, "High (>=10)"))
lung_cancer_genes <- readRDS(paste0(path_to_rds,"lung_cancer_genes.rds"))

pdf(paste0(path_to_rds,"waterfall_ADC_PDX.pdf"), width = 30, height = 15)
waterfall(oncoADC,
          plotGenes = lung_cancer_genes,
          mainDropMut = TRUE,
          rmvSilent = TRUE,
          mainLabelCol = "variant_type",
          clinData = clinicalADC,
          fileType = "Custom",
          variant_class_order = most_deleterious,
          clinVarCol = c("High (>=10)" =  "darkblue",
                         "Low (<10)" = "lightblue"),          
          mainPalette = colors,
          mainXlabel = TRUE,
          section_heights = c(0,1,0.05),
          main_geneLabSize = 15,
          plotMutBurden = FALSE,
          mainLabelSize = 5,
          mainGrid = F
)
dev.off()

oncoSCC <-  Aberrations %>% filter(Sample %in% SCCnames) %>% 
  select(sample=Sample, gene=Gene, variant_type=ExtraInfo, variant_class=VariantEffect) %>% as.data.frame

clinicalSCC <- Clinical %>%filter(Sample %in% SCCnames)%>% select(sample=Sample, "TMB" = "TMB Consensus Status", Histology) %>% gather(key="variable", value= value,-sample )
clinicalSCC<- clinicalSCC %>% mutate(value = replace(value, value==0, "Low (<10)")) %>% 
  mutate(value = replace(value, value==1, "High (>=10)"))
pdf(paste0(path_to_rds,"waterfall_SCC_PDX.pdf"), width = 30, height = 15)
waterfall(oncoSCC,
          plotGenes = lung_cancer_genes,
          mainDropMut = TRUE,
          rmvSilent = TRUE,
          mainLabelCol = "variant_type",
          clinData = clinicalSCC,
          fileType = "Custom",
          variant_class_order = most_deleterious,
          clinVarCol = c("High (>=10)" =  "darkblue",
                         "Low (<10)" = "lightblue"),          
          mainPalette = colors,
          mainXlabel = TRUE,
          section_heights = c(0,1,0.05),
          main_geneLabSize = 15,
          plotMutBurden = FALSE,
          mainLabelSize = 4.5, #4,5
          mainGrid = F
)

dev.off()


sessionInfo()
