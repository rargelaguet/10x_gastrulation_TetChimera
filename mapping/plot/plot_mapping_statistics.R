library(data.table)
library(purrr)
library(ggplot2)

sample_metadata <- fread("/Users/ricard/data/scnmt_gastrulation_TetKO/rna/mapping_10x/sample_metadata_mapping_mnn.txt") %>%
  .[!lineage10x_2%in%c("ExE ectoderm","ExE endoderm")]
  .[!is.na(lineage10x_2) ]

colors_lineage10x_2 <- c(
  "Epiblast" = "grey50",
  "Ectoderm" = "steelblue",
  # "ExE ectoderm" = "steelblue",
  "Primitive streak" = "sandybrown",
  "Mesoderm" = "#CD3278",
  # "Primitive endoderm" = "#2E8B57",
  # "ExE endoderm" = "#2E8B57",
  "Endoderm" = "#43CD80",
  "Blood" = "black"
)

colors_lineage10x <- c("Epiblast" = "#635547",
                      "Primitive Streak" = "#DABE99",
                      "Caudal epiblast" = "#9e6762",
                      
                      "PGC" = "#FACB12",
                      
                      "Anterior Primitive Streak" = "#c19f70",
                      "Notochord" = "#0F4A9C",
                      "Def. endoderm" = "#F397C0",
                      "Gut" = "#EF5A9D",
                      
                      "Nascent mesoderm" = "#C594BF",
                      "Mixed mesoderm" = "#DFCDE4",
                      "Intermediate mesoderm" = "#139992",
                      "Caudal Mesoderm" = "#3F84AA",
                      "Paraxial mesoderm" = "#8DB5CE",
                      "Somitic mesoderm" = "#005579",
                      "Pharyngeal mesoderm" = "#C9EBFB",
                      "Cardiomyocytes" = "#B51D8D",
                      "Allantois" = "#532C8A",
                      "ExE mesoderm" = "#8870ad",
                      "Mesenchyme" = "#cc7818",
                      
                      "Haematoendothelial progenitors" = "#FBBE92",
                      "Endothelium" = "#ff891c",
                      "Blood progenitors 1" = "#f9decf",
                      "Blood progenitors 2" = "#c9a997",
                      "Erythroid1" = "#C72228",
                      "Erythroid2" = "#f79083",
                      "Erythroid3" = "#EF4E22",
                      
                      "NMP" = "#8EC792",
                      
                      "Rostral neurectoderm" = "#65A83E",
                      "Caudal neurectoderm" = "#354E23",
                      "Neural crest" = "#C3C388",
                      "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                      "Spinal cord" = "#CDE088",
                      
                      "Surface ectoderm" = "#f7f79e",
                      
                      "Visceral endoderm" = "#F6BFCB",
                      "ExE endoderm" = "#7F6874",
                      "ExE ectoderm" = "#989898",
                      "Parietal endoderm" = "#1A1A1A"
                      
)

to.plot <- sample_metadata[,.N, by=c("stage","lineage10x_2","genotype")] %>% 
  # .[, lineage10x_2:=stringr::str_replace_all( lineage10x_2,"_"," ")] %>%
  .[, lineage10x_2:=factor( lineage10x_2,levels=names(colors_lineage10x_2))]
  
stopifnot(all(to.plot$lineage10x_2 %in% names(colors_lineage10x_2)))

p <- ggplot(to.plot, aes(x= lineage10x_2, y=N)) +
  geom_bar(aes(fill= lineage10x_2), stat="identity", color="black") +
  scale_fill_manual(values=colors_lineage10x_2) +
  facet_wrap(~genotype, nrow=1, scales="fixed") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(1.3)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.text.x = element_text(size=rel(1.1), color="black")
  )

# pdf("/Users/ricard/data/gastrulation/mapping_10x/mapping_stats/mapping_lineage10x_2.pdf", width=11, height=5)
print(p)
# dev.off()






to.plot <- sample_metadata[,.N, by=c("stage","lineage10x","genotype","embryo")] %>% 
  .[complete.cases(.)] %>%
  .[, lineage10x:=factor( lineage10x,levels=names(colors_lineage10x))]

stopifnot(all(to.plot$lineage10x %in% names(colors_lineage10x)))

p <- ggplot(to.plot, aes(x= lineage10x, y=N)) +
  geom_bar(aes(fill=lineage10x), stat="identity", color="black") +
  scale_fill_manual(values=colors_lineage10x) +
  facet_wrap(~genotype+embryo, nrow=1, scales="fixed") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(1.3)),
    axis.title.x = element_text(color="black", size=rel(1.1)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.text.x = element_text(size=rel(1.1), color="black")
  )

# pdf("/Users/ricard/data/gastrulation/mapping_10x/mapping_stats/mapping_lineage10x_2.pdf", width=11, height=5)
print(p)
# dev.off()
