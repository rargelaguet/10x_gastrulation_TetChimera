#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
io$outdir <- paste0(io$basedir,"/results/soupX")

#######################
## Load soupX output ##
#######################

io$dir <- "/Users/ricard/data/10x_gastrulation_TetChimera/original/td_tomato"
opts$batches <- c(
  "E75_TET_TKO",
  "E75_WT_Host",
  "E85_Rep1_TET_TKO",
  "E85_Rep1_WT_Host",
  "E85_Rep2_TET_TKO",
  "E85_Rep2_WT_Host",
  "SIGAE4_E105_3_TET123_Chimera_Host",
  "SIGAF4_E105_3_TET123_Chimera_TKO",
  "SIGAG4_E105_5_TET123_Chimera_Host",
  "SIGAH4_E105_5_TET123_Chimera_TKO"
)

soup.dt <- opts$batches %>% map(function(i) {
  file <- sprintf("%s/%s/soup/soup.tsv.gz",io$dir,i)
  fread(file) %>% .[,batch:=i]
}) %>% rbindlist

# Save
fwrite(soup.dt, sprintf("%s/soupX_estimates.txt.gz",io$outdir), sep="\t")


##########
## Plot ##
##########

to.plot <- soup.dt %>% 
  setorder(batch,-est) %>%
  split(.$batch) %>% map(~ head(.,n=25)) %>% rbindlist

p <- ggbarplot(to.plot, x="gene", y="est", fill="gray70") +
  coord_flip() +
  facet_wrap(~batch, scale="fixed", nrow=1) +
  labs(x="", "Soup relative (?) abundance") +
  theme(
    axis.text = element_text(size=rel(0.74))
  )

pdf(paste0(io$outdir,"/pdf/soup_abundance_E8.5.pdf"), width=8, height=7)
print(p)
dev.off()