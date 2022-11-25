#Set disease id and date of run, change when running different diseases and updated pipelines
#In this example we are using schizophrenia
Disease <- 'Schizophrenia'
Date <- '091222'

# Load necessary libraries
library(tidyverse)
library(topGO)
library(reshape2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(STRINGdb)
string_db <- STRINGdb$new(version="11", species=9606, score_threshold=0, input_directory="" )
library('jsonlite')

#Annotation resources, most of these files can be generated using the 'UpdatingPipelineResources.R' code
DiseaseGenes <- readRDS(paste0('./PipelineData/',Disease,'Genes.rds'))
GTExbraindata <- readRDS('./PipelineData/GTExbraindata.rds')
drugtargets <- readRDS('./PipelineData/drugtargets.rds')
PharmGKB <- readRDS('./PipelineData/PharmGKB.rds')
HumanGenes_Mouseortholog <- readRDS(paste0('./PipelineData/HumanGeneMouseOrthologs_', Date, '.rds'))
hsproteincodinggenes <-  readRDS('./PipelineData/hsproteincodinggenes.rds')
mouseproteincodinggenes <- readRDS('./PipelineData/mmproteincodinggenes.rds')
ACMGGenes <- readRDS('./PipelineData/ACMGgenelist.rds')
#Read in GOtoMP mappings downloaded from https://www.ebi.ac.uk/spot/oxo/search# and mapped to top level IMPC terms
GOtoMP_mappings <- read.csv('./PipelineData/GOtoMP_mappings_MPToplevel.csv', header=TRUE)

# Initiate lists with results
DiseaseGene.Results <- NULL

# 1. Drug Targets and pharmacogenomics:
DiseaseGenes.DrugTargets_proteins <-
  left_join(DiseaseGenes, drugtargets[c('UniProtID', 'idgTDL')], by = 'UniProtID')# results by gene and keep protein

DiseaseGenes.DrugTargets.summary <- DiseaseGenes.DrugTargets_proteins %>%
  distinct(UniProtID, .keep_all = TRUE) %>%
  group_by(idgTDL) %>%
  summarise(n = n()) # Summarized results

PharmGKB <- dplyr::select(PharmGKB, c(1:14))
DiseaseGenes.PharmGKB <-
  left_join(DiseaseGenes, PharmGKB, by = c('Gene.Name'='Gene.Symbol')) # results by Gene.Name since these are what is stored in PharmGKB
DiseaseGenes.PharmGKB$Variant <- as.factor(DiseaseGenes.PharmGKB$Variant)
DiseaseGenes.PharmGKB$Chemical <- as.factor(DiseaseGenes.PharmGKB$Chemical)
DiseaseGenes.PharmGKB.summary <- DiseaseGenes.PharmGKB %>%
  filter(Significance=='yes') %>%
  group_by(EntrezID) %>%
  mutate(nVars = n_distinct(Variant), nDrugs = n_distinct(Chemical), nDrugVarCombo=n_distinct(Variant, Chemical)) %>%
  ungroup() %>%
  distinct(EntrezID, Gene.Name, nVars, nDrugs, nDrugVarCombo) # Summarized results

# 2. Gene Expression:
#Create dataframe with Disease candidate genes expressed in human brain
Diseasebraingenes <- left_join(DiseaseGenes, GTExbraindata, by=c('EnsID'='EnsID.GTEx'))
Diseasebraingenes <- dplyr::select(Diseasebraingenes, #DisGeNETScore, 
                                   Gene.Name, EntrezID, EnsID, UniProtID, Brodman1909.area24, Brodmann1909.area9, 
                                   Amygdala, Caudate.Nucleus, Cerebellar.Hemisphere, Cerebellum, Cerebral.Cortex, Hippocampus.Proper,
                                   Hypothalamus, Nucleus.Accumbens, Pituitary.Gland, Putamen, Substantia.Nigra)

#Prepare GTEx data for summaries
meltGTEx <- melt(GTExbraindata, by = 'EnsID.GTEx')

DiseaseGenes.Brain <-
  left_join(DiseaseGenes, meltGTEx, by = c('EnsID'= 'EnsID.GTEx')) # results by gene

DiseaseGenes.Brain.summary <- DiseaseGenes.Brain %>%
  filter(!is.na(value)) %>%
  distinct() %>%
  group_by(EnsID) %>%
  summarise(
    nRegions = n(),
    meanexp = mean(value),
    sdexp = sd(value)
  ) # Summarized results by gene

DiseaseGenes.Brain.totals <- DiseaseGenes.Brain %>%
  filter(!is.na(value)) %>%
  distinct() %>%
  group_by(variable) %>%
  summarise(nGenes = n())

#3. Gene Set Overrepresentation
all_genes <-
  unique(as.character(hsproteincodinggenes$EnsID))
geneUniverse <- rep(0, length(all_genes))
names(geneUniverse) <- all_genes

DiseaseGenes_topGO <- DiseaseGenes %>%
  filter(!is.na(EnsID))
genesOfInterest <- DiseaseGenes_topGO$EnsID
geneUniverse[names(geneUniverse) %in% genesOfInterest] <- 1

GOBPdata_DiseaseGenes <- new(
  "topGOdata",
  description = "testdata",
  ontology = "BP",
  allGenes = geneUniverse,
  geneSel = function(p)
    p == 1,
  annot = annFUN.org,
  ID = "ensembl",
  mapping = "org.Hs.eg.db"
)

# Run Overrepresentation Test
BPGSA_DiseaseGenes <-
  runTest(GOBPdata_DiseaseGenes,
          algorithm = "parentchild",
          statistic = "fisher")

#Note p-value threshold is based on bonferroni of 0.05/(GO DAG topologies tests in TopGO)
testednodes <- n_distinct(GOBPdata_DiseaseGenes@graph@nodes)

BPGSA_DiseaseGenes.results_DiseaseGenes <-
  GenTable(GOBPdata_DiseaseGenes,
           classicFisher = BPGSA_DiseaseGenes,
           topNodes = testednodes)

#Pull significant results bonferroni correcting for number of nodes
#Add full p-values and calculate fold enrichment
Topscores <- as.data.frame(BPGSA_DiseaseGenes@score)
write.csv(Topscores, './Results/Topscores.csv')
Fullpvalues <- read.csv('./Results/Topscores.csv')
BPs_DiseaseGenes <- left_join(BPGSA_DiseaseGenes.results_DiseaseGenes, Fullpvalues, by=c('GO.ID'='X'))
BPs_DiseaseGenes <- BPs_DiseaseGenes %>%
  mutate(FoldEnrichment=(Significant/Expected)) %>%
  dplyr::select(GO.ID, Term, Annotated, Significant, Expected, classicFisher=BPGSA_DiseaseGenes.score, FoldEnrichment)
#Pull significant results bonferroni correcting for number of nodes
BPs_DiseaseGenes <- BPs_DiseaseGenes %>%
  filter(classicFisher < (0.05))

# 4. Mouse KO:
#Run GSOA for mouse gene orthologs
DiseaseMouseGenes_topGO <- inner_join(mouseproteincodinggenes[c('external_gene_name', 'ensembl_gene_id')], 
                                      HumanGenes_Mouseortholog[c('Mouse.Symbol')], by =c('external_gene_name'='Mouse.Symbol'))

DiseaseMouseGenes_topGO <- DiseaseMouseGenes_topGO %>%
  mutate('EnsID'=ensembl_gene_id) %>%
  filter(!is.na(EnsID))

genesOfInterest <- DiseaseMouseGenes_topGO$EnsID

all_genes <-
  unique(as.character(mouseproteincodinggenes$ensembl_gene_id))
geneUniverse <- rep(0, length(all_genes))
names(geneUniverse) <- all_genes
geneUniverse[names(geneUniverse) %in% genesOfInterest] <- 1

GOBPdata_DiseaseGenes_mice <- new(
  "topGOdata",
  description = "testdata",
  ontology = "BP",
  allGenes = geneUniverse,
  geneSel = function(p)
    p == 1,
  annot = annFUN.org,
  ID = "ensembl",
  mapping = "org.Mm.eg.db"
)

# Run Overrepresentation Test
BPGSA_DiseaseGenes_mice <-
  runTest(GOBPdata_DiseaseGenes_mice,
          algorithm = "parentchild",
          statistic = "fisher")

#Note p-value threshold is based on bonferroni of 0.05/(GO DAG topologies tests in TopGO)
testednodes <- n_distinct(GOBPdata_DiseaseGenes_mice@graph@nodes)

BPGSA_classicfisher.results_DiseaseGenes_mice <-
  GenTable(GOBPdata_DiseaseGenes_mice,
           classicFisher = BPGSA_DiseaseGenes_mice,
           topNodes = testednodes)

#Pull significant results bonferroni correcting for number of nodes
#Add full p-values and calculate fold enrichment
Topscores <- as.data.frame(BPGSA_DiseaseGenes_mice@score)
write.csv(Topscores, './Results/Topscores_mm.csv')
Fullpvalues <- read.csv('./Results/Topscores_mm.csv')
BPs_DiseaseGenes_mice <- left_join(BPGSA_classicfisher.results_DiseaseGenes_mice, Fullpvalues, by=c('GO.ID'='X'))
BPs_DiseaseGenes_mice <- BPs_DiseaseGenes_mice %>%
  mutate(FoldEnrichment=(Significant/Expected)) %>%
  dplyr::select(GO.ID, Term, Annotated, Significant, Expected, classicFisher=BPGSA_DiseaseGenes_mice.score, FoldEnrichment)
#Pull significant results bonferroni correcting for number of nodes
BPs_DiseaseGenes_mice <- BPs_DiseaseGenes_mice %>%
  filter(classicFisher < (0.05 / testednodes))

#Find overlap between BPs overrepresented for human Disease genes and for mice orthologues 
HumanMouse_BPs <- inner_join(BPs_DiseaseGenes, BPs_DiseaseGenes_mice,
                             by  = "GO.ID")
saveRDS(HumanMouse_BPs, paste0('./Results/HumanMouse_BPs_',Disease,'.rds'))

#Pull genes that are assigned to significantly enriched Disease GO terms
HumanMouse_BPs.GOIDs <- c(HumanMouse_BPs$GO.ID)
Genes.in.DiseaseGOBPs <-
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = HumanMouse_BPs.GOIDs,
    columns = c("ALIAS", "ENSEMBL", "ENTREZID"),
    keytype = "GOALL"
  )

Genes.in.DiseaseGOBPs <- Genes.in.DiseaseGOBPs %>%
  mutate(EntrezID=as.integer(ENTREZID)) %>%
  dplyr::select("ALIAS", "EntrezID") %>%
  distinct(EntrezID, .keep_all = TRUE)
Genes.in.DiseaseGOBPs <- inner_join(Genes.in.DiseaseGOBPs, hsproteincodinggenes[c('Human.GeneID', 'UniProtID')], by = c('EntrezID'='Human.GeneID'))
Genes.in.DiseaseGOBPs <- filter(Genes.in.DiseaseGOBPs, UniProtID!="")
saveRDS(Genes.in.DiseaseGOBPs, paste0('./PipelineData/Genes.in.DiseaseGOBPs_',Disease,'.rds'))

DiseaseGenes.in.DiseaseGOBPs <- inner_join(DiseaseGenes, Genes.in.DiseaseGOBPs[c('EntrezID')],
                                           by = c('EntrezID'='EntrezID'))
DiseaseGenes.in.DiseaseGOBPsTotal <- DiseaseGenes.in.DiseaseGOBPs %>%
  mutate(EntrezID=as.factor(EntrezID)) %>%
  distinct(EntrezID)
BPs_DiseaseGenes <- c(DiseaseGenesinBPs=n_distinct(DiseaseGenes.in.DiseaseGOBPsTotal$EntrezID), SigBPs=n_distinct(HumanMouse_BPs.GOIDs))

##################################Glitch point for pipeline automation (waiting for IMPC API for term mapping)###############################
#Identify genes that when knocked out of mice are associated (p<0.05) with Disease-related symptoms
#Note: must pull 'top level terms' from IMPC
GOtoMP_mappings_DiseaseGeneBPs <- inner_join(GOtoMP_mappings[c('curie_id', 'mapped_curie', 'mapped_label', 'top_level_mp_id', 'top_level_mp_term')], HumanMouse_BPs,
                                             by = c('mapped_curie'='GO.ID'))
#Read in results of batch query
MPToplevel <- GOtoMP_mappings_DiseaseGeneBPs %>%
  distinct(curie_id, mapped_curie, mapped_label, top_level_mp_id, top_level_mp_term)
MPToplevel <- MPToplevel  %>%
  distinct(top_level_mp_term, .keep_all = TRUE)
i <- seq(1, nrow(MPToplevel), 1)
MPToplevel <- MPToplevel  %>%
  mutate(MPNUM=paste0("MP", i)) 
MPToplevel$top_level_mp_id[MPToplevel$top_level_mp_id=='info not available'] <- 
  MPToplevel$curie_id[MPToplevel$top_level_mp_id=='info not available']

write.table(GOtoMP_mappings_DiseaseGeneBPs$curie_id, paste0('./Results/',Disease,'GeneBPs_hsmm.csv'),
            row.names = F, col.names = F, quote = F)

GO_MouseKO <- list()

for (i in 1:nrow(MPToplevel)) {
  GO_MouseKO[[paste0("MP", i)]] <- fromJSON(
    paste0('https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=top_level_mp_term_id:%22', 
           MPToplevel$top_level_mp_id[[i]], '%22AND%20p_value:%5b0%20TO%200.05%5d&rows=100000&wt=json&indent=1'),
    flatten=TRUE)
}

#Reformat lists to pull mouse gene names and phenotypic consequences in KOs
#Note must pull the number of columns in the list
nMPs <- n_distinct(GO_MouseKO)
nMP <- NULL
for (i in 1:nMPs){
  nMP[[i]] <- GO_MouseKO[[i]]$response$numFound
}

phenotype <- NULL
for (i in 1:nMPs) {
phenotype[[i]] <- GO_MouseKO[[i]]$response$docs$mp_term_name[c(1:nMP[[i]])]
}

genes <- NULL
for (i in 1:nMPs) {
  genes[[i]] <- GO_MouseKO[[i]]$response$docs$marker_symbol[c(1:nMP[[i]])]
}

MPNUM <- NULL
for (i in 1:nMPs) {
  MPNUM[[i]] <- paste0('MP',i)
}

MouseGenes_MPs <- list('phenotype'=phenotype, 'genes'=genes, 'MPNUM'=MPNUM)
#Note may have to manually search for some MPs if they do not pull via the API but have associated genes
#https://www.mousephenotype.org/data/phenotypes/MP:0005380#genesAssociations
#MouseGenes_MP3 <- read.csv('./PipelineData/embryo_phenotype.csv', header=TRUE)
#Generate function to pull results from large list of lists of Random Gene Results
MouseGenes_MP1 <- data.frame(c(phenotype=phenotype[1],  
                             genes=genes[1], 
                             MPNUM=MPNUM[1]),
                             stringsAsFactors = FALSE)
MouseGenes_MP2 <- data.frame(c(phenotype=phenotype[2],  
                             genes=genes[2], 
                             MPNUM=MPNUM[2]),
                             stringsAsFactors = FALSE)
MouseGenes_MP3 <- data.frame(c(phenotype=phenotype[3],  
                             genes=genes[3], 
                             MPNUM=MPNUM[3]),
                             stringsAsFactors = FALSE)

MouseKODiseasetraits <- rbind(MouseGenes_MP1, MouseGenes_MP3)

MouseKODiseasetraits <- distinct(MouseKODiseasetraits)

DiseaseBPKOTraitGeneMouseOrthologs <- inner_join(MouseKODiseasetraits, 
                                                 HumanGenes_Mouseortholog,
                                                 by = c('genes'='Mouse.Symbol'))
DiseaseBPKOTraitGeneMouseOrthologs <- DiseaseBPKOTraitGeneMouseOrthologs %>%
  filter(Rank=='high')
MouseKODiseasetraits <- distinct(MouseKODiseasetraits)
MouseKODiseasetraits <- left_join(MouseKODiseasetraits, MPToplevel, by='MPNUM')
AllMouseGenes_DiseaseTraits <- left_join(MouseKODiseasetraits, HumanGenes_Mouseortholog,
                                         by = c('genes'='Mouse.Symbol'))
AllMouseGenes_DiseaseTraits <- AllMouseGenes_DiseaseTraits %>%
  filter(!is.na(Human.Symbol))
saveRDS(AllMouseGenes_DiseaseTraits, paste0('./PipelineData/AllMouseGenes_',Disease,'Traits.rds'))

#Identify Disease risk genes that have a relevant Disease-related phenotype when knocked out of mice
MouseKODiseasetraits.hsortho <- inner_join(AllMouseGenes_DiseaseTraits, DiseaseGenes.in.DiseaseGOBPs, 
                                           by=c('Human.Gene.ID'='EntrezID'))
MouseKODiseasetraits.hsortho <- MouseKODiseasetraits.hsortho %>%
  filter(!is.na(Human.Gene.ID)) %>%
  distinct()
saveRDS(MouseKODiseasetraits.hsortho, paste0('./Results/',Disease,'GenesMouseKO',Disease,'traits.rds'))

DiseaseGenes.MouseKOPhenotypes <- MouseKODiseasetraits.hsortho
DiseaseGenes.MouseKOPhenotypes.summary <- MouseKODiseasetraits.hsortho %>%
  mutate(MPNUM=as.factor(MPNUM), mp_id=as.factor(curie_id)) %>%
  group_by(genes) %>%
  mutate(nMouseKOTraitsperGene = n_distinct(mp_id)) %>%
  ungroup() %>%
  group_by(MPNUM) %>%
  mutate(nGenesbyKOTrait = n_distinct(genes)) %>%
  ungroup()

DiseaseGenes.MouseKOPhenotypes.totals <- MouseKODiseasetraits.hsortho %>%
  summarise(n = n_distinct(genes))

# 5. ACMG-PPIs
DiseaseGenesnames <- DiseaseGenes %>%
  dplyr::select(EntrezID) %>%
  distinct()
ACMG_DiseaseGenes <-rbind(ACMGGenes, DiseaseGenesnames)
ACMG_DiseaseGenes <- distinct(ACMG_DiseaseGenes, EntrezID)
DiseaseACMG_mapped <- string_db$map(ACMG_DiseaseGenes, "EntrezID", removeUnmappedRows = TRUE )
ACMG_Disease_PPI <- string_db$get_interactions(DiseaseACMG_mapped$STRING_id)
ACMG_Disease_PPI_mediumconfidence <- ACMG_Disease_PPI %>%
  filter(combined_score>=400)
ACMG_Disease_PPI_mediumconfidence.from <- left_join(ACMG_Disease_PPI_mediumconfidence, 
                                                    DiseaseACMG_mapped, by=c('from'='STRING_id'))
ACMG_Disease_PPI_mediumconfidence <- left_join(ACMG_Disease_PPI_mediumconfidence.from, 
                                               DiseaseACMG_mapped, by=c('to'='STRING_id'))
ACMG_Disease_PPI_mediumconfidence <- ACMG_Disease_PPI_mediumconfidence %>%
  mutate(from_gene=EntrezID.x, to_gene=EntrezID.y) %>%
  dplyr::select(from_gene, to_gene, combined_score) %>%
  distinct()

ACMG_Disease_PPI_mediumconfidence$from_Diseasegene <- ACMG_Disease_PPI_mediumconfidence$from_gene %in% DiseaseGenesnames$EntrezID
ACMG_Disease_PPI_mediumconfidence$to_ACMGgene <- ACMG_Disease_PPI_mediumconfidence$to_gene %in% ACMGGenes$EntrezID
ACMG_Disease_PPI_mediumconfidence$from_ACMGgene <- ACMG_Disease_PPI_mediumconfidence$from_gene %in% ACMGGenes$EntrezID
ACMG_Disease_PPI_mediumconfidence$to_Diseasegene <- ACMG_Disease_PPI_mediumconfidence$to_gene %in% DiseaseGenesnames$EntrezID
ACMG_Disease_PPI_mediumconfidence$Disease.ACMGconnection <- (ACMG_Disease_PPI_mediumconfidence$from_Diseasegene=='TRUE' & ACMG_Disease_PPI_mediumconfidence$to_ACMGgene=='TRUE') | (ACMG_Disease_PPI_mediumconfidence$from_ACMGgene=='TRUE' & ACMG_Disease_PPI_mediumconfidence$to_Diseasegene=='TRUE')

ACMG_Disease_PPI_mediumconfidence <- ACMG_Disease_PPI_mediumconfidence %>%
  filter(Disease.ACMGconnection)

Disease_Connection.Count.from <- ACMG_Disease_PPI_mediumconfidence %>%
  mutate(DiseaseGene=from_gene) %>%
  dplyr::select(DiseaseGene)
Disease_Connection.Count.to <- ACMG_Disease_PPI_mediumconfidence %>%
  mutate(DiseaseGene=to_gene) %>%
  dplyr::select(DiseaseGene)
Disease_Connection.Count <- rbind(Disease_Connection.Count.from, Disease_Connection.Count.to)
DiseasetoACMG_PPIs <- Disease_Connection.Count %>%
  group_by(DiseaseGene) %>%
  summarise(nConnections=n()) %>%
  ungroup()

#Average number of connections made by genes in disease set
DiseasetoACMG_PPIs.MeanSD <- summarise(DiseasetoACMG_PPIs, DiseaseGene_MeanConnectionsACMG=mean(nConnections), DiseaseGene_SDofMeanConnectionsACMG=sd(nConnections))

# Save genes in the list as character vector
DiseaseGene.Results[["DiseaseList_"]][['gene_list']] <-
  DiseaseGenes

# Summaries
DiseaseGene.Results[["DiseaseList_"]][['summary']][['gene_expression']] <-
  DiseaseGenes.Brain.summary
DiseaseGene.Results[["DiseaseList_"]][['totals']][['gene_expression']] <-
  DiseaseGenes.Brain.totals
DiseaseGene.Results[["DiseaseList_"]][['totals']][['GO_BPs']] <-
  BPs_DiseaseGenes
DiseaseGene.Results[["DiseaseList_"]][['summary']][['mouse_KO']] <-
  DiseaseGenes.MouseKOPhenotypes.summary
DiseaseGene.Results[["DiseaseList_"]][['totals']][['mouse_KO']] <-
  DiseaseGenes.MouseKOPhenotypes.totals
DiseaseGene.Results[["DiseaseList_"]][['summary']][['drug_ontology']] <-
  DiseaseGenes.DrugTargets.summary
DiseaseGene.Results[["DiseaseList_"]][['summary']][['PharmGKB']] <-
  DiseaseGenes.PharmGKB.summary
DiseaseGene.Results[[paste0("DiseaseList_")]][['summary']][['ACMG_PPIs']] <-
  DiseasetoACMG_PPIs

# Individual gene annotations (merge all possible annotations. Example: RandomGenes.DrugTargets and others)
DiseaseGene.Results[["DiseaseList_"]][['individual_annotations']][['gene_expression']] <-
  DiseaseGenes.Brain
DiseaseGene.Results[["DiseaseList_"]][['individual_annotations']][['Disease.GOBPs']] <-
  DiseaseGenes.in.DiseaseGOBPs
DiseaseGene.Results[["DiseaseList_"]][['individual_annotations']][['mouse_KO']] <-
  DiseaseGenes.MouseKOPhenotypes
DiseaseGene.Results[["DiseaseList_"]][['individual_annotations']][['drug_ontology']] <-
  DiseaseGenes.DrugTargets_proteins
DiseaseGene.Results[["DiseaseList_"]][['individual_annotations']][['PharmGKB']] <-
  DiseaseGenes.PharmGKB
DiseaseGene.Results[[paste0("DiseaseList_")]][['individual_annotations']][['ACMG_PPIs']] <-
  DiseasetoACMG_PPIs.MeanSD

saveRDS(DiseaseGene.Results, paste0('./Results/',Disease,'Gene.Results_',Date,'.rds'))

#Clear workspace
rm(list = ls(all.names = TRUE))
