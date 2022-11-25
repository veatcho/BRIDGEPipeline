#Install libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("STRINGdb")
#BiocManager::install("topGO")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("ExpressionAtlas")

#load necessary libraries
library(tidyverse)
library(jsonlite)
library(httr)
library(biomaRt)
library(ghql)

#Identify Disease CUIs
#1) Go to MeSH browser: https://meshb.nlm.nih.gov/treeView, select entity of interest (for example 'Mental Disorders [F03]')
#2) Select disease/disorder of interest to pull and identify Unique ID under details tab here
#3) Search MedGen for unique identifier and select concept unique identifier (https://www.ncbi.nlm.nih.gov/medgen)
#4) Cross reference with manual search by typing in 'preferred name" from MeSH browser in https://www.disgenet.org/search
#5) other useful resource: https://monarchinitiative.org/; https://www.ebi.ac.uk/ols/ontologies/mondo
#In past I have also gone to OMIM for all potential CUIs from DisGeNET and manually reviewed details (see ASD, NF1, LateOnsetAD)
#for ASD="^Autism|MACROCEPHALY/AUTIS|^AUTIS|High-functioning autism"
#ASD=C1510586, C1845336, C3840214, C1854416, C3275438, C3150693, C1839707, C3552491, C2677504, C1845540 C3550875, C1845539, C3554373, C0795888, C3150677
#for NF1="^Neurofibromatosis 1|Recklinghausen|^NF1|Peripheral Neurofibromatosis"
#for LateOnsetAD=two step with negation if pulling from the download file: filter(grepl("Alz", diseaseName)) %>% filter(!grepl('Early|Familial|familial|type 1', diseaseName))
#LOAD = C0002395, C0494463, C0750900, C1735366, C1851958, C1853555, C3665464, C1970209, C1837149, C1970144
#for Schizophrenia: https://www.ncbi.nlm.nih.gov/medgen/48574, C0036341
#for Bipolar disorder: https://www.ncbi.nlm.nih.gov/medgen/?term=D001714, C0005586
#for Intellectual disability: https://www.ncbi.nlm.nih.gov/medgen/?term=D008607, C3714756|C0025363

#Annotation resources
#Check for updates to resources (poll monthly)
#Access DisGeNET data to identify all genes with disease-associated variants, last update October 2021
#Note this resource will require you to create an account and then update code included below on line 40 to include your email and password
#Gene/Disease Associations (https://www.disgenet.org/api/#/GDA/gdaByDisease)

#This example is applied to Schizophrenia
Disease <- 'Schizophrenia'
CUI <- 'C0036341'

auth_params <- list("email"="oveatch@kumc.edu","password"="BRIDGE_2021")

api_host <- "https://www.disgenet.org/api"

api_key <- NULL

r <-POST(paste(api_host,"/auth/", sep=""), body=auth_params)
if(status_code(r)==200){
  #Lets store the api key in a new variable and use it again in new requests
  json_response <- content(r, "parsed")
  api_key <- json_response$token
  print(paste(api_key," This is your user API key.",sep="")) #Comment this line if you don't want your API key to show up in the terminal
}else{
  print(status_code(r))
  print(content(r, "parsed"))
}
if(!is.null(api_key)){
  #Store the api key into a named list with the following format, see next line.
  authorization_headers <- c(Authorization=paste("Bearer ",api_key, sep=""))
  #Lets get all the diseases associated with CUIs of interest
  gda_response <- GET(paste(api_host,"/gda/disease/", CUI, sep=""), query=list(source="ALL", format="json"), add_headers(.headers=authorization_headers)) #All protected endpoints require the add_headers and the authorization token in order to be accessed.
  GDA_parsed <- print(content(gda_response, "parsed"))
}

#Pull elements of list that correspond to results of interest and output in a dataframe
geneid <- unlist(sapply(GDA_parsed, "[[", 1))
gene_symbol <- unlist(sapply(GDA_parsed, "[[", 2))
uniprotid <- sapply(GDA_parsed, "[[", 3)
is.na(uniprotid) <- lengths(uniprotid) ==0
uniprotid <- unlist(uniprotid)
diseaseid <- unlist(sapply(GDA_parsed, "[[", 9))
disease_name <- unlist(sapply(GDA_parsed, "[[", 10))
score <- unlist(sapply(GDA_parsed, "[[", 15))
diseasecuis <- data.frame(
  'geneid'=geneid, 'gene_symbol'=gene_symbol, 'uniprotid'=uniprotid, 'diseaseid'=diseaseid, 
  'disease_name'=disease_name, 'score'=score
  )
diseasecuis <- diseasecuis %>%
  arrange(desc(score)) %>%
  distinct(geneid, .keep_all = TRUE)
  
saveRDS(diseasecuis, paste0('./PipelineData/', Disease, '.rds'))

#Update protein coding genes in humans
#Pull names for all protein coding genes known in humans that are included in Ensembl (GRCh38.p13)
humangenes <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#listFilters(humangenes)
#listAttributes(humangenes)
ensIDswithproteins<-getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprot_gn_id"), 
                          mart=humangenes, uniqueRows=TRUE)
ensIDswithproteins_distinctgenes <- ensIDswithproteins %>%
  filter(!uniprot_gn_id=="") %>%
  distinct(entrezgene_id, .keep_all=TRUE)

hsproteincodinggenes <- ensIDswithproteins_distinctgenes %>%
  mutate(Gene.Name = as.character(external_gene_name), Human.GeneID = as.integer(entrezgene_id), EnsID = as.character(ensembl_gene_id), UniProtID = as.character(uniprot_gn_id)) %>%
  dplyr::select(Gene.Name, Human.GeneID, EnsID, UniProtID)

saveRDS(hsproteincodinggenes, './PipelineData/hsproteincodinggenes.rds')

#Pull Ensembl IDs for risk gene symbols
diseasecuis <- inner_join(diseasecuis, hsproteincodinggenes[c('Gene.Name', 'EnsID')], by = c('gene_symbol'='Gene.Name'))
DiseaseGenes <- diseasecuis %>%
  mutate(EntrezID=geneid, Gene.Name=gene_symbol, UniProtID=uniprotid, DisGeNETScore=score) %>%
  filter(!is.na(UniProtID)) %>%
  distinct(Gene.Name, .keep_all = TRUE) %>%
  dplyr::select(EntrezID, Gene.Name, UniProtID, EnsID, diseaseid, disease_name, DisGeNETScore)

#save list of genes named using broad disease category of interest for pipeline
saveRDS(DiseaseGenes, paste0('./PipelineData/',Disease,'Genes.rds'))

#Export gene curation results from ClinGen for adult onset conditions use: https://actionability.clinicalgenome.org/ac/Adult/ui/summ
#for pediatric onset use: https://actionability.clinicalgenome.org/ac/Pediatric/ui/summ
#To interpret scoring criteria see: https://www.clinicalgenome.org/site/assets/files/2180/actionability_sq_metric.png
#Download American College of Medical Genetics (ACMG) list of recommended genes from: https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/
ACMGgenelist <- read.csv('./PipelineData/ACMGlist.txt', sep = ',', header = F)
ACMGgenes <- inner_join(ACMGgenelist, hsproteincodinggenes[c('Gene.Name', 'Human.GeneID')],
                        by = c('V1'='Gene.Name'))
ACMGgenes <- ACMGgenes %>%
  dplyr::select(EntrezID=Human.GeneID) %>%
  distinct()
saveRDS(ACMGgenes, './PipelineData/ACMGgenelist.rds')

###################Gene Expression Data###########################
# To download GTEx data in entirety (most recent release=Data Source: GTEx Analysis Release V8 (dbGaP Accession phs000424.v8.p2))
#library('ExpressionAtlas')
#Find dataset of interest
#Atlasdata <- searchAtlasExperiments("organism part", species = "Homo sapiens")
#Create dataframe with Accession number and experiment title
#Atlasdata_search <- data.frame(Accession=Atlasdata$Accession, Experiment=Atlasdata$Title)
#Find line number of dataset of interest
#GTExaccession_line <- agrep('GTEx', Atlasdata_search$Experiment)
#Find accession number of experiment, using line number, for which data are to be pulled
#GTExaccession <- Atlasdata_search[paste(GTExaccession_line),]
#Pull data
#GTEx <- getAtlasExperiment(paste(GTExaccession$Accession))
#GTExdata <- GTEx$rnaseq
#GTExcounts <- assays(GTExdata)$counts
#samplestuff <- colData(GTExdata)
#tissues <- data.frame(attribute=samplestuff@rownames, tissue=samplestuff@listData$organism_part)
#Extract RNA-seq counts in brain samples
#brain <- tissues %>%
#  filter(tissue=='amygdala' | tissue=='Brodmann (1909) area 24' | tissue=='Brodmann (1909) area 9' | tissue=='caudate nucleus' |
#           tissue=='cerebellar hemisphere' | tissue=='cerebellum' | tissue=='cerebral cortex' | tissue=='hippocampus proper' |
#           tissue=='hypothalamus' | tissue=='nucleus accumbens' | tissue=='pituitary gland' | tissue=='putamen' |
#           tissue=='substantia nigra')
#Extract counts where brain sample accession ids matches count data column names (i.e., sample accession ids)
#GTExbraindata <- cbind(brain, GTExcounts[match(brain$attribute, colnames(GTExcounts)),c(1:ncol(GTExcounts))])
#GTExbraindata <- GTExbraindata %>%
#  mutate(EnsID.GTEx=row.names(GTExbraindata), tissue=as.factor(tissue)) %>%
#  group_by(EnsID.GTEx, tissue) %>%
#  mutate(Average=mean(c(3:18738)))

#Download Expression value across all gene (TPM) from https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5214/Downloads on 04012021
#GTExdata <- read.table('./dowload/E-MTAB-5214-query-results.tpms.tsv', header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote="")
#GTExbraindata <- GTExdata %>%
#  dplyr::select(EnsID.GTEx=Gene.ID, amygdala, Brodmann1909.area24=Brodmann..1909..area.24, Brodmann1909.area9=Brodmann..1909..area.9, caudate.nucleus,
#         cerebellar.hemisphere, cerebellum, cerebral.cortex, hippocampus.proper, hypothalamus, nucleus.accumbens,
#         pituitary.gland, putamen, substantia.nigra)
saveRDS(GTExbraindata, './PipelineData/GTExbraindata.rds')

###################Human to Mouse Mappings##########################
#Update MP to GO mappings
#Go to https://www.ebi.ac.uk/spot/oxo/datasources/MP and click on GO then download results and save as .csv
#Ontology ID: mp
#Version: 11-11-2019
#Number of terms: 39543
#Last loaded: Fri Jan 29 11:10:47 GMT 2021
#Read in GOtoMP mappings downloaded from https://www.ebi.ac.uk/spot/oxo/search#
GOtoMP_mappings <- read.csv('./PipelineData/GOtoMP_mappings.csv', header=TRUE)
#Pull 'top level terms' from IMPC for 'mappings.csv', to find top level terms go to https://www.mousephenotype.org/data/batchQuery
#Select 'MP' as ID type and customize output to add 'top level mp id' and 'top level mp term'
GOtoMP_mappings_MPToplevel <- read.table('./PipelineData/batch_query_dataset.tsv', header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote="")
GOtoMP_mappings_MPToplevel <- inner_join(GOtoMP_mappings, GOtoMP_mappings_MPToplevel, by = c('curie_id'='mp_id'))
write.csv(GOtoMP_mappings_MPToplevel, './PipelineData/GOtoMP_mappings_MPToplevel.csv', row.names = F)
#Update protein coding genes in mice
#Pull names for all protein coding genes known in mice that are included in Ensembl (GRCm38.p6)
micegenes <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#listFilters(micegenes)
#listAttributes(micegenes)
mouseensIDswithproteins<-getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "uniprot_gn_id", "external_gene_name"), filters = list('with_protein_id'=TRUE), mart=micegenes)
mouseensIDswithproteins_distinctgenes <- distinct(mouseensIDswithproteins, entrezgene_id, .keep_all=TRUE)
saveRDS(mouseensIDswithproteins_distinctgenes, './PipelineData/mmproteincodinggenes.rds')

#The human uniprotid list was exported and batch searched using all prediction tools available via the DRSC Integrative Ortholog Prediction Tool (DIOPT; Version 7.1 March 2018; https://www.flyrnai.org/diopt) 
#with filtering to return only best match when there is more than one match per input gene or protein
hsproteincodinggenes_UniProtID <- distinct(hsproteincodinggenes, UniProtID)
write.table(hsproteincodinggenes_UniProtID$UniProtID, './hsproteins.txt', row.names = F, col.names = F, quote = F)

#Pull results from human/mouse ortholog query using DIOPT, version 8.0 (Aug 2019)
#Ran on uniprot ids using 'Return only best match when there is more than one match per input gene or protein' filter on May 11, 2021 
#Results:
#18972  query symbols mapped to 18580  genes
#210 of them had no orthologs
Humangenes.msortho <- read.table(file = './PipelineData/hsmmorthologs_DIOPTquery051121.txt', 
                                 header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
#Remove low or none ranks
Humangenes.msortho_distinct <- Humangenes.msortho %>%
  filter(Rank!='low' & Rank!='None') %>%
  arrange(desc(Weighted.Score)) %>%
  distinct(Human.Symbol, .keep_all = TRUE)
saveRDS(Humangenes.msortho_distinct, paste0('./PipelineData/HumanGeneMouseOrthologs_',Date,'.rds'))
#If issues with time for DIOPT can also perform analyses to identify orthologs of disease genes instead
#write.table(DiseaseGenes$UniProtID, paste0('./PipelineData/hsproteincodinggenes_',Disease,'.csv'),
    #        sep = '\n', row.names = FALSE, col.names = FALSE, quote = FALSE)

###################Pharos##########################
#Drug Ontology (current release v3.4.3, Update: 2020-02-26)
token <-Sys.getenv("GITHUB_GRAPHQL_TOKEN")
con <- GraphqlClient$new(url = 'https://pharos-api.ncats.io/graphql')

#Perform query to pull gene symbols for Tclin, Tchem targets
qry <- Query$new()

qry$query('mydata', '{
  targets(filter:{facets:[{
    facet: "Target Development Level",
    values: ["Tclin", "Tchem"]
  }]
  })
  {
    targets(top: 18319) {
      sym
      tdl}
  }
}')

x <- con$exec(qry$queries$mydata)
idg <- jsonlite::fromJSON(x)
drugtargets <- data.frame(Gene.Name=idg$data$targets$targets$sym, idgTDL=idg$data$targets$targets$tdl)
drugtargets <- inner_join(drugtargets, hsproteincodinggenes[c('Gene.Name', 'Human.GeneID', 'UniProtID')],
                          by = 'Gene.Name')
saveRDS(drugtargets, './PipelineData/drugtargets.rds')

###################PharmGKB##########################
#Pull data from PharmGKB  , downloaded the 'var_drug_ann.tsv' file which contains
#associations in which the variant affects a drug dose, response, metabolism, etc
#file last modified March 5, 2021
VariantDrugAssociations <- read.table('./PipelineData/var_drug_ann.tsv', header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote="")
VariantDrugAssociations.1 <- VariantDrugAssociations %>%
  filter(!grepl(',', Gene))
#Find variants with more than one gene and reformat
#Count number of commas to determine the largest number of splits required
VariantDrugAssociations.multi <- VariantDrugAssociations %>%
  filter(grepl(',', Gene))
n_pharmgkbids <- str_count(VariantDrugAssociations.multi$Gene, ',')
n_pharmgkbids <- max(n_pharmgkbids)

for (i in 1:n_pharmgkbids) {
  
  VariantDrugAssociations.multi <- separate(VariantDrugAssociations.multi, col = 'Gene', 
                                            into = c(paste0('Gene.Symbol_PharmGKBID_', i)),
                                            sep = ",", remove = FALSE)
}

VariantDrugAssociations.m <- NULL

for (i in 1:n_pharmgkbids) {
  VariantDrugAssociations.m[[i]] <- VariantDrugAssociations.multi %>%
    filter(!is.na(paste0('Gene.Symbol_PharmGKBID_',i))) %>%
    mutate(Gene=paste0('Gene.Symbol_PharmGKBID_',i)) %>%
    dplyr::select('Annotation.ID', 'Variant', 'Gene',
                  'Chemical', 'PMID', 'Phenotype.Category', 'Significance', 'Notes', 'Sentence', 'StudyParameters',
                  'Alleles', 'Chromosome')
}

VariantDrugAssociations.m <- bind_rows(VariantDrugAssociations.m)

VariantDrugAssociations <- rbind(VariantDrugAssociations.m, VariantDrugAssociations.1)

VariantDrugAssociations <- separate(VariantDrugAssociations, col = 'Gene', 
                                    into = c('Gene.Symbol', 'PharmGKBID'), sep = " ", remove = TRUE)
VariantDrugAssociations$Gene.Symbol <- gsub('"', "", VariantDrugAssociations$Gene.Symbol)

VariantDrugAssociations <- left_join(VariantDrugAssociations, hsproteincodinggenes, by=c('Gene.Symbol'='Gene.Name'))
VariantDrugAssociations <- VariantDrugAssociations %>%
  filter(!is.na(PharmGKBID)) %>%
  distinct()
saveRDS(VariantDrugAssociations, './PipelineData/PharmGKB.rds')

#Clear workspace
rm(list = ls(all.names = TRUE))