###counting GRPs in the human proteome

proteins <- read.csv("/Volumes/Sandisk/human_oneline.txt", sep = "\t",
                     header = FALSE,
                     col.names = c("Gene", "Length", "GeneID", "TranscriptID", "Name", "AAseq"))

#aa freq in proteome
library(stringr)
# aa_list <- unlist(str_split(proteins$AAseq, pattern = ""))
# aa_freq <- as.data.frame(prop.table(table(aa_list)))
# #probability of sequence ocurring randomly in protein of x length: 
# #P[R or P or K]*P[G]*P[R]*P[P]*P[R or P or K]*length 
# aa_random <- (aa_freq$Freq[aa_freq$aa_list=="R"]+
#                  aa_freq$Freq[aa_freq$aa_list=="P"]+
#                  aa_freq$Freq[aa_freq$aa_list=="K"])*
#                  (aa_freq$Freq[aa_freq$aa_list=="G"]*
#                   aa_freq$Freq[aa_freq$aa_list=="R"]*
#                   aa_freq$Freq[aa_freq$aa_list=="P"])*
#                   (aa_freq$Freq[aa_freq$aa_list=="R"]+
#                   aa_freq$Freq[aa_freq$aa_list=="P"]+
#                   aa_freq$Freq[aa_freq$aa_list=="K"])
# 
# 
# #protein lengths
# proteins$lengths <- str_count(proteins$AAseq)
# #multiply
# proteins$expected_grp_counts <- proteins$lengths*aa_random

#how many proteins have how many matches to grp histogram
library(dplyr)
# proteins <- proteins %>%
#   mutate(result= +(str_count(AAseq, "[RPK]GRP[RPK]")))
library(ggplot2)
# log_freq_grpcount <- data.frame(log10(table(proteins$result)))
# ggplot(data=log_freq_grpcount, aes(x=Var1, y=Freq, )) +
#   geom_bar(stat="identity", color="black", fill="grey") +xlab("Occurrences of [RPK]GRP[RPK]")+ylab("Log10 Proportion of Proteins")+theme_classic()


######################################
proteins1 <- read.csv("/Volumes/Sandisk/uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2022.11.10-22.12.05.62.tsv", sep = "\t",
                     header = TRUE)
#IF NO CANONICAL NAME, DROP
proteins <- proteins1[!is.na(proteins1$Gene.Names..primary.) & proteins1$Gene.Names..primary. != "",]
#select longest isoform for each
proteins <- proteins1 %>% group_by(Gene.Names..primary.) %>% 
  filter(!is.na(Gene.Ontology..cellular.component.) & Gene.Ontology..cellular.component. != "") %>% 
  slice(which.max(Length))
#nogo cellular compartment
hpa <- read.csv("/Volumes/Sandisk/has_protein_data_in_Cell.tsv", sep = "\t")
hpa1<-hpa  %>% select(-contains(c("RNA" , "prognostics" , "Blood")))
hpa1<- hpa1[!(is.na(hpa1$Uniprot) | hpa1$Uniprot==""), ]
nogo <- proteins1[which(proteins1$Gene.Ontology..cellular.component. == "") , ]
nogo1 <- inner_join(hpa1,nogo,  by = c("Gene" = "Gene.Names..primary."))
nogo1[nogo1$Gene =="SRRM5","Subcellular.location"] = "Nuclear"
nogo1_filt <- nogo1[grepl("Nuclear",nogo1$Subcellular.location), ]
nogo1_filt <- nogo1_filt[!duplicated(nogo1_filt[c("Gene")]), ]

#filter non-nuclear by GO, some have multiple go categories so ill split and consider if they have nuclear
library(tidyr)
proteins <- proteins %>% 
   mutate(Gene.Ontology..cellular.component. = strsplit(as.character(Gene.Ontology..cellular.component.), ";")) %>%
  unnest(Gene.Ontology..cellular.component.)
proteins_filt<- proteins[grepl("nuclear|chromat|complex|splic|DNA|RNA|histone|transcription|nucleus", proteins$Gene.Ontology..cellular.component., ignore.case=TRUE),]
proteins_filt<- proteins_filt[!grepl("glycan|transport|pore|mito|membrane|proteasome|ribosome|
                          phosphatase|kinase|integrin|peptid|glyco|phospho|adhes|muscle|neur|
                                    lamin|immuno|tRNA|rRNA|acety|glycero|acti|cyto
                           ", proteins_filt$Gene.Ontology..cellular.component., ignore.case=TRUE),]
proteins_filt <- proteins_filt[!duplicated(proteins_filt[c("Gene.Names..primary.")]),]
# total_proteins
names(nogo1_filt)[names(nogo1_filt) == 'Gene'] <- 'Gene.Names..primary.'
nogo2 <- nogo1_filt %>% select(c("Entry", "Entry.Name", "Coiled.coil", "Domain..FT.", "Motif",
               "Protein.families", "Zinc.finger","Length", "DNA.binding", 
                "Gene.Names..primary.", "Gene.Ontology..cellular.component.", "Sequence"))
proteins_filt2 <- proteins_filt  %>% select(c("Entry", "Entry.Name", "Coiled.coil", "Domain..FT.", "Motif",
                                             "Protein.families", "Zinc.finger","Length", "DNA.binding", 
                                             "Gene.Names..primary.", "Gene.Ontology..cellular.component.", "Sequence"))
total <- rbind(nogo2, proteins_filt2)

#do they have grp
library(tidyr)
x = total %>%
  mutate(result= str_extract_all(Sequence, ".{5}GRP.{5}")) %>%
  unnest(result, keep_empty = TRUE) %>%
  mutate(location = str_locate(Sequence, result)) %>% 
  drop_na(result) %>%
  mutate(counts = str_count(result, "G|R|P|K")) %>%
  filter(counts >= 7) 

#  unnest(result, keep_empty = TRUE) %>% drop_na(result) %>%
#filter
# x<-x %>% mutate(location = gregexpr(pattern ="[RPK]GRP[RPK]",Sequence))
# y<- x %>% select(c("Entry", "Entry.Name", "Coiled.coil", "Domain..FT.", "Motif", "Protein.families", "Zinc.finger","Length", "DNA.binding", "result", "Gene.Names..primary.", "Gene.Ontology..cellular.component."))
# y$coils <- gsub("*/evidence=ECO:0000256*|SAM:Coils; | /evidence=ECO:0000255; | *COILED* | */evidence=ECO:0000255* | *SAM:Coils*", "", y$Coiled.coil)         
# y<- y %>% mutate(AT_hook_annot = str_detect(DNA.binding, 'A.T hook'))
# y <- y[!grepl("envelope|calcium|cytoplasm|interleukin|potassium|kinesin|sacchar|ubiquitin|glutamate|connexin", y$Gene.Ontology..cellular.component., ignore.case=TRUE), ]
# df <- apply(y,2,as.character)
# write.table(df, row.names = F, file = "/Volumes/Sandisk/potential_athooks_29112022.csv", sep = "\t", quote=FALSE)

#batch wget from alphafold lol
for (i in x$Entry) {
  print(paste("https://alphafold.ebi.ac.uk/files/AF-", i, "-F1-model_v4.pdb", sep = ""))
}

#with alphafold structures
aps_names <- read.table("/Volumes/Sandisk/list_pdbs")

#ggplot lengths
library(ggplot2)
x_mod <- x
x_mod$ap <- x_mod$Entry %in% aps_names$V3
x_mod$Length[x_mod$Length > 5000 ] <- 5000
ggplot(x_mod, aes(Length, fill=ap)) +
  geom_histogram( alpha = 0.5, position = 'identity', color="black", breaks=c(seq(0, 5000, by=50))) +
  scale_x_continuous(limits=c(0, 5000), breaks=c(seq(0, 5000, by=50)), labels=c(seq(0,4950, by=50), "5000+"))+
  xlab("Length")+
  ylab("Protein Isoform Count")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept  =2650,  linetype=2,colour="red")

#how many known at hooks
x<- x %>% mutate(AT_hook_annot = str_detect(DNA.binding, 'A.T hook'))
athook<- x[x$AT_hook_annot=="TRUE"  | (x$Protein.families == "HMGA family"), ]

#protein family association?
x<- x[!grepl("gluta|chloride|cyto|calcium|golgi|macrophage|voltage|ubiquitin|collagen|phosphatase|dystrophin|connexin|channel|MHC|interleukin|receptor|lamin|junction|
                           ", x$Gene.Ontology..cellular.component., ignore.case=TRUE),]
x<- x[!grepl("collagen|kinase|ubiquitin|synapsin|sclerostin|phosphatase|phage|peptidase|MINDY|
              ribosom|dystrophin|copine|adaptor|gelsolin|oxygenase
                           ", x$Protein.families, ignore.case=TRUE),]
x$Protein.families <- sub("^$", "No Family", x$Protein.families) 
x_unique <- unique(x[, c('Gene.Names..primary.', 'Protein.families','AT_hook_annot')])
counts <- data.table(table(x_unique$Protein.families, x_unique$AT_hook_annot))
ggplot(counts, aes(N)) + geom_histogram()

#file accessibilities
# acc<- read.csv("/Volumes/Sandisk/9606.accessibility.tdt", header = F, sep = "\t")
acc_mine <- read.csv("/Volumes/Sandisk/asa_ss_29112022.tsv", header = T, sep = "\t")
joined <-left_join(x, acc_mine, by = c("Entry" = "X"))
joined_inner <-inner_join(x, acc_mine, by = c("Entry" = "X"))
#remainder accessibilities
joined_rem1 <-subset(x, !(Entry %in% joined_inner$Entry)) %>% distinct()
joined_rem<- joined_rem1 %>%
  mutate(result= str_extract_all(Sequence, "(.{5}GRP.{5})")) %>%
  unnest(c(result)) %>%
  mutate(location = str_locate(Sequence, result)) %>% 
  drop_na(result) %>%
  mutate(counts = str_count(result, "G|R|P|K")) %>%
  filter(counts >= 7) %>% distinct()
joined_rem<- joined_rem%>%  
  mutate(start_context = location[,1] - 130,
         end_context = location[,2] + 130,
         Length = as.numeric(as.character(Length))) %>%
  mutate(start_context= if_else(start_context < 0, 0, start_context),
         end_context= if_else(end_context > Length, Length, end_context)) %>%
  mutate(substring = substr(Sequence, start_context, end_context))
joined_rem<- joined_rem %>%
  group_by(Entry) %>%
  mutate(col_suffixed = paste0(Entry, "_", row_number()))

# #to fasta
# con <- file("/Volumes/Sandisk/sequences.fasta", "w")
# # Loop through the rows of the data frame
# for (i in 1:nrow(joined_rem)) {
#   # Write the ID of the sequence to the FASTA file
#   write(paste0(">", joined_rem$col_suffixed[i]), sep = "",con, append = TRUE)
#   # Write the sequence to the FASTA file
#   write(joined_rem$substring[i], con, append = TRUE)
# }
# 
# # Close the file connection
# close(con)
#accessibilities of remaining at hooks
acc_rem <- read.csv("/Volumes/Sandisk/athook_rem_acc/asa_ss_12122022.tsv", header = T, sep = "\t")
merged_accrem <- merge(acc_rem, joined_rem %>% select(col_suffixed, result, Entry, substring), by.x = "X", by.y = "col_suffixed")

joined1 <- joined %>% 
  left_join(merged_accrem, by =c("Entry", "result")) %>% 
  mutate(asa = coalesce(asa.x, asa.y),
         sec_structure = coalesce(sec_structure.x, sec_structure.y),
         sequence_grp = coalesce( substring, Sequence), 
         new_len = nchar(sequence_grp)) %>% 
  select(-asa.x, -asa.y, sec_structure.x, sec_structure.y, X)
  

#accesbilities of known AT hooks?
athook<- joined1[joined1$AT_hook_annot=="TRUE"  | (joined1$Protein.families == "HMGA family"), ]
library(tidyr)
athook <- athook %>% 
  unnest(location)


#scratch
athook <- athook %>%
  mutate(location_corr = str_locate(sequence_grp, result),
         range = paste0(location_corr[,1] , "..", location_corr[,2]))
# athook$start = lapply(sapply(str_split(athook$range,  '\\.\\.'), `[`, 1), as.numeric)
# athook$end = lapply(sapply(str_split(athook$range,  '\\.\\.'), `[`, 2), as.numeric)
# athook$substring = sapply(strsplit(athook$asa, ","), as.numeric)
# athook$act_len = lapply(athook$substring, length)
# athook$substring = data.frame(athook$sequence_grp[athook$start:athook$end])

func_name <- function (range, string, delim) {
  x = str_split(range, delim)
  start = as.numeric(x[[1]][1])
  end = as.numeric(x[[1]][2])
  substring = c(as.numeric(strsplit(string,',')[[1]]))
  print(length(substring))
  substring = data.frame(substring[start:end])
  return(substring)
}

func_namev <- Vectorize(func_name)
f = func_namev(athook$range, athook$asa, '\\.\\.')

library(data.table)
athook_acc <- rbindlist(lapply(f, function(x) data.frame(t(x))), fill = T)
thresh <- quantile( apply(athook_acc, 1, mean, na.rm=T), .50)



#putative athook accessibilities
put_athook<- joined1[joined1$AT_hook_annot!="TRUE" & joined1$Protein.families != "HMGA family",]
put_athook  <- put_athook %>% 
  unnest(location)
put_athook <- put_athook %>%
  mutate(location_corr = str_locate(sequence_grp, result),
         range = paste0(location_corr[,1] , "..", location_corr[,2]))


# joined <-inner_join(put_athook, acc_mine, by = c("Entry" = "X"))
f = func_namev(put_athook$range, put_athook$asa, '\\.\\.')
put_athook$mean_acc <- apply(rbindlist(lapply(f, function(x) data.frame(t(x))), fill = T), 1, mean) 
put_athook_acc1 = put_athook[put_athook$mean_acc > thresh,]

#joined ss
func_name_ss <- function (range, string, delim) {
  x = str_split(range, delim)
  start = as.numeric(x[[1]][1])
  end = as.numeric(x[[1]][2])
  substring = c(strsplit(string,',')[[1]])[start:end]
  # substring = data.frame(substring[start:end])
  substr = paste(substring, collapse="")
  return(substr)
}
func_name_ss_v <- Vectorize(func_name_ss)
##########known at hooks
# f = func_name_ss_v(athook$range, athook$sec_structure, '\\.\\.')
# # athook_ss <- rbindlist(lapply(f, function(x) data.frame(t(x))), fill = T)
# athook$ss <- f
athook <- athook %>% drop_na(Entry)
athook$start = unlist(lapply(sapply(str_split(athook$range,  '\\.\\.'), `[`, 1), as.numeric))
athook$end = unlist(lapply(sapply(str_split(athook$range,  '\\.\\.'), `[`, 2), as.numeric))
athook$substring = unlist(lapply(sapply(strsplit(athook$sec_structure, ","), as.character),paste, collapse=""))
athook$ss = substr(athook$substring,athook$start, athook$end)

#put at hook ss
put_athook_acc <- put_athook_acc1 %>% drop_na(Entry)
put_athook_acc$start = unlist(lapply(sapply(str_split(put_athook_acc$range,  '\\.\\.'), `[`, 1), as.numeric))
put_athook_acc$end = unlist(lapply(sapply(str_split(put_athook_acc$range,  '\\.\\.'), `[`, 2), as.numeric))
put_athook_acc$substring = unlist(lapply(sapply(strsplit(put_athook_acc$sec_structure, ","), as.character),paste, collapse=""))
put_athook_acc$ss = substr(put_athook_acc$substring,put_athook_acc$start, put_athook_acc$end)


#filter out those with more than 25% as H
put_athook_acc <- put_athook_acc %>%
  filter(str_count(ss, "H") <3, 
         str_count(ss, 'E|G|I|B|C|S')<3)



                                
