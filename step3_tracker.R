# sudo apt-get install libudunits2-dev
# sudo apt-get install libgdal-dev
# Install development version from GitHub
# devtools::install_github("outbreak-info/R-outbreak-info")

library(jsonlite)
library(plotly)
library(ggplot2)
library(readr)
library(reactable)
library(outbreakinfo)
library(dplyr)
library(Biostrings)


# run mutation.ipynb
# run ~/pyro-cov/backtest_usher.ipynb derived from https://github.com/broadinstitute/pyro-cov
# download WHO weekly report https://www.who.int/emergencies/diseases/novel-coronavirus-2019/situation-reports

recent = 30  # define 'recent' days to look back
min_ratio = 0.02 # minimum threshold of prevalent variant/lineage
min_freq = 0.1  # minimum threshold of mutation frequency in a variant/lineage
data_folder='COVID/data/'
this_run=Sys.Date()

who_date='2023-04-27'
fitness_result='~/pyro-cov/backtest/result_2022-12-22/strains.tsv'
mutation_result='~/pyro-cov/backtest/result_2022-12-22/mutations.tsv'

last_run = tryCatch(
  expr = {readRDS(paste0(data_folder, '/recent_run.RDS'))},
  error = function(e){'NO_DATE'}
)

saveRDS(this_run, paste0(data_folder, '/recent_run.RDS'))

#authenticateUser()


## Variants of concern (VOCs)

# Retrieve a curated list of Variants of Concern, Variants of Interest, 
# Variants Under Monitoring, and de-escalated variants maintained by 
# the [outbreak.info team](https://outbreak.info/situation-reports)
voc_report = getCuratedLineages()
voc = voc_report[voc_report$variantType %in% c('Variant of Concern','Previously Circulating Variant of Concern'),]

# get the earliest VOC designation from multiple institutes
voc$voc_name=ifelse(is.na(voc$short_name), voc$who_name, voc$short_name)

voc$earliest_designation_date=unlist(lapply(voc$classifications, function(x){
  min(x$dateModified[x$variantType=='VOC'], na.rm = TRUE)
}))

outbreak_date<-as.character(as.Date(voc$classificationTable$VOC$outbreak$label, format="%d %b %Y"))
outbreak_date<-pmax(voc$earliest_designation_date, outbreak_date, na.rm = TRUE)
voc$earliest_designation_date<-outbreak_date

# voc$earliest_designation_date[grep('BA.5', voc$voc_name)]='2022-07-20'
# voc$earliest_designation_date[grep('BA.2.75', voc$voc_name)]='2022-07-20'
# voc$earliest_designation_date[grep('XBB', voc$voc_name)]='2022-11-10'
# voc$earliest_designation_date[grep('BQ.1', voc$voc_name)]='2022-11-10'
# voc$earliest_designation_date[grep('BF.7', voc$voc_name)]='2022-11-10'
# voc$earliest_designation_date[grep('XBB.1.5', voc$voc_name)]='2023-01-09'
# voc$earliest_designation_date[grep('XBB.1.16', voc$voc_name)]='2023-01-09'


voc_display =voc[,c("voc_name", "pangolin_lineage", "nextstrain_clade",
                    "gisaid_clade",  "location_first_identified",
                    "earliest_designation_date", "variantType")]

reactable(voc_display, 
          defaultSorted = "earliest_designation_date", defaultSortOrder="desc")

saveRDS(voc_report, 
        paste0(data_folder, '/voc_report_', Sys.Date(), '.RDS'))



voc_prevalence=NULL

for(i in 1:nrow(voc)){
  des=voc$pango_descendants[i][[1]]
  #des=des[!grepl('^X',des)]
  size=length(des)
  round=floor(size/500)
  prevalence<-NULL
   for(j in 0:round){
     if(j==round){
       prevalence_temp=getPrevalence(
         pangolin_lineage=paste(des[(j*500+1):size],collapse=' OR '),
         logInfo = FALSE
       )
     }else{
       prevalence_temp=getPrevalence(
         pangolin_lineage=paste(des[(j*500+1):(j*500+500)],collapse=' OR '),
         logInfo = FALSE
       )
     }
     prevalence<-rbind(prevalence, prevalence_temp)
   }
  prevalence$lineage=voc$voc_name[i]
  prevalence_sum<-prevalence %>% 
    select(date, lineage, proportion) %>%
    group_by(date, lineage) %>%
    summarise(proportion_s=sum(proportion))
  
  voc_prevalence=rbind(voc_prevalence, prevalence_sum)
  
}
voc_prevalence=voc_prevalence %>% dplyr::rename(proportion=proportion_s) %>% ungroup()


plot_ly(voc_prevalence, x = ~date, y = ~proportion, color = ~lineage) %>%
  add_lines() %>% layout(title = 'Daily prevalance of the designated variants of concern')

saveRDS(voc_prevalence, 
        paste0(data_folder, '/voc_prevalence_', Sys.Date(), '.RDS'))

## Recently prevalent VOC(s)


voc_prevalence_recent =voc_prevalence[
  voc_prevalence$date>=Sys.Date()-recent, ]

prevalent_voc=unique(
  voc_prevalence_recent$lineage[voc_prevalence_recent$proportion>min_ratio])

for(i in 1:length(prevalent_voc)){
  print(prevalent_voc[i])
}

prevalent_voc = prevalent_voc[!is.na(prevalent_voc)]

#if(as.character(last_run)=='NO_DATE'){
#  prevalent_voc_last=NULL
#}else{
#  prevalent_voc_last=readRDS(paste0(data_folder, '/voc_recent_', last_run, '.RDS'))
#}

#print(setdiff(prevalent_voc, prevalent_voc_last))


#print(setdiff(prevalent_voc_last, prevalent_voc))



# get the mutations of the recently prevalent lineage(s) among VOCs
prevalent_mutation=NULL
sub_lins=NULL

for(i in 1:length(prevalent_voc)){
  des=voc$pango_descendants[voc$voc_name==prevalent_voc[i]][[1]]
  #des=des[!grepl('^X',des)]
  
  size=length(des)
  round=floor(size/500)
  mutations<-NULL
  for(j in 0:round){
    if(j==round){
      mutations_temp=getMutationsByLineage(
        pangolin_lineage=paste(des[(j*500+1):size],collapse=' OR '),
        frequency=min_freq,
        logInfo = FALSE)
    }else{
      mutations_temp=getMutationsByLineage(
        pangolin_lineage=paste(des[(j*500+1):(j*500+500)],collapse=' OR '),
        frequency=min_freq,
        logInfo = FALSE
      )
    }
    mutations<-rbind(mutations, mutations_temp)
  }
  
  mutations$lineage=prevalent_voc[i]
  
  mutations_sum<-mutations %>% 
    select(mutation, codon_num, gene, lineage, mutation_count, lineage_count) %>%
    group_by(lineage, gene, mutation, codon_num) %>%
    summarise(mutation_count_s=sum(mutation_count),
              lineage_count_s=sum(lineage_count))
  
  prevalent_mutation=rbind(prevalent_mutation, mutations_sum)
  sub_lins=c(sub_lins, des)
}

prevalent_mutation <- prevalent_mutation %>% 
  mutate(prevalence=mutation_count_s/lineage_count_s) %>%
  select(-c(mutation_count_s,lineage_count_s)) %>%
  ungroup()


sub_lins=unique(sub_lins)
prevalent_mutation_sub1 = getMutationsByLineage(pangolin_lineage=sub_lins, 
                                                frequency=min_freq, logInfo = FALSE)
prevalent_mutation_sub1 <- prevalent_mutation_sub1 %>% 
  select(lineage, gene, mutation, codon_num, prevalence)

prevalent_mutation_sub1_s = rbind(prevalent_mutation_sub1[prevalent_mutation_sub1$gene=='S',], 
                                  prevalent_mutation[prevalent_mutation$gene=='S',])
prevalent_mutation_sub1_s$mutation = toupper(prevalent_mutation_sub1_s$mutation)

plot_ly(prevalent_mutation_sub1_s,
        x=~reorder(mutation, codon_num), y=~lineage, z=~prevalence, 
        type = "heatmap", colors = "Greys",
        height= 18*length(unique(prevalent_mutation_sub1_s$lineage))+100,
        width = 12*length(unique(prevalent_mutation_sub1_s$mutation))+100) %>%
  layout(margin=list(l=50, r=50, t=50, b=50, pad = 4),showlegend=FALSE,
         title = 'Mutation frequency of the recently prevalent VOC(s) and all its sublineages')


saveRDS(prevalent_mutation_sub1_s, 
        paste0(data_folder, '/voc_lineage_mutation_', Sys.Date(), '.RDS'))

# get the recently prevalent sublineage(s) among recent VOCs

sub_prevalence=NULL
for(i in 1:length(prevalent_voc)){
  prevalence=getPrevalence(
    pangolin_lineage=voc$pango_descendants[voc$voc_name==prevalent_voc[i]][[1]],
    logInfo = FALSE
  )
  prevalence$who_name=prevalent_voc[i]
  sub_prevalence=rbind(sub_prevalence, prevalence)
}

# prevalent sublineage(s) in most recent week
sub_prevalence_recent =sub_prevalence[sub_prevalence$date>=Sys.Date()-recent, ]

prevalent_sub=unique(
  sub_prevalence_recent$lineage[sub_prevalence_recent$proportion >min_ratio])

prevalent_sub_df=sub_prevalence[sub_prevalence$lineage %in% prevalent_sub,]

prevalent_sub_df_unique=unique(prevalent_sub_df[,c("lineage","who_name")])

reactable(prevalent_sub_df_unique)

plot_ly(prevalent_sub_df, x = ~date, y = ~proportion, color = ~lineage) %>%
  add_lines() %>% layout(title = 'Daily prevalance of the recently prevalent sublineages')

saveRDS(prevalent_sub_df_unique, 
        paste0(data_folder, '/voc_sublineage_recent_', Sys.Date(), '.RDS'))
saveRDS(prevalent_sub_df, 
        paste0(data_folder, '/voc_sublineage_prevalence_', Sys.Date(), '.RDS'))

  
#  if(as.character(last_run)=='NO_DATE'){
#    prevalent_sub_last=NULL
#  }else{
#    prevalent_sub_last=readRDS(paste0(data_folder, '/voc_sublineage_recent_', last_run, '.RDS'))
#  }
#print(setdiff(prevalent_sub_df_unique$lineage, prevalent_sub_last$lineage))

#print(setdiff(prevalent_sub_last$lineage, prevalent_sub_df_unique$lineage))


# get the mutations of the recently prevalent sublineage(s) among recent VOCs
prevalent_mutation_sub2_s=prevalent_mutation_sub1_s[prevalent_mutation_sub1_s$lineage %in% prevalent_sub,]

prevalent_mutation_sub2_s = rbind(prevalent_mutation_sub2_s[prevalent_mutation_sub2_s$gene=='S',], 
                                  prevalent_mutation[prevalent_mutation$gene=='S',])
prevalent_mutation_sub2_s$mutation = toupper(prevalent_mutation_sub2_s$mutation)

prevalent_mutation_sub2_s_max_freq = aggregate(prevalent_mutation_sub2_s$prevalence, 
                                               by=list(prevalent_mutation_sub2_s$mutation), max)

prevalent_mutation_sub2_s<-prevalent_mutation_sub2_s[prevalent_mutation_sub2_s$mutation %in% prevalent_mutation_sub2_s_max_freq$Group.1[prevalent_mutation_sub2_s_max_freq$x>min_freq],]

plot_ly(prevalent_mutation_sub2_s,
        x=~reorder(mutation, codon_num), y=~lineage, z=~prevalence, 
        type = "heatmap", colors = "Greys",
        height= 25*length(unique(prevalent_mutation_sub2_s$lineage))+100,
        width = 12*length(unique(prevalent_mutation_sub2_s$mutation))+100) %>%
  layout(margin=list(l=50, r=50, t=50, b=50, pad = 4),showlegend=FALSE,
         title = 'Mutation frequency of the recently prevalent sublineages')

saveRDS(prevalent_mutation_sub2_s, 
        paste0(data_folder, '/voc_sublineage_mutation_', Sys.Date(), '.RDS'))


fitness_lineage=read.delim(fitness_result, header=T, stringsAsFactors = F)
fitness_lineage=fitness_lineage[,c(1,2,5)]
names(fitness_lineage)=c('rank','lineage','fitness score')
saveRDS(fitness_lineage,paste0(data_folder,'/fitness_lineage_', Sys.Date(), '.RDS'))

fitness_mutation=read.delim(mutation_result, header=T, stringsAsFactors = F)
fitness_mutation=fitness_mutation[,c(1,2,8)]
names(fitness_mutation)=c('rank','mutation','fitness score')
saveRDS(fitness_mutation,paste0(data_folder,'/fitness_mutation_', Sys.Date(), '.RDS'))



#date_voc<-c('2022-06-17','2022-06-24')
date_voc=readRDS(paste0(data_folder,'date_voc.rds'))
date_voc=c(date_voc, as.character(this_run))
date_voc=sort(unique(date_voc), decreasing = T)
saveRDS(date_voc, paste0(data_folder,'date_voc.rds'))

#date_who<-c('2022-06-15', '2022-06-22')
date_who=readRDS(paste0(data_folder,'date_who.rds'))
date_who=c(date_who, who_date)
date_who=sort(unique(date_who), decreasing = T)
saveRDS(date_who, paste0(data_folder,'date_who.rds'))


# ref<-'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*'
# 
# mutation_list=prevalent_mutation_sub2_s[prevalent_mutation_sub2_s$prevalence>0.5,]
# 
# lineages=unique(prevalent_mutation_sub2_s$lineage)
# 
# seqs=NULL
# 
# for(l in lineages){
#   mutations=mutation_list[mutation_list$lineage==l,]
#   alt=ref
#   for(i in 1:nrow(mutations)){
#     start=mutations$codon_num[i]
#     if(mutations$type[i]=='substitution'){
#       stopifnot(mutations$ref_aa[i]==substr(ref, start, start) )
#       substr(alt, start, start)=mutations$alt_aa[i]
#     }else{
#       end=as.numeric(mutations$codon_end[i])
#       substr(alt, start, end)=paste0(rep('-',end-start+1),collapse = '')
#     }
#   }
#   seqs=c(seqs,alt)
#   cat(paste0('>',l,"\n", collapse = ''))
#   cat(alt)
#   cat("\n")
# }




# fitness_lineage=read.delim(fitness_result, header=T, stringsAsFactors = F)
# fitness_lineage_top=fitness_lineage[fitness_lineage$rank<=20, c(1,2,5)]
# 
# names(fitness_lineage_top)=c('rank','lineage','fitness score')
# 
# reactable(fitness_lineage_top)
# 
# fit_prevalence=getPrevalence(pangolin_lineage=fitness_lineage_top$lineage, logInfo = FALSE)
# 
# plot_ly(fit_prevalence, x = ~date, y = ~proportion, color = ~lineage) %>%
#   add_lines() %>% layout(title = 'Daily prevalance of the top 20 lineages based on predicted fitness (some may not show up)')
# 
# fit_mutations = getMutationsByLineage(
#   pangolin_lineage=fitness_lineage_top$lineage, 
#   frequency=min_freq,
#   logInfo = FALSE)
# 
# 
# plot_ly(fit_mutations,
#         x=~reorder(mutation, codon_num), y=~lineage, z=~prevalence, 
#         type = "heatmap", colors = "Greys",
#         height= 25*length(unique(prevalent_mutation_sub2_s$lineage))+100,
#         width = 12*length(unique(prevalent_mutation_sub2_s$mutation))+100) %>%
#   layout(margin=list(l=50, r=50, t=50, b=50, pad = 4),showlegend=FALSE,
#          title = 'Mutation frequency of the top 20 lineages based on predicted fitness (some may not show up)')


#system('wget https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv -O COVID/data/lineages.csv')

#system('wget https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json -O COVID/data/alias_key.json')

#lineages=read.csv('COVID/data/lineages.csv')
#lineages$taxon=gsub(' ','_', lineages$taxon)

#prot<-read.csv('COVID/data/lineage_prot.csv')





ref_nuc='ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACATTCAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTATACATGTCTCTGGGACCAATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGTCTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAACGCTACTAATGTTGTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTGGGTGTTTATTACCACAAAAACAACAAAAGTTGGATGGAAAGTGAGTTCAGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATGTCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAATTTGTGTTTAAGAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATTAATTTAGTGCGTGATCTCCCTCAGGGTTTTTCGGCTTTAGAACCATTGGTAGATTTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTACTTGCTTTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGTGCTGCAGCTTATTATGTGGGTTATCTTCAACCTAGGACTTTTCTATTAAAATATAATGAAAATGGAACCATTACAGATGCTGTAGACTGTGCACTTGACCCTCTCTCAGAAACAAAGTGTACGTTGAAATCCTTCACTGTAGAAAAAGGAATCTATCAAACTTCTAACTTTAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAATTACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTAGTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCTACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGGCACAGGTGTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGACATTGCTGACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGTTCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTTCTTTATCAGGATGTTAACTGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACTCCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTAATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGCAGGTATATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAACTCTATTGCCATACCCACAAATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCTATGACCAAGACATCAGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGAATGCAGCAATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCTTTAACTGGAATAGCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAATTTACAAAACACCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGATCCATCAAAACCAAGCAAGAGGTCATTTATTGAAGATCTACTTTTCAACAAAGTGACACTTGCAGATGCTGGCTTCATCAAACAATATGGTGATTGCCTTGGTGATATTGCTGCTAGAGACCTCATTTGTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGAAATGATTGCTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAGGTGCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAGTTACACAGAATGTTCTCTATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTATTGGCAAAATTCAAGACTCACTTTCTTCCACAGCAAGTGCACTTGGAAAACTTCAAGATGTGGTCAACCAAAATGCACAAGCTTTAAACACGCTTGTTAAACAACTTAGCTCCAATTTTGGTGCAATTTCAAGTGTTTTAAATGATATCCTTTCACGTCTTGACAAAGTTGAGGCTGAAGTGCAAATTGATAGGTTGATCACAGGCAGACTTCAAAGTTTGCAGACATATGTGACTCAACAATTAATTAGAGCTGCAGAAATCAGAGCTTCTGCTAATCTTGCTGCTACTAAAATGTCAGAGTGTGTACTTGGACAATCAAAAAGAGTTGATTTTTGTGGAAAGGGCTATCATCTTATGTCCTTCCCTCAGTCAGCACCTCATGGTGTAGTCTTCTTGCATGTGACTTATGTCCCTGCACAAGAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACACTTTCCTCGTGAAGGTGTCTTTGTTTCAAATGGCACACACTGGTTTGTAACACAAAGGAATTTTTATGAACCACAAATCATTACTACAGACAACACATTTGTGTCTGGTAACTGTGATGTTGTAATAGGAATTGTCAACAACACAGTTTATGATCCTTTGCAACCTGAATTAGACTCATTCAAGGAGGAGTTAGATAAATATTTTAAGAATCATACATCACCAGATGTTGATTTAGGTGACATCTCTGGCATTAATGCTTCAGTTGTAAACATTCAAAAAGAAATTGACCGCCTCAATGAGGTTGCCAAGAATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTATGAGCAGTATATAAAATGGCCATGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGACAATTATGCTTTGCTGTATGACCAGTTGCTGTAGTTGTCTCAAGGGCTGTTGTTCTTGTGGATCCTGCTGCAAATTTGATGAAGACGACTCTGAGCCAGTGCTCAAAGGAGTCAAATTACATTACACATAA'

ref_prot='MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*'

lineage_max_date_nuc=read.csv('COVID/data/lineage_max_date_nuc.csv')
lineage_min_date_nuc=read.csv('COVID/data/lineage_min_date_nuc.csv')
lineage_freq_max_nuc=read.csv('COVID/data/lineage_freq_max_nuc.csv')
#lineage_max_date_prot=read.csv('data/lineage_max_date_prot.csv')
#lineage_min_date_prot=read.csv('data/lineage_min_date_prot.csv')
lineage_freq_max_prot=read.csv('COVID/data/lineage_freq_max_prot.csv')
lineages=unique(lineage_freq_max_prot$lineage)
lineages=c('Reference',lineages)

lineage_max_date_nuc = lineage_max_date_nuc %>% 
  rbind(c("Reference","2019-12-30","Wuhan/WIV04/2019", 
          "Spike|hCoV-19/Wuhan/WIV04/2019|2019-12-30|EPI_ISL_402124|Original|hCoV-19^^Hubei|Human|Wuhan Jinyintan Hospital|Wuhan Institute of Virology|Shi|China",
          ref_nuc, "Wuhan/WIV04/2019"))
lineage_min_date_nuc = lineage_min_date_nuc %>% 
  rbind(c("Reference","2019-12-30","Wuhan/WIV04/2019", 
          "Spike|hCoV-19/Wuhan/WIV04/2019|2019-12-30|EPI_ISL_402124|Original|hCoV-19^^Hubei|Human|Wuhan Jinyintan Hospital|Wuhan Institute of Virology|Shi|China",
          ref_nuc, "Wuhan/WIV04/2019"))
lineage_freq_max_nuc = lineage_freq_max_nuc %>% 
  rbind(c("Reference",1,ref_nuc))
lineage_freq_max_nuc$size=as.numeric(lineage_freq_max_nuc$size)
lineage_freq_max_prot = lineage_freq_max_prot %>% 
  rbind(c("Reference",1,ref_prot))
lineage_freq_max_prot$size=as.numeric(lineage_freq_max_prot$size)

saveRDS(ref_nuc,paste0(data_folder,'/ref_nuc.rds'))
saveRDS(ref_prot,paste0(data_folder,'/ref_prot.rds'))
saveRDS(lineages,paste0(data_folder,'/lineages.rds'))
saveRDS(lineage_max_date_nuc,paste0(data_folder,'/lineage_max_date_nuc.rds'))
saveRDS(lineage_min_date_nuc,paste0(data_folder,'/lineage_min_date_nuc.rds'))
saveRDS(lineage_freq_max_nuc,paste0(data_folder,'/lineage_freq_max_nuc.rds'))
saveRDS(lineage_freq_max_prot,paste0(data_folder,'/lineage_freq_max_prot.rds'))


library(seqinr)
lineage_prot=read.csv(paste0(data_folder, 'lineage_prot.csv'))
lineage_prot = lineage_prot %>% distinct(lineage, sequence, .keep_all=T)
write.fasta(sequences=as.list(lineage_prot$sequence),
            names=paste0(lineage_prot$lineage,':', lineage_prot$taxon),
            file.out=paste0(data_folder, "lineage_prot.fasta"))

lineage_nuc=read.csv(paste0(data_folder, 'lineage_nuc.csv'))
lineage_nuc = lineage_nuc %>% distinct(lineage, sequence, .keep_all=T)
write.fasta(sequences=as.list(lineage_nuc$sequence),
            names=paste0(lineage_nuc$lineage,':', lineage_nuc$taxon),
            file.out=paste0(data_folder, "lineage_nuc.fasta"))

system('makeblastdb -in COVID/data/lineage_prot.fasta -out COVID/data/lineage_prot -dbtype prot')
system('makeblastdb -in COVID/data/lineage_nuc.fasta -out COVID/data/lineage_nuc -dbtype nucl')


save.image(paste0(data_folder, '/all_variable_', Sys.Date(), '.RData'))

# scp -i ~/rvac.pem -r covid_variant_tracker_predictor/COVID/ ubuntu@44.208.199.237:~/

