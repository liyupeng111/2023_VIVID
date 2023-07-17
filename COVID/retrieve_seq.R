library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(jsonlite)
library(plotly)
library(ggplot2)
library(readr)
library(reactable)
library(dplyr)
library(Biostrings)
library(seqinr)

ref_nuc=readRDS('data/ref_nuc.rds')
ref_prot=readRDS('data/ref_prot.rds')
lineages=readRDS('data/lineages.rds')
lineage_max_date_nuc=readRDS('data/lineage_max_date_nuc.rds')
lineage_min_date_nuc=readRDS('data/lineage_min_date_nuc.rds')
lineage_freq_max_nuc=readRDS('data/lineage_freq_max_nuc.rds')
lineage_freq_max_prot=readRDS('data/lineage_freq_max_prot.rds')

get_descendent = function(lineage, descendent="No"){
  lineage_df=NULL
  for(l in lineage){
    if(l %in% c('A','B','Reference') | descendent=='No'){
      tmp=data.frame(lineage=l, descendent=l)
    }else{
      tmp=data.frame(lineage=l, descendent=lineages[grepl(paste0('^',l),lineages)])
    }
    lineage_df=rbind(lineage_df, tmp)
  }
 return(lineage_df)
}


retrieve_sequence = function(lineage, descendent, generation_method, seq_type){

  if(is.null(lineage)){return('')}
  lineage_df=get_descendent(lineage, descendent)
  
  if(generation_method=='Earliest'){
    
    seq_df=merge(lineage_min_date_nuc, lineage_df, by.x='lineage', by.y='descendent', all=F)
    lin_df=aggregate(seq_df$date, by=list(seq_df$lineage.y), FUN='min')
    seq_df=merge(seq_df,lin_df, by.x=c('lineage.y','date'), by.y=c('Group.1','x'), all=F)
    seq_df_nuc = seq_df %>% distinct(lineage.y, sequence, .keep_all=T)
    seq_df_rna = seq_df_nuc
    seq_df_rna$sequence=gsub('T','U', seq_df_rna$sequence)
    seq_df_prot = seq_df_nuc
    seq_df_prot$sequence=unlist(lapply(seq_df_nuc$sequence, 
                                       FUN=function(x){
                                         as.character(Biostrings::translate(DNAString(x)))
                                         }
                                       )
                                )
    seq_df_rna$num <- ave(seq_df_rna$sequence, seq_df_rna$lineage, FUN = seq_along)
    seq_df_nuc$num <- ave(seq_df_nuc$sequence, seq_df_nuc$lineage, FUN = seq_along)
    seq_df_prot$num <- ave(seq_df_prot$sequence, seq_df_prot$lineage, FUN = seq_along)
    
  }else if(generation_method=='Latest'){
  
    seq_df=merge(lineage_max_date_nuc, lineage_df, by.x='lineage', by.y='descendent', all=F)
    lin_df=aggregate(seq_df$date, by=list(seq_df$lineage.y), FUN='max')
    seq_df=merge(seq_df,lin_df, by.x=c('lineage.y','date'), by.y=c('Group.1','x'), all=F)
    seq_df_nuc = seq_df %>% distinct(lineage.y, sequence, .keep_all=T)
    seq_df_rna = seq_df_nuc
    seq_df_rna$sequence=gsub('T','U', seq_df_rna$sequence)
    seq_df_prot = seq_df_nuc
    seq_df_prot$sequence=unlist(lapply(seq_df_nuc$sequence, FUN=function(x){
      as.character(Biostrings::translate(DNAString(x)))
      }
    ))
    seq_df_rna$num <- ave(seq_df_rna$sequence, seq_df_rna$lineage, FUN = seq_along)
    seq_df_nuc$num <- ave(seq_df_nuc$sequence, seq_df_nuc$lineage, FUN = seq_along)
    seq_df_prot$num <- ave(seq_df_prot$sequence, seq_df_prot$lineage, FUN = seq_along)
    
  }else if(generation_method=='Most frequent DNA'){
    
    seq_df=merge(lineage_freq_max_nuc, lineage_df, by.x='lineage', by.y='descendent', all=F)
    lin_df_sum=aggregate(seq_df$size, by=list(seq_df$lineage.y, seq_df$sequence), FUN='sum')
    lin_df_max=aggregate(lin_df_sum$x, by=list(lin_df_sum$Group.1), FUN='max')
    seq_df=merge(lin_df_max, lin_df_sum, all=F)
    seq_df_nuc = seq_df %>% distinct(Group.1, Group.2, .keep_all=T)
    names(seq_df_nuc)=c('lineage','freq','sequence')
    seq_df_rna = seq_df_nuc
    seq_df_rna$sequence=gsub('T','U', seq_df_rna$sequence)
    seq_df_prot = seq_df_nuc
    seq_df_prot$sequence=unlist(lapply(seq_df_nuc$sequence, FUN=function(x){
      as.character(Biostrings::translate(DNAString(x)))
    }
    ))
    seq_df_rna$num <- ave(seq_df_rna$sequence, seq_df_rna$lineage, FUN = seq_along)
    seq_df_nuc$num <- ave(seq_df_nuc$sequence, seq_df_nuc$lineage, FUN = seq_along)
    seq_df_prot$num <- ave(seq_df_prot$sequence, seq_df_prot$lineage, FUN = seq_along)
    
  }else if(generation_method=='Most frequent protein'){
    
    seq_df=merge(lineage_freq_max_prot, lineage_df, by.x='lineage', by.y='descendent', all=F)
    lin_df_sum=aggregate(seq_df$size, by=list(seq_df$lineage.y, seq_df$sequence), FUN='sum')
    lin_df_max=aggregate(lin_df_sum$x, by=list(lin_df_sum$Group.1), FUN='max')
    seq_df=merge(lin_df_max, lin_df_sum, all=F)
    seq_df_prot = seq_df %>% distinct(Group.1, Group.2, .keep_all=T)
    names(seq_df_prot)=c('lineage','freq','sequence')
    seq_df_prot$num <- ave(seq_df_prot$sequence, seq_df_prot$lineage, FUN = seq_along)
    
  }
  
  final_string='<span style="background-color:#FFFFFF;padding:20px;width:70%; word-wrap:break-word; display:inline-block;font-family: monospace;">'
  
  if(generation_method=='Most frequent protein'){
    final_string=paste0(final_string, '<h3>Spike gene protein sequences</h3>')
    for(i in 1:nrow(seq_df_prot)){
      final_string=paste0(final_string, paste0('>', seq_df_prot$lineage[i],'',
                                          '_Spike_protein frequency:', 
                                          seq_df_prot$freq[i], '<br>',
                                          seq_df_prot$sequence[i],'<br>'))
    }
  }else {
    if('cDNA' %in% seq_type){
      final_string=paste0(final_string, paste0('<h3>Spike gene cDNA sequences</h3>'))
      if (generation_method=='Most frequent DNA'){
        for(i in 1:nrow(seq_df_nuc)){
          final_string=paste0(final_string, paste0('>', seq_df_nuc$lineage[i],'_',seq_df_nuc$num[i],
                                              '_Spike_cDNA' ,' frequency:', 
                                              seq_df_nuc$freq[i], '<br>',
                                              seq_df_nuc$sequence[i],'<br>'))
        }
      }else{
        for(i in 1:nrow(seq_df_nuc)){
          final_string=paste0(final_string, paste0('>', seq_df_nuc$lineage.y[i],'_',seq_df_nuc$num[i],
                                              '_Spike_cDNA',' date:', 
                                              seq_df_nuc$date[i],
                                              ' taxon:', seq_df_nuc$taxon[i],
                                              '<br>',
                                              seq_df_nuc$sequence[i],'<br>'))
        }
      }
    }
    
    if('mRNA' %in% seq_type){
      final_string=paste0(final_string, paste0('<h3>Spike gene mRNA sequences</h3>'))
      if (generation_method=='Most frequent DNA'){
        for(i in 1:nrow(seq_df_rna)){
          final_string=paste0(final_string, paste0('>', seq_df_rna$lineage[i],'_',seq_df_rna$num[i],
                                              '_Spike_mRNA' ,' frequency:', 
                                              seq_df_rna$freq[i], '<br>',
                                              seq_df_rna$sequence[i],'<br>'))
        }
      }else{
        for(i in 1:nrow(seq_df_rna)){
          final_string=paste0(final_string, paste0('>', seq_df_rna$lineage.y[i],'_',seq_df_rna$num[i],
                                              '_Spike_mRNA',' date:', 
                                              seq_df_rna$date[i],
                                              ' taxon:', seq_df_rna$taxon[i],
                                              '<br>',
                                              seq_df_rna$sequence[i],'<br>'))
        }
      }
    }
    
    if('Protein' %in% seq_type){
      final_string=paste0(final_string, paste0('<h3>Spike gene protein sequences</h3>'))
      if (generation_method=='Most frequent DNA'){
        for(i in 1:nrow(seq_df_prot)){
          final_string=paste0(final_string, paste0('>', seq_df_prot$lineage[i],'_',seq_df_prot$num[i],
                                              '_Spike_protein' ,' frequency:', 
                                              seq_df_prot$freq[i], '<br>',
                                              seq_df_prot$sequence[i],'<br>'))
        }
      }else{
        for(i in 1:nrow(seq_df_prot)){
          final_string=paste0(final_string, paste0('>', seq_df_prot$lineage.y[i],'_',seq_df_prot$num[i],
                                              '_Spike_protein',' date:', 
                                              seq_df_prot$date[i],
                                              ' taxon:', seq_df_prot$taxon[i],
                                              '<br>',
                                              seq_df_prot$sequence[i],'<br>'))
        }
      }
    }
  }
  
  final_string=paste0(final_string,'</span>')
  
  return(HTML(final_string))
}