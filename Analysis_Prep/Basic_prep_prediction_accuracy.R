library(tidyr)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)


# variance quantiles for ascertainment in Europe
quant=c(0,0.75,0.95,0.96,0.97,0.98,0.99)

# ascertainment column names
ascertainment_cols = c(
  "popA_top_100percent_var" ,
  "popA_maf.01", "popA_maf.001" ,"popA_top_25percent_var","popA_top_5percent_var", "popA_top_4percent_var", "popA_top_3percent_var", "popA_top_2percent_var","popA_top_1percent_var")

output_path="path/Prepped_Data/prediction_accuracy/"


#### BASIC ####

input_path="path/Basic/divergence_output/mut_track/"

files=list.files(path=input_path)

# get the files without referring to pop A or pop B
files_general = unique(sapply(files,function(f) gsub("_popA.txt|_popB.txt","",f) ))

# make one big data set containing every w and run for the population split times
Basic=do.call(rbind,lapply(files_general,function(f){

  data_popA = read.csv(paste0(input_path,f,"_popA.txt")) 
  data_popB = read.csv(paste0(input_path,f,"_popB.txt")) 
  
  split_df=full_join(data_popA,data_popB,by=setdiff(colnames(data_popA),"freq"),suffix=c("_A","_B")) %>%
    mutate(freq_A=replace_na(freq_A,0),
           freq_B=replace_na(freq_B,0))


  # add info about w and run, maybe add artificially added mutation data if the data is neutral
  w=strsplit(f,split="_")[[1]][1]
  w=as.numeric(gsub("w","",w))
  burn_run=strsplit(f,split="_")[[1]][3]
  burn_run=as.numeric(gsub("burnRun","",burn_run))
  div_run=strsplit(f,split="_")[[1]][4]
  div_run=as.numeric(gsub("divRun|.txt","",div_run))
  run=paste(burn_run,div_run,sep="-")
  corr=1
  nTraits=1
  sd_a=strsplit(f,split="_")[[1]][2]
  sd_a=as.numeric(gsub("sd","",sd_a))
  sim_type="Basic"
  opt_shift=0
  opt_shift_start_gen=NA


  all_df=cbind(data.frame(sim_type=sim_type,w=w,run=run,sd_a=sd_a,corr=corr,nTraits=nTraits, opt_shift=opt_shift, opt_shift_start_gen=opt_shift_start_gen),split_df)

  all_df = all_df %>% mutate(mut_type=ifelse(mut_type==1,"selected","neutral"))

  ### ENTER PREDICTION ACCURACY

  # rename IDs for uniqueness
  df = all_df %>% filter(mut_type=="selected") %>%
    mutate(mut_ID=paste0("w",w,"_r",run,"_m",mut_ID,"_p",position)) %>%
    mutate(mut_ID=gsub(" ","",mut_ID))

  # label ancestrally segregating SNPs (anc_seg = T or F)
  anc_seg_mutIDs = df$mut_ID[df$generation==1 & df$freq_A>0] # all ancestrally seg mutations in population of interest
  df = df %>% mutate(anc_seg = mut_ID %in% anc_seg_mutIDs)

  # get total variance of ancestrally segregating sites
  df_anc=df %>% filter(anc_seg & generation==1)
  anc_var = ( df_anc %>% summarize(Va=sum(effect_size^2*2*freq_A*(1-freq_A))))$Va

  # work generation by generation
  df = df %>% group_by(generation) %>% group_split()

  do.call("rbind", lapply(df, function(df_g){

    # make a list of mutations in each category that we care about
    ascertainment = list(
      popA_maf.01 = df_g$mut_ID[df_g$freq_A>=0.01 & df_g$freq_A<=0.99],
      popA_maf.001 = df_g$mut_ID[df_g$freq_A>=0.001 & df_g$freq_A<=0.999]
    )
    sub_df_g=df_g %>% filter(freq_A>0 & effect_size!=0) %>% mutate(var_popA = effect_size^2*2*freq_A*(1-freq_A))
    quant_limits=quantile(sub_df_g$var_popA,probs=quant) # get the var setting each quantile boundary specified above
    qu=lapply(quant_limits,function(q) sub_df_g$mut_ID[sub_df_g$var_popA>= q] )
    names(qu)=paste0("popA_top_",100*(1-quant), "percent_var")
    ascertainment=append(ascertainment,qu)

    get_numerator=function(mut_ID,freq,effect_size){
      data.frame(t(sapply(ascertainment,function(ascertained_mut){
        keep=mut_ID %in% ascertained_mut
        sum(2*effect_size[keep]^2*freq[keep]*(1-freq[keep]))
      })))

    }

    # in pop B contemporary with pop A
    space =df_g %>% group_by(sim_type,w,sd_a,corr,nTraits,opt_shift, opt_shift_start_gen ,run,generation) %>%
      group_modify(~ get_numerator(.x$mut_ID,.x$freq_B,.x$effect_size))

    space = space %>% pivot_longer(cols=all_of(ascertainment_cols),
                                   names_to="ascertainment",
                                   values_to="cov")
    var = sum(df_g$effect_size^2*2*df_g$freq_B*(1-df_g$freq_B))
    space$prediction_r2 = space$cov/var
    space$ascertainment_type="space"

    # in pop A where GWAS sample was drawn
    GWAS_pop =df_g %>% group_by(sim_type,w,sd_a,corr,nTraits,opt_shift, opt_shift_start_gen ,run,generation) %>%
      group_modify(~ get_numerator(.x$mut_ID,.x$freq_A,.x$effect_size))
    GWAS_pop = GWAS_pop %>% pivot_longer(cols=all_of(ascertainment_cols),
                                   names_to="ascertainment",
                                   values_to="cov")
    var_A = sum(df_g$effect_size^2*2*df_g$freq_A*(1-df_g$freq_A))
    GWAS_pop$prediction_r2 = GWAS_pop$cov/var_A
    GWAS_pop$ascertainment_type="GWAS_pop"

    # in ancestor to pop A
    time=data.frame(t(sapply(ascertainment,function(ascertained_mut){
      (df_anc %>% ungroup() %>% filter(mut_ID %in% ascertained_mut) %>%
         summarise(cov=sum(effect_size^2*2*freq_A*(1-freq_A))))$cov }))) %>%
      pivot_longer(cols=all_of(ascertainment_cols),
                   names_to="ascertainment",
                   values_to="cov")%>%
      mutate(prediction_r2=cov/anc_var,
             ascertainment_type="time")

    time=cbind(df_g %>% group_by(sim_type,w,sd_a,corr,nTraits,opt_shift, opt_shift_start_gen,run,generation) %>% group_keys(), time)

    bind_rows(bind_rows(space, GWAS_pop),time)

  }))




}))

saveRDS(Basic,paste0(output_path,"Basic.RDS"))
