library(tidyr)
library(dplyr)
library(stringr)
options(stringsAsFactors = F)


# variance quantiles for ascertainment in pop A
quant=c(0,0.75,0.95,0.96,0.97,0.98,0.99)

# ascertainment column names
ascertainment_cols = c("all",
  "popA_top_100percent_var" ,
  "popA_maf.01", "popA_maf.001" ,"popA_top_25percent_var","popA_top_5percent_var", "popA_top_4percent_var", "popA_top_3percent_var", "popA_top_2percent_var","popA_top_1percent_var")

output_path="output/Prepped_Data/mean_scores/"

#### DirSel ####

input_path="input/DirSel/divergence_output/mut_track/"

files=list.files(path=input_path)

# get the files without referring to pop A or pop B
files_general = unique(sapply(files,function(f) gsub("_popA.txt|_popB.txt","",f) ))

# make one big data set containing every w and run for the population split times
DirSel=do.call(rbind,lapply(files_general,function(f){

  data_popA = read.csv(paste0(input_path,f,"_popA.txt")) %>% distinct(mut_ID,generation,.keep_all=T)
  data_popB = read.csv(paste0(input_path,f,"_popB.txt")) %>% distinct(mut_ID,generation,.keep_all=T)

  split_df=full_join(data_popA,data_popB,by=setdiff(colnames(data_popA),"freq"),suffix=c("_A","_B")) %>%
    mutate(freq_A=replace_na(freq_A,0),
           freq_B=replace_na(freq_B,0))


  # add info about w and run, maybe add artificially added mutation data if the data is neutral
  w=strsplit(f,split="_")[[1]][1]
  w=as.numeric(gsub("w","",w))
  burn_run=strsplit(f,split="_")[[1]][5]
  burn_run=as.numeric(gsub("burnRun","",burn_run))
  div_run=strsplit(f,split="_")[[1]][6]
  div_run=as.numeric(gsub("divRun|.txt","",div_run))
  run=paste(burn_run,div_run,sep="-")
  corr=1
  nTraits=1
  sd_a=strsplit(f,split="_")[[1]][2]
  sd_a=as.numeric(gsub("sd","",sd_a))
  sim_type="DirSel"
  opt_shift=strsplit(f,split="_")[[1]][3]
  opt_shift=as.numeric(gsub("sdShift","",opt_shift))
  opt_shift_start_gen=strsplit(f,split="_")[[1]][4]
  opt_shift_start_gen=as.numeric(gsub("startGen","",opt_shift_start_gen))


  all_df=cbind(data.frame(sim_type=sim_type,w=w,run=run,sd_a=sd_a,corr=corr,nTraits=nTraits, opt_shift=opt_shift, opt_shift_start_gen=opt_shift_start_gen),split_df)

  # filter to have only QTL sites
  all_df = all_df %>% filter(mut_type==1)

  # rename IDs for uniqueness
  d = all_df %>%
    mutate(mut_ID=paste0("w",w,"_r",run,"_m",mut_ID,"_p",position)) %>%
    mutate(mut_ID=gsub(" ","",mut_ID))

  # analyze generation by generation
  d = d %>% group_by(generation) %>% group_split()

  do.call("rbind",lapply(d,function(df){

    # make a list of mutations in each category that we care about
    ascertainment = list(
      all = df$mut_ID,
      popA_maf.01 = df$mut_ID[df$freq_A>=0.01 & df$freq_A<=0.99],
      popA_maf.001 = df$mut_ID[df$freq_A>=0.001 & df$freq_A<=0.999]

    )
    # add ascertainment by variance contributed to ascertainment list
    sub_df=df %>% filter(freq_A>0 ) %>% mutate(var_popA = effect_size^2*2*freq_A*(1-freq_A))
    quant_limits=quantile(sub_df$var_popA,probs=quant) # get the var setting each quantile boundary specified above
    qu=lapply(quant_limits,function(q) sub_df$mut_ID[sub_df$var_popA>= q] )
    names(qu)=paste0("popA_top_",100*(1-quant), "percent_var")
    ascertainment=append(ascertainment,qu)


    # for each ascertainment category, get mean score
    get_score=function(mut_ID,freq,effect_size){
      data.frame(t(sapply(ascertainment,function(ascertained_mut){
        # identify indices of mutations belonging to ascertainment category
        keep=mut_ID %in% ascertained_mut

        # make new subsetted vectors that only correspond to those mutations
        a = effect_size[keep]
        f = freq[keep]


        sum(2*a*f)

      })))

    }

    # return data frame where for each ascertainment category (column) the row is the mean score
    mean_A = df %>% group_by(sim_type,w,sd_a,corr,nTraits,opt_shift,opt_shift_start_gen,run,generation) %>%
      group_modify(~ get_score(.x$mut_ID,.x$freq_A,.x$effect_size)) %>%
      pivot_longer(cols=all_of(ascertainment_cols),names_to="ascertainment",values_to="mean_A")


    mean_B = df %>% group_by(sim_type,w,sd_a,corr,nTraits,opt_shift,opt_shift_start_gen,run,generation) %>%
      group_modify(~ get_score(.x$mut_ID,.x$freq_B,.x$effect_size)) %>%
      pivot_longer(cols=all_of(ascertainment_cols),names_to="ascertainment",values_to="mean_B")

    # bind columns together onto mean_A, the df we return
    mean_A$mean_B = mean_B$mean_B

    return(mean_A)

  }))



}))

saveRDS(DirSel,paste0(output_path,"DirSel.RDS"))
