library(tidyr)
library(dplyr)
library(stringr)

# variance quantiles for ascertainment in Europe
quant=c(0,0.75,0.95,0.96,0.97,0.98,0.99)

# function to get mean frequency of allele in two populations with dplyr
get_mean=Vectorize(function(f1,f2){
  mean(c(f1,f2))
})

output_path="path/Prepped_Data/Qst_Fst/"

#### GxE ####

input_path="path/GxE/divergence_output/mut_track/"

files=list.files(path=input_path)

# get the files without referring to pop A or pop B
files_general = unique(sapply(files,function(f) gsub("_popA.txt|_popB.txt","",f) ))

# make one big data set containing every w and run for the population split times
GxE=do.call(rbind,lapply(files_general,function(f){

  data_popA=read.csv(paste0(input_path,f,"_popA.txt"))
  data_popB=read.csv(paste0(input_path,f,"_popB.txt"))

  split_df=full_join(data_popA,data_popB,by=setdiff(colnames(data_popA),"popA_freq")) %>%
    rename(freq_A=popA_freq, freq_B=popB_freq) %>%
    mutate(freq_A=replace_na(freq_A,0),
           freq_B=replace_na(freq_B,0))


  # add info about w and run, maybe add artificially added mutation data if the data is neutral
  w=strsplit(f,split="_")[[1]][1]
  w=as.numeric(gsub("w","",w))
  burn_run=strsplit(f,split="_")[[1]][4]
  burn_run=as.numeric(gsub("burnRun","",burn_run))
  div_run=strsplit(f,split="_")[[1]][5]
  div_run=as.numeric(gsub("divRun|.txt","",div_run))
  run=paste(burn_run,div_run,sep="-")
  corr=strsplit(f,split="_")[[1]][3]
  corr=as.numeric(gsub("corr","",corr))
  nTraits=1
  sd_a=strsplit(f,split="_")[[1]][2]
  sd_a=as.numeric(gsub("sd","",sd_a))
  sim_type="GxE"
  opt_shift=0
  opt_shift_start_gen=NA

  all_df=cbind(data.frame(sim_type=sim_type,w=w,run=run,sd_a=sd_a,corr=corr,nTraits=nTraits,opt_shift=opt_shift,opt_shift_start_gen=opt_shift_start_gen),split_df)

  all_df = all_df %>% mutate(mut_type=ifelse(mut_type==1,"selected","neutral"))

  ### ENTER QST FST

  # rename IDs for uniqueness
  d = all_df %>%
    mutate(mut_ID=paste0("w",w,"_r",run,"_m",mut_ID,"_p",position)) %>%
    mutate(mut_ID=gsub(" ","",mut_ID))

  # for now, use single Fst value btwn pops for a generation

  # get total population values
  # get the "total" population frequencies estimated from present day populations (not true ancestral)
  d = d %>% mutate(tot_freq=get_mean(freq_B,freq_A))

  # analyze generation by generation
  d = d %>% group_by(generation) %>% group_split()

  do.call("rbind",lapply(d,function(df){
    # get list of mutations to ascertain for Qst
    df_sel = df %>% filter(mut_type=="selected")
    # add ascertainment by variance contributed to ascertainment list
    sub_df_sel=df_sel %>% filter(freq_A>0 ) %>% mutate(var_popA = popA_effect_size^2*2*freq_A*(1-freq_A))
    quant_limits=quantile(sub_df_sel$var_popA,probs=quant) # get the var setting each quantile boundary specified above
    qu=lapply(quant_limits,function(q) sub_df_sel$mut_ID[sub_df_sel$var_popA>= q] )
    names(qu)=paste0("popA_top_",100*(1-quant), "percent_var")
    ascertainment=qu #append(ascertainment,qu)
    rm(df_sel)

    d=2 #number of demes
    df_neut = df %>%filter(mut_type=="neutral")
    fst_sigma_b = sum((1/(d-1))*((df_neut$freq_A-df_neut$tot_freq)^2 + (df_neut$freq_B-df_neut$tot_freq)^2))
    fst_sigma_denom = sum(df_neut$tot_freq*(1-df_neut$tot_freq))
    Fst = fst_sigma_b / fst_sigma_denom

    # for each ascertainment category, get Qst btwn Eur & EA
    get_Qst=function(mut_ID,freq_A,freq_B,tot_freq,effect_size_A,effect_size_B){
      tmp=data.frame(t(sapply(names(ascertainment),function(asc_name){

        ascertained_mut=ascertainment[[asc_name]]

        # identify indices of mutations belonging to ascertainment category
        keep=mut_ID %in% ascertained_mut

        if(asc_name=="all"){
          # make new subsetted vectors that only correspond to those mutations
          a_A = effect_size_A[keep]
          a_B = effect_size_B[keep]
          A_f = freq_A[keep]
          B_f = freq_B[keep]
          tot_f = tot_freq[keep]


          mean_A = sum(2*a_B*B_f)
          mean_B = sum(2*a_A*A_f)
          tot_var = sum(a_A^2*A_f + a_B^2*B_f - 2*a_A*a_B*A_f*B_f)

          c(mean_A=mean_A,mean_B=mean_B,tot_var=tot_var)
        } else{
          # make new subsetted vectors that only correspond to those mutations
          a = effect_size_A[keep]
          A_f = freq_A[keep]
          B_f = freq_B[keep]
          tot_f = tot_freq[keep]


          mean_A = sum(2*a*B_f)
          mean_B = sum(2*a*A_f)
          tot_var = sum(2*a^2*tot_f*(1-tot_f))

          c(mean_A=mean_A,mean_B=mean_B,tot_var=tot_var)
        }


      })))
      tmp$ascertainment = rownames(tmp)
      rownames(tmp)=NULL
      return(tmp)

    }

    # return data frame where for each ascertainment category (column) the row is the value of Qx
    Q= df %>% group_by(sim_type,w,sd_a,corr,nTraits,opt_shift,opt_shift_start_gen,run,generation) %>%
      group_modify(~ get_Qst(.x$mut_ID,.x$freq_A,.x$freq_B,.x$tot_freq,.x$popA_effect_size, .x$popB_effect_size)) %>%
      mutate(Qst = ((mean_A- get_mean(mean_A,mean_B))^2 + (mean_B-get_mean(mean_A,mean_B))^2)/(2*tot_var),
             Fst = Fst,
             Qx = Qst / Fst,
             sim_type = paste0(sim_type,"_wrongEffects"))



    return(Q)



  }))




}))

saveRDS(GxE,paste0(output_path,"GxE.RDS"))
