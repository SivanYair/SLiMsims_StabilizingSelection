library(tidyr)
library(dplyr)
library(stringr)

# variance quantiles for ascertainment in pop A
quant=c(0,0.75,0.95,0.96,0.97,0.98,0.99)

# ascertainment column names
ascertainment_cols = c("all",
  "popA_top_100percent_var" ,
  "popA_maf.01", "popA_maf.001" ,"popA_top_25percent_var","popA_top_5percent_var", "popA_top_4percent_var", "popA_top_3percent_var", "popA_top_2percent_var","popA_top_1percent_var")

# function to get mean frequency of allele in two populations with dplyr
get_mean=Vectorize(function(f1,f2){
  mean(c(f1,f2))
})

output_path="path/Prepped_Data/mean_scores/"

#### GxE ####

input_path="path/GxE/divergence_output/mut_track/"

files=list.files(path=input_path)
files=setdiff(files,c("mut_track.tar.gz","stay.txt"))

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
    sub_df=df %>% filter(freq_A>0 ) %>% mutate(var_popA = popA_effect_size^2*2*freq_A*(1-freq_A))
    quant_limits=quantile(sub_df$var_popA,probs=quant) # get the var setting each quantile boundary specified above
    qu=lapply(quant_limits,function(q) sub_df$mut_ID[sub_df$var_popA>= q] )
    names(qu)=paste0("popA_top_",100*(1-quant), "percent_var")
    ascertainment=append(ascertainment,qu)


    # for each ascertainment category, get mean scores of pop A and pop B
    get_score_A=function(mut_ID,freq,effect_size){
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
      group_modify(~ get_score_A(.x$mut_ID,.x$freq_A,.x$popA_effect_size)) %>%
      pivot_longer(cols=all_of(ascertainment_cols),names_to="ascertainment",values_to="mean_A")

    get_score_B=function(mut_ID,freq,effect_size_A,effect_size_B){
      data.frame(t(sapply(names(ascertainment),function(ascertainment_cat){

        ascertained_mut=ascertainment[[ascertainment_cat]]

        # identify indices of mutations belonging to ascertainment category
        keep=mut_ID %in% ascertained_mut

        # if you ascertain all SNPs
        if(ascertainment_cat=="all"){
          a = effect_size_B[keep]
          f = freq[keep]
          return(sum(2*a*f))
        } else{
          a = effect_size_A[keep]
          f = freq[keep]
          return(sum(2*a*f))
        }

      })))

    }


    mean_B = df %>% group_by(sim_type,w,sd_a,corr,nTraits,opt_shift,opt_shift_start_gen,run,generation) %>%
      group_modify(~ get_score_B(.x$mut_ID,.x$freq_B,.x$popA_effect_size,.x$popB_effect_size)) %>%
      pivot_longer(cols=all_of(ascertainment_cols),names_to="ascertainment",values_to="mean_B")

    # bind columns together onto mean_A, the df we return
    mean_A$mean_B = mean_B$mean_B

    return(mean_A)
  }))

}))

saveRDS(GxE,paste0(output_path,"GxE.RDS"))
