# load libraries
lib <- c("magrittr", "tidyverse", 
         "data.table", "fastmatch",
         "beepr", "mailR", "kableExtra",
         "lme4", "janitor", "grid", "gridExtra")
lapply(lib, require, character.only = TRUE)
rm(lib)

# load homemade funs
source("homemade_funs.R")

#### local funs ####
adj_comb_complete <- function(Data) {
  A <- lapply(2:(length(Data)-1), sequence)
  B <- lapply(rev(vapply(A, length, 1L))-1, function(x) c(0, sequence(x)))
  c(unlist(lapply(seq_along(A), function(x) {
    lapply(B[[x]], function(y) Data[A[[x]]+y])
  }), recursive = FALSE, use.names = FALSE), list(Data))
}

seq2_5 <- function(df) {
  df %>%
    mutate(split = str_split(phon, "_") %>%
             sapply(adj_comb_complete)) %>%
    mutate(seq2 = split %>% 
             sapply(function(x) {
               x[lapply(x, length) == 2]
             }) %>%
             sapply(function(y) {
               sapply(y, function(j) {
                 paste(j, collapse = "_")
               })
             }),
           seq3 = split %>% 
             sapply(function(x) {
               x[lapply(x, length) == 3]
             }) %>%
             sapply(function(y) {
               sapply(y, function(j) {
                 paste(j, collapse = "_")
               })
             }),
           seq4 = split %>% 
             sapply(function(x) {
               x[lapply(x, length) == 4]
             }) %>%
             sapply(function(y) {
               sapply(y, function(j) {
                 paste(j, collapse = "_")
               })
             }),
           seq5 = split %>% 
             sapply(function(x) {
               x[lapply(x, length) == 5]
             }) %>%
             sapply(function(y) {
               sapply(y, function(j) {
                 paste(j, collapse = "_")
               })
             })
    )
}

extract_gram_types <- function(df, cat_type) {
  df %>%
    select(baby, section, types_root, cat) %>%
    apply(1, function(x) {
      tibble(id = x$baby, section = x$section, 
             word = unlist(x$types_root), cat = unlist(x$cat))
    }) %>% 
    rbindlist() %>%
    inner_join(., mot_phon, "word") %>%
    mutate(syllable = str_split(phon, "_") %>% 
             sapply(function(x) {
               x %in% vowels %>% 
                 sum()
             })) %>%
    filter(cat == cat_type)
}

assign_learned <- function(df) {
  df %<>%
    mutate(learned = unique(id) %>%
             sapply(function(x) {
               df %>%
                 filter(id == x) %>%
                 .$word %>%
                 unlist() %>%
                 {. %in% (chi_uni_cat %>%
                            filter(baby == x) %>%
                            .$types_root %>%
                            unlist())}
             }) %>%
             unlist())
}

restage_learned_words <- function(df) {
  chi_learned <- chi_uni_cat %>% 
    apply(1, function(x) {
      tibble(id = x$baby, section = x$section, 
             word = unlist(x$types_root))
    }) %>%
    rbindlist() 
  
  for (ids in unique(df$id)) {
    for (w in filter(df, id == ids)$word) {
      df$section[which(df$id == ids & df$word == w)] <-
        chi_learned %>%
        filter(id == ids, word == w) %>%
        .$section
    }
  }
  
  df
}

frame_analysis <- function(first_df, syll_num, cattype) {
  df1 <- extract_gram_types(mot_uni_cat, cattype)
  
  df2 <- df1 %>%
    filter(syllable == syll_num)
  
  df2 %<>%
    assign_learned()
  
  df3 <- df2 %>%
    filter(learned == TRUE)
  
  df3 %<>%
    restage_learned_words()
  
  df4 <- df2 %>%
    filter(learned == FALSE) %>%
    inner_join(., df3 %>%
                 group_by(id) %>%
                 summarize(avg_section = section %>%
                             mean() %>%
                             round(0)), "id") %>%
    select(-section) %>%
    rename(section = avg_section)
  
  df3 %<>%
    filter(str_detect(phon, "_")) %>%
    seq2_5() %>%
    mutate(seq2 = sapply(seq2, unique))
  
  df4 %<>%
    filter(str_detect(phon, "_")) %>%
    seq2_5() %>%
    mutate(seq2 = sapply(seq2, unique))
  
  list(df1, df2, df3, df4)
}

count_frames <- function(df, df_count, seq_num) {
  sapply(1:nrow(df), function(i) {
    df_count %>%
      filter(phon != df$phon[i],
             baby == df$id[i],
             section %in% 0:df$section[i]) %>%
      .[[seq_num]] %>%
      unlist() %>%
      str_count(paste("^", df[[seq_num]][[i]], "$", collapse = "|", sep="")) %>%
      sum()
  })
}

summary_frame_analysis <- function(df, df_count) {
  df %>%
    mutate(seq2_count = count_frames(df, df_count, "seq2"),
           seq3_count = count_frames(df, df_count, "seq3"),
           seq4_count = count_frames(df, df_count, "seq4"),
           seq5_count = count_frames(df, df_count, "seq5")) %>%
    group_by(id) %>%
    summarise(mean_seq2 = mean(seq2_count, na.rm = TRUE),
              mean_seq3 = mean(seq3_count, na.rm = TRUE),
              mean_seq4 = mean(seq4_count, na.rm = TRUE),
              mean_seq5 = mean(seq5_count, na.rm = TRUE)) %>%
    adorn_percentages() %>%
    bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(., na.rm = TRUE) else "MEAN"))) %>%
    bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(., na.rm = TRUE) else "SD")))
}

freq_up_stage <- function(df) {
  sapply(1:nrow(df), function(i) {
    mot_na_baby_section %>%
      filter(word == df$word[i],
             baby == df$id[i],
             section %in% 0:df$section[i]) %>%
      nrow()
  })
}

summary_freq <- function(df, name_var) {
  df %>%
    group_by(id) %>%
    summarize(!!name_var := mean(freq))
}

equalize_mean <- function(df, id_name, var_name) {
  df %>% 
    filter(id == id_name) %>%
    arrange(desc(freq)) %>%
    (function(x) {
      mean_equal = FALSE
      row = 1
      while (mean_equal == FALSE) {
        x_temp <- x[-c(1:row),]
        
        if (mean(x_temp$freq) < (summary_freq_sets[[var_name]][which(summary_freq_sets$id == id_name)] + 0)) {
          mean_equal = TRUE
        } else {
          row = row + 1
        }
      }
      if (row - 1 == 0) {
        x
      } else {
        x[-c(1:(row-sample(c(-2,3), 1))),]
      }
    })
}

match_df <- function(summary_freq_sets, temp) {
  matched_dfs <- list()
  for (rows in 1:nrow(temp)) {
    matched_dfs[[paste(temp$dfs[rows], "match", sep = "_")]] <-
      unique(get(temp$dfs[rows])$id) %>%
      lapply(function(x) {
        equalize_mean(get(temp$dfs[rows]), x, temp$vars[rows])
      }) %>%
      rbindlist()
  }
  
  matched_dfs
}

#### reference dfs ####
# mot_tokens
mot_tokens <- mot_na_baby_section %>%
  select(baby:word) %>%
  inner_join(., mot_phon, "word") %>%
  (function(df) {
    reference <- df[!duplicated(df$phon),] %>%
      mutate(syllable = str_split(phon, "_") %>% 
               sapply(function(x) {
                 x %in% vowels %>% 
                   sum()
               }),
             phonemes = str_split(phon, "_") %>% 
               sapply(function(x) {
                 x %>% 
                   length()
               })) %>%
      select(phon:phonemes)
    
    df %>%
      inner_join(., reference, "phon")
  }) %>%
  filter(str_detect(phon, "_")) %>%
  (function(y) {
    reference <- y[!duplicated(y$phon),] %>%
      mutate(split = str_split(phon, "_") %>%
               sapply(function(x) {
                 adj_comb_complete(x)
               })) %>%
      mutate(seq2 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 2]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq3 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 3]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq4 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 4]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq5 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 5]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               })
      ) 
    
    y %>%
      inner_join(., reference, "phon")
  }) %>%
  select(baby.x:phonemes.x, split:seq5) %>%
  (function(x) {
    colnames(x) <- str_remove(colnames(x), "[.]{1}.*$")
    x
  }) %>%
  mutate(seq2 = sapply(seq2, unique))

mot_types <- mot_tokens %>%
  arrange(baby, section) %>%
  group_by(baby) %>%
  distinct(phon, .keep_all = TRUE) 

chi_phon <- phon
colnames(chi_phon) <- c("word", "phon")
rm(phon)

chi_na_id_sec <- mot_chi_na %>%
  filter(id == "CHI") %>%
  arrange(baby, hour, half_hour) %>%
  mutate(section = c(rep(1:20, 1038)[order(rep(1:20, 1038))], rep(20, 11),
                     rep(1:20, 884)[order(rep(1:20, 884))], rep(20, 12),
                     rep(1:20, 1222)[order(rep(1:20, 1222))],
                     rep(1:20, 1293)[order(rep(1:20, 1293))], rep(20, 13),
                     rep(1:20, 1100)[order(rep(1:20, 1100))], rep(20, 16),
                     rep(1:20, 866)[order(rep(1:20, 866))], 20,
                     rep(1:20, 940)[order(rep(1:20, 940))], rep(20, 15),
                     rep(1:20, 685)[order(rep(1:20, 685))], rep(20, 11),
                     rep(1:20, 835)[order(rep(1:20, 835))], rep(20, 2),
                     rep(1:20, 902)[order(rep(1:20, 902))], 20,
                     rep(1:20, 1046)[order(rep(1:20, 1046))], rep(20, 8),
                     rep(1:20, 862)[order(rep(1:20, 862))], rep(20, 13))) %>%
  select(-id, -half_hour, -hour, -mor) %>%
  rename(id = baby) %>%
  select(id, section, string) %>%
  apply(., 1, function(x) {
    x %>%
      (function(x) {
        tibble(id = x$id, 
               section = x$section, 
               word = unlist(x$string))
      })
  }) %>%
  rbindlist() %>%
  inner_join(., chi_phon, "word") 

chi_tokens <- chi_na_id_sec %>%
  filter(str_detect(phon, "_")) %>%
  (function(y) {
    reference <- y[!duplicated(y$phon),] %>%
      mutate(split = str_split(phon, "_") %>%
               sapply(function(x) {
                 adj_comb_complete(x)
               })) %>%
      mutate(seq2 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 2]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq3 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 3]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq4 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 4]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               }),
             seq5 = split %>% 
               sapply(function(x) {
                 x[lapply(x, length) == 5]
               }) %>%
               sapply(function(y) {
                 sapply(y, function(j) {
                   paste(j, collapse = "_")
                 })
               })
      ) 
    
    y %>%
      inner_join(., reference, "phon")
  }) %>%
  select(id.x:phon, split:seq5) %>%
  (function(x) {
    colnames(x) <- str_remove(colnames(x), "[.]{1}.*$")
    x
  }) %>%
  mutate(seq2 = sapply(seq2, unique)) %>%
  rename(baby = id)

chi_types <- chi_tokens %>%
  arrange(baby, section) %>%
  group_by(baby) %>%
  distinct(phon, .keep_all = TRUE) 

#### frames analysis  ####
# (input as tokens or types)
# mono nouns
temp <- frame_analysis(mot_uni_cat, 1, "N")
mot_nouns <- temp[[1]]
mot_nouns_mono <- temp[[2]]
mot_nouns_mono_learned <- temp[[3]]
mot_nouns_mono_not_learned <- temp[[4]] ; rm(temp)

summary_seq_tokens_mono_learned <- summary_frame_analysis(mot_nouns_mono_learned, mot_tokens)
summary_seq_tokens_mono_not_learned <- summary_frame_analysis(mot_nouns_mono_not_learned, mot_tokens)

summary_seq_types_mono_learned <- summary_frame_analysis(mot_nouns_mono_learned, mot_types)
summary_seq_types_mono_not_learned <- summary_frame_analysis(mot_nouns_mono_not_learned, mot_types)

summary_seq_tokens_mono_learned_chi <- summary_frame_analysis(mot_nouns_mono_learned, chi_tokens)
summary_seq_tokens_mono_not_learned_chi <- summary_frame_analysis(mot_nouns_mono_not_learned, chi_tokens)

summary_seq_types_mono_learned_chi <- summary_frame_analysis(mot_nouns_mono_learned, chi_types)
summary_seq_types_mono_not_learned_chi <- summary_frame_analysis(mot_nouns_mono_not_learned, chi_types)

# bis nouns
temp <- frame_analysis(mot_uni_cat, 2, "N")
mot_nouns <- temp[[1]]
mot_nouns_bis <- temp[[2]]
mot_nouns_bis_learned <- temp[[3]]
mot_nouns_bis_not_learned <- temp[[4]] ; rm(temp)

summary_seq_tokens_learned <- summary_frame_analysis(mot_nouns_bis_learned, mot_tokens)
summary_seq_tokens_not_learned <- summary_frame_analysis(mot_nouns_bis_not_learned, mot_tokens)

summary_seq_types_learned <- summary_frame_analysis(mot_nouns_bis_learned, mot_types)
summary_seq_types_not_learned <- summary_frame_analysis(mot_nouns_bis_not_learned, mot_types)

summary_seq_tokens_learned_chi <- summary_frame_analysis(mot_nouns_bis_learned, chi_tokens)
summary_seq_tokens_not_learned_chi <- summary_frame_analysis(mot_nouns_bis_not_learned, chi_tokens)

summary_seq_types_learned_chi <- summary_frame_analysis(mot_nouns_bis_learned, chi_types)
summary_seq_types_not_learned_chi <- summary_frame_analysis(mot_nouns_bis_not_learned, chi_types)

# tri nouns
temp <- frame_analysis(mot_uni_cat, 3, "N")
mot_nouns <- temp[[1]]
mot_nouns_tri <- temp[[2]]
mot_nouns_tri_learned <- temp[[3]]
mot_nouns_tri_not_learned <- temp[[4]] ; rm(temp)

summary_seq_tokens_tri_learned <- summary_frame_analysis(mot_nouns_tri_learned, mot_tokens)
summary_seq_tokens_tri_not_learned <- summary_frame_analysis(mot_nouns_tri_not_learned, mot_tokens)

summary_seq_types_tri_learned <- summary_frame_analysis(mot_nouns_tri_learned, mot_types)
summary_seq_types_tri_not_learned <- summary_frame_analysis(mot_nouns_tri_not_learned, mot_types)

summary_seq_tokens_tri_learned_chi <- summary_frame_analysis(mot_nouns_tri_learned, chi_tokens)
summary_seq_tokens_tri_not_learned_chi <- summary_frame_analysis(mot_nouns_tri_not_learned, chi_tokens)

summary_seq_types_tri_learned_chi <- summary_frame_analysis(mot_nouns_tri_learned, chi_types)
summary_seq_types_tri_not_learned_chi <- summary_frame_analysis(mot_nouns_tri_not_learned, chi_types)

# mono verbs
temp <- frame_analysis(mot_uni_cat, 1, "V")
mot_verbs <- temp[[1]]
mot_verbs_mono <- temp[[2]]
mot_verbs_mono_learned <- temp[[3]]
mot_verbs_mono_not_learned <- temp[[4]] ; rm(temp)

summary_seq_tokens_mono_learned_verbs <- summary_frame_analysis(mot_verbs_mono_learned, mot_tokens)
summary_seq_tokens_mono_not_learned_verbs <- summary_frame_analysis(mot_verbs_mono_not_learned, mot_tokens)

summary_seq_types_mono_learned_verbs <- summary_frame_analysis(mot_verbs_mono_learned, mot_types)
summary_seq_types_mono_not_learned_verbs <- summary_frame_analysis(mot_verbs_mono_not_learned, mot_types)

summary_seq_tokens_mono_learned_chi_verbs <- summary_frame_analysis(mot_verbs_mono_learned, chi_tokens)
summary_seq_tokens_mono_not_learned_chi_verbs <- summary_frame_analysis(mot_verbs_mono_not_learned, chi_tokens)

summary_seq_types_mono_learned_chi_verbs <- summary_frame_analysis(mot_verbs_mono_learned, chi_types)
summary_seq_types_mono_not_learned_chi_verbs <- summary_frame_analysis(mot_verbs_mono_not_learned, chi_types)

# bis verbs
temp <- frame_analysis(mot_uni_cat, 2, "V")
mot_verbs <- temp[[1]]
mot_verbs_bis <- temp[[2]]
mot_verbs_bis_learned <- temp[[3]]
mot_verbs_bis_not_learned <- temp[[4]] ; rm(temp)

summary_seq_tokens_learned_verbs <- summary_frame_analysis(mot_verbs_bis_learned, mot_tokens)
summary_seq_tokens_not_learned_verbs <- summary_frame_analysis(mot_verbs_bis_not_learned, mot_tokens)

summary_seq_types_learned_verbs <- summary_frame_analysis(mot_verbs_bis_learned, mot_types)
summary_seq_types_not_learned_verbs <- summary_frame_analysis(mot_verbs_bis_not_learned, mot_types)

summary_seq_tokens_learned_chi_verbs <- summary_frame_analysis(mot_verbs_bis_learned, chi_tokens)
summary_seq_tokens_not_learned_chi_verbs <- summary_frame_analysis(mot_verbs_bis_not_learned, chi_tokens)

summary_seq_types_learned_chi_verbs <- summary_frame_analysis(mot_verbs_bis_learned, chi_types)
summary_seq_types_not_learned_chi_verbs <- summary_frame_analysis(mot_verbs_bis_not_learned, chi_types)

# tri verbs
temp <- frame_analysis(mot_uni_cat, 3, "V")
mot_verbs <- temp[[1]]
mot_verbs_tri <- temp[[2]]
mot_verbs_tri_learned <- temp[[3]]
mot_verbs_tri_not_learned <- temp[[4]] ; rm(temp)

summary_seq_tokens_tri_learned_verbs <- summary_frame_analysis(mot_verbs_tri_learned, mot_tokens)
summary_seq_tokens_tri_not_learned_verbs <- summary_frame_analysis(mot_verbs_tri_not_learned, mot_tokens)

summary_seq_types_tri_learned_verbs <- summary_frame_analysis(mot_verbs_tri_learned, mot_types)
summary_seq_types_tri_not_learned_verbs <- summary_frame_analysis(mot_verbs_tri_not_learned, mot_types)

summary_seq_tokens_tri_learned_chi_verbs <- summary_frame_analysis(mot_verbs_tri_learned, chi_tokens)
summary_seq_tokens_tri_not_learned_chi_verbs <- summary_frame_analysis(mot_verbs_tri_not_learned, chi_tokens)

summary_seq_types_tri_learned_chi_verbs <- summary_frame_analysis(mot_verbs_tri_learned, chi_types)
summary_seq_types_tri_not_learned_chi_verbs <- summary_frame_analysis(mot_verbs_tri_not_learned, chi_types)

#### frequency up to stage ####
# nouns
mot_nouns_mono_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_nouns_mono_not_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_nouns_bis_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_nouns_bis_not_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_nouns_tri_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_nouns_tri_not_learned %<>%
  mutate(freq = freq_up_stage(.))

# verbs
mot_verbs_mono_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_verbs_mono_not_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_verbs_bis_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_verbs_bis_not_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_verbs_tri_learned %<>%
  mutate(freq = freq_up_stage(.))

mot_verbs_tri_not_learned %<>%
  mutate(freq = freq_up_stage(.))

#### dfs matched on frequency ####
# nouns
summary_freq_sets <- summary_freq(mot_nouns_mono_learned, "Mono_learned_freq") %>%
  cbind(summary_freq(mot_nouns_mono_not_learned, "Mono_not_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_bis_learned, "Bis_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_bis_not_learned, "Bis_not_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_tri_learned, "Tri_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_tri_not_learned, "Tri_not_learned_freq")[,2]) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD")))

temp <- tibble(dfs = c("mot_nouns_mono_learned",
                       "mot_nouns_bis_learned",
                       "mot_nouns_tri_learned"),
               vars = colnames(summary_freq_sets)[-1] %>% .[grepl("not",.)])

matched_dfs <- match_df(summary_freq_sets, temp)

# verbs
summary_freq_sets <- summary_freq(mot_verbs_mono_learned, "Mono_learned_freq") %>%
  cbind(summary_freq(mot_verbs_mono_not_learned, "Mono_not_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_verbs_bis_learned, "Bis_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_verbs_bis_not_learned, "Bis_not_learned_freq")[,2]) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD")))

temp <- tibble(dfs = c("mot_verbs_mono_learned",
                       "mot_verbs_bis_learned"),
               vars = colnames(summary_freq_sets)[-1] %>% .[grepl("not",.)])

matched_dfs_verbs <- match_df(summary_freq_sets, temp)

rm(summary_freq_sets, temp)

summary_seq_tokens_mono_learned_match <- summary_frame_analysis(matched_dfs[["mot_nouns_mono_learned_match"]], mot_tokens) 
summary_seq_tokens_bis_learned_match <- summary_frame_analysis(matched_dfs[["mot_nouns_bis_learned_match"]], mot_tokens) 
summary_seq_tokens_tri_learned_match <- summary_frame_analysis(matched_dfs[["mot_nouns_tri_learned_match"]], mot_tokens) 

summary_seq_tokens_mono_learned_match_chi <- summary_frame_analysis(matched_dfs[["mot_nouns_mono_learned_match"]], chi_tokens) 
summary_seq_tokens_bis_learned_match_chi <- summary_frame_analysis(matched_dfs[["mot_nouns_bis_learned_match"]], chi_tokens) 
summary_seq_tokens_tri_learned_match_chi <- summary_frame_analysis(matched_dfs[["mot_nouns_tri_learned_match"]], chi_tokens) 

summary_seq_types_mono_learned_match <- summary_frame_analysis(matched_dfs[["mot_nouns_mono_learned_match"]], mot_types) 
summary_seq_types_bis_learned_match <- summary_frame_analysis(matched_dfs[["mot_nouns_bis_learned_match"]], mot_types) 
summary_seq_types_tri_learned_match <- summary_frame_analysis(matched_dfs[["mot_nouns_tri_learned_match"]], mot_types) 

summary_seq_types_mono_learned_match_chi <- summary_frame_analysis(matched_dfs[["mot_nouns_mono_learned_match"]], chi_types) 
summary_seq_types_bis_learned_match_chi <- summary_frame_analysis(matched_dfs[["mot_nouns_bis_learned_match"]], chi_types) 
summary_seq_types_tri_learned_match_chi <- summary_frame_analysis(matched_dfs[["mot_nouns_tri_learned_match"]], chi_types) 

summary_seq_tokens_mono_learned_match_verbs <- summary_frame_analysis(matched_dfs_verbs[["mot_verbs_mono_learned_match"]], mot_tokens) 
summary_seq_tokens_mono_learned_match_chi_verbs <- summary_frame_analysis(matched_dfs_verbs[["mot_verbs_mono_learned_match"]], chi_tokens) 
summary_seq_types_mono_learned_match_verbs <- summary_frame_analysis(matched_dfs_verbs[["mot_verbs_mono_learned_match"]], mot_types) 
summary_seq_types_mono_learned_match_chi_verbs <- summary_frame_analysis(matched_dfs_verbs[["mot_verbs_mono_learned_match"]], chi_types) 
