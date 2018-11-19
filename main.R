# load libraries
lib <- c("magrittr", "tidyverse", 
         "data.table", "fastmatch",
         "beepr", "mailR", "kableExtra",
         "lme4", "janitor", "grid", "gridExtra")
lapply(lib, require, character.only = TRUE)
rm(lib)

# load homemade funs
source("homemade_funs.R")

#### MOTHERS ####
# create mot_uni_cat
mot_uni_cat <- mot_na %>%
  (function(x) {
    for (i in 1:nrow(x)) {
      x$section2[[i]] <- rep(x$section[i], length(x$string[[i]]))
      x$baby2[[i]] <- rep(x$baby[i], length(x$string[[i]]))
    }
    x
  }) %>%
  select(-baby, -section) %>%
  rename(section = section2,
         baby = baby2) %>%
  (function(x) {
    tibble(baby = unlist(x$baby),
           section = unlist(x$section),
           word = unlist(x$string),
           mor = unlist(x$mor))
  }) %>%
  na.omit() %>%
  group_by(baby, section) %>%
  summarise(uni_diff = list(unique(word))) %>%
  ungroup() %>%
  (function(x) {
    x$uni_diff <- setdif_list(x, "uni_diff")
    x
  })

# clean mot_mor_uni from motVSspoken_bnc_study
lab_var <- function(x, new_lab, old_lab, new_lab2 = NULL) {
  if (length(new_lab2) == 0) {
    ifelse(x == old_lab, new_lab, x)
  } else {
    ifelse(x == old_lab, list(c(new_lab, new_lab2)), x)
  }
}

mot_mor_uni %<>%
  mutate(cat = cat %>%
           lab_var("N", "N:PROP") %>%
           lab_var("N", "N|-N-CL") %>%
           lab_var("N", "N-CL") %>%
           lab_var("N", "N:LET") %>%
           lab_var("N", "N:PROP|-N-CL") %>%
           lab_var("PRO", "PRO:INDEF") %>%
           lab_var("V", "V:AUX") %>%
           lab_var("N", "N:PROP|-N-CL|V", "V") %>%
           lab_var("N", "N|-N-CL|V", "V") %>%
           lab_var("ADV", "WH:ADV") %>%
           lab_var("PRO", "PRO:DEM") %>%
           lab_var("PRO", "PRO:REFL") %>%
           lab_var("PRO", "PRO:DEM|-N-CL|V", "V") %>%
           lab_var("PRO", "PRO:POSS") %>%
           lab_var("PRO", "PRO|-N-CL") %>%
           lab_var("V", "V|-V-CL|NEG") %>%
           lab_var("PRO", "WH:PRO") %>%
           lab_var("ADV", "ADV:INT") %>%
           lab_var("ADV", "ADV|-N-CL|V", "V") %>%
           lab_var("N", "N|-PL-N-CL") %>%
           lab_var("V", "INF") %>%
           lab_var("PRO", "PRO:INDEF|-N-CL|V", "V") %>%
           lab_var("PRO", "PRO|-N-CL|V", "V") %>%
           lab_var("V", "V|-N-CL|V") %>%
           lab_var("PRO", "WH:PRO|-N-CL|V", "V") %>%
           lab_var("ADV", "ADV|-N-CL") %>%
           lab_var("N", "N:PREP|-N-CL") %>%
           lab_var("N", "N:PROP|-PL-N-CL") %>%
           lab_var("N", "N|-DIM-N-CL|V", "V") %>%
           lab_var("N", "N|-V-CL", "V") %>%
           lab_var("PRO", "PRO:POSS|-N-CL|V", "V") %>%
           lab_var("V", "V|-CL|PRO", "PRO") %>%
           lab_var("V", "V|-PROG~INF") %>%
           lab_var("ADV", "WH:ADV|-N-CL|V", "V") ) 

# extract most frequent cat converting plural common nouns and verbs to root
type_cat_root <- function(type) {
  if (type %in% mot_mor_filter$word) {
    mot_mor_filter %>%
      filter(word == type) %>%
      group_by(mor_raw) %>%
      summarise(n = n()) %>%
      arrange(desc(n)) %>%
      filter(row_number() == 1) %>%
      (function(z) {z$mor_raw}) %>%
      (function(z) {mot_mor_filter %>%
          filter(mor_raw == z & word == type) %>%
          filter(row_number() == 1) %>%
          (function(t) {
            if (t$pl == TRUE | t$verb == TRUE) {
              t$root %>% unlist() 
            } else {
              t$word
            }
          })})
  } else {
    NA
  }
}

# types root from uni_diff in mot_uni_cat
type_to_root <- function(df) {
  df %>%
    mutate(types_root = uni_diff %>%
             sapply(function(x) {
               sapply(x, type_cat_root)
             }) %>%
             sapply(function(x) {
               x %>%
                 unlist() %>%
                 unique() %>%
                 sort()
             }))
} 

mot_uni_cat %<>%
  type_to_root()

# redo unique types in types_root
mot_uni_cat %<>%
  (function(x) {
    x$types_root <- setdif_list(x, "types_root")
    x
  })

# assign cat from mot_mor_uni
assign_cat <- function(df) {
  df %>%
    mutate(cat = types_root %>%
             sapply(function(x) {
               mot_mor_uni$cat[fmatch(x, mot_mor_uni$word)] %>%
                 (function(y) {
                   y[unlist(lapply(y , is.null))] <- NA 
                   y %>% 
                     unlist()
                 })
             })) %>%
    mutate(cat_prop = lapply(cat, function(x) {
      x <- table(x)
      perc <- prop.table(x) * 100
      
      cats <- c("N", "V", "ADJ", "ADV", "PRO")
      
      x[names(x) %in% cats] %>%
        sort(decreasing = TRUE) %>%
        (function(x) {
          cats <- names(x)
          
          x <- tibble(cat = names(x), freq = x, perc = perc[cats])
          x[x == "N"] <- "NOUN"
          x[x == "V"] <- "VERB"
          x[x == "PRO"] <- "PRON"
          x
        }) %>%
        mutate(freq = as.numeric(freq),
               perc = as.numeric(perc))
    }))  %>%
    ungroup() %>%
    mutate(noun_perc = sapply(cat_prop, function(x) {x$perc[which(x$cat == "NOUN")]}),
           verb_perc = sapply(cat_prop, function(x) {x$perc[which(x$cat == "VERB")]}),
           pron_perc = sapply(cat_prop, function(x) {x$perc[which(x$cat == "PRON")]}),
           adj_perc = sapply(cat_prop, function(x) {x$perc[which(x$cat == "ADJ")]}),
           adv_perc = sapply(cat_prop, function(x) {x$perc[which(x$cat == "ADV")]})) %>%
    select(-cat_prop) %>%
    mutate(adv_perc = sapply(adv_perc, function(x) {
      ifelse(length(x) == 0, 0, x)
    }),
    pron_perc = sapply(pron_perc, function(x) {
      ifelse(length(x) == 0, 0, x)
    }),
    adj_perc = sapply(adj_perc, function(x) {
      ifelse(length(x) == 0, 0, x)
    }),
    verb_perc = sapply(verb_perc, function(x) {
      ifelse(length(x) == 0, 0, x)
    }),
    noun_perc = sapply(noun_perc, function(x) {
      ifelse(length(x) == 0, 0, x)
    })) %>%
    group_by(baby) %>%
    mutate(noun_perc_cum = cumsum(noun_perc),
           verb_perc_cum = cumsum(verb_perc),
           adj_perc_cum = cumsum(adj_perc),
           adv_perc_cum = cumsum(adv_perc),
           pron_perc_cum = cumsum(pron_perc)) %>%
    group_by(baby) %>%
    mutate(noun_perc_cum = ( noun_perc_cum / seq(100, 2000, by = 100) ) * 100,
           verb_perc_cum = ( verb_perc_cum / seq(100, 2000, by = 100) ) * 100,
           adj_perc_cum = ( adj_perc_cum / seq(100, 2000, by = 100)  ) * 100,
           adv_perc_cum = ( adv_perc_cum / seq(100, 2000, by = 100) ) * 100,
           pron_perc_cum = ( pron_perc_cum / seq(100, 2000, by = 100) ) * 100)
  
}

mot_uni_cat %<>%
  assign_cat()

table_cat <- function(df) {
  df %>%
    filter(section == 20) %>%
    ungroup() %>%
    select(baby, noun_perc_cum:pron_perc_cum) %>%
    rename(`Noun (%)`=noun_perc_cum,
           `Verb (%)`=verb_perc_cum,
           `Pron (%)`=pron_perc_cum,
           `Adj (%)`=adj_perc_cum,
           `Adv (%)`=adv_perc_cum) %>%
    (function(x) {
      x %>%
        rbind(tibble(baby = "MEAN",
                     `Noun (%)` = mean(x$`Noun (%)`),
                     `Verb (%)` = mean(x$`Verb (%)`),
                     `Pron (%)` = mean(x$`Pron (%)`),
                     `Adj (%)` = mean(x$`Adj (%)`),
                     `Adv (%)` = mean(x$`Adv (%)`)))
    }) %>%
    round_df(2)
  
}

mot_uni_cat_table <- mot_uni_cat %>% 
  table_cat()
  
#### CHILDREN ####
# plurals and verbs to root
chi_uni_cat <- chi_uni %>%
  select(baby, section, uni_diff) %>%
  type_to_root()

# redo unique types
chi_uni_cat %<>%
  (function(x) {
    x$types_root <- setdif_list(x, "types_root")
    x
  })

# assign cats
chi_uni_cat %<>%
  assign_cat()

chi_uni_cat_table <- chi_uni_cat %>%
  table_cat()

#### MODEL ####
# assign phon to mot_na_baby_section
mot_na_baby_section %<>%
  mutate(phon = mot_phon$phon[fmatch(word, mot_phon$word)])

# mod phons to orts
mod_uni_cat <- mod %>%
  select(baby:phon)

phon_to_ort <- function(df) {
  df$uni_diff <- df$phon
  
  for (i in seq_along(df$section)) {
    for (j in seq_along(df$phon[[i]])) {
      ort <- df$phon[[i]][j] %>%
        (function(x) {
          mot_na_baby_section %>%
            filter(phon == x) %>%
            (function(y) {
              y$word %>%
                table() %>%
                sort(decreasing = T) %>%
                (function(z) {
                  names(z[1])
                })
            })
        })
      
      if (length(ort) == 0) {
        df$uni_diff[[i]][j] <- NA
      } else {
        df$uni_diff[[i]][j] <- ort
      }
    }
  }
  
  df
}

mod_uni_cat %<>%
  phon_to_ort()

# re run mot/chi script
mod_uni_cat %<>%
  select(baby, section, uni_diff) %>%
  type_to_root()

# redo unique types
mod_uni_cat %<>%
  (function(x) {
    x$types_root <- setdif_list(x, "types_root")
    x
  })

# assign cats
mod_uni_cat %<>%
  assign_cat()

mod_uni_cat_table <- mod_uni_cat %>%
  table_cat()

#### save mot nouns not in chi and vice versa ####
types_root <- function(df) {
  df %>%
    (function(x) {
      new_df <- tibble(types_root = unlist(x$types_root) %>% str_replace("\\\\", ""), 
                       cat = unlist(x$cat)) %>%
        filter(cat == "N")
      new_df[!duplicated(new_df$types_root),] %>%
        arrange(types_root) %>%
        select(-cat) %>%
        unlist()
    })
}

intersect_types_root <- tibble(mot_nouns_root = types_root(mot_uni_cat),
       chi_nouns_root = types_root(mot_uni_cat) %in% types_root(chi_uni_cat)) %>%
  mutate(chi_nouns_root = ifelse(chi_nouns_root == TRUE, mot_nouns_root, NA))

write_csv(intersect_types_root, "mot_chi_nouns.csv", na = "")

write_lines((intersect_types_root %>%
              filter(is.na(chi_nouns_root)))$mot_nouns_root, 
            "mot_NOT_learned.txt")

write_lines((intersect_types_root %>%
               filter(!is.na(chi_nouns_root)))$mot_nouns_root, 
            "mot_learned.txt")

#### Mot bisyllabic nouns analysis ####
# extract (1) MOT bisyllabic noun types (bis) and (2) MOT monosyllabic (mono) types
# convert (1) and (2) to phonetic
mot_nouns_bis <- mot_uni_cat %>%
  select(baby, types_root, cat) %>%
  group_by(baby) %>%
  summarise(types_root = list(unlist(types_root)),
            cat = list(unlist(cat))) %>%
  apply(1, function(x) {
    tibble(id = x$baby, word = unlist(x$types_root), cat = unlist(x$cat))
  }) %>% 
  rbindlist() %>%
  inner_join(., mot_phon, "word") %>%
  mutate(syllable = str_split(phon, "_") %>% 
           sapply(function(x) {
             x %in% vowels %>% 
               sum()
           })) %>%
  filter(cat == "N") %>%
  filter(syllable == 2)

mot_mono <- mot_uni %>%
  select(baby, uni_diff) %>%
  group_by(baby) %>%
  summarise(uni_diff = list(unlist(uni_diff))) %>%
  apply(1, function(x) {
    tibble(id = x$baby, word = unlist(x$uni_diff))
  }) %>% 
  rbindlist() %>%
  inner_join(., mot_phon, "word") %>%
  mutate(syllable = str_split(phon, "_") %>% 
           sapply(function(x) {
             x %in% vowels %>% 
               sum()
           })) %>%
  filter(syllable == 1) %>%
  arrange(id, phon)

# get rid of deplicates in phon
mot_mono %<>% 
  group_by(id, phon) %>%
  slice(1) %>%
  ungroup()

# split (1) into learned (1a) and not-learned (1b)
mot_nouns_bis_learned <- mot_nouns_bis %>%
  (function(x) {
    lapply(unique(x$id), function(y) {
      reference <- (tibble(word = (chi_uni_cat %>% filter(baby == y))$types_root %>% unlist(), 
                          cat = (chi_uni_cat %>% filter(baby == y))$cat %>% unlist()))$word
      
      mot_nouns_bis %>%
        filter(id == y) %>%
        filter(word %in% reference)
    }) %>% rbindlist()
  })

mot_nouns_bis_not_learned <- mot_nouns_bis %>%
  (function(x) {
    lapply(unique(x$id), function(y) {
      reference <- (tibble(word = (chi_uni_cat %>% filter(baby == y))$types_root %>% unlist(), 
                           cat = (chi_uni_cat %>% filter(baby == y))$cat %>% unlist()))$word
      
      mot_nouns_bis %>%
        filter(id == y) %>%
        filter(!word %in% reference)
    }) %>% rbindlist()
  })

# assign to each bis how many mono contains (grouped by id)
mot_nouns_bis_learned %<>%
  mutate(monosyll_type = sapply(unique(id), function(name) {
  sapply((mot_nouns_bis_learned %>% filter(id == name))$phon, function(x) {
    x %>%
      str_detect((mot_mono %>% filter(id == name))$phon) %>%
      sum()
    })
    }) %>%
    unlist())

mot_nouns_bis_not_learned %<>%
  mutate(monosyll_type = sapply(unique(id), function(name) {
    sapply((mot_nouns_bis_not_learned %>% filter(id == name))$phon, function(x) {
      x %>%
        str_detect((mot_mono %>% filter(id == name))$phon) %>%
        sum()
    })
  }) %>%
    unlist())

# sum number of mono in (1a) and (1b) and take average by id
summary_bysillables <- function(df1, df2, type = TRUE) {
 if (type == TRUE) {
   df1 %>%
     mutate(monosyll = monosyll_type %>% (function(x) {
       names(x) <- NULL
       x
     })) %>%
     group_by(id) %>%
     summarize(tot_mono_learned = mean(monosyll)) %>%
     bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
     bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD"))) %>%
     cbind(tot_mono_not_learned = df2 %>%
             mutate(monosyll = monosyll_type %>% (function(x) {
               names(x) <- NULL
               x
             })) %>%
             group_by(id) %>%
             summarize(tot_mono_not_learned = mean(monosyll)) %>%
             bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
             bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD"))) %>%
             select(tot_mono_not_learned))
 } else {
   df1 %>%
     mutate(monosyll = monosyll_tokens %>% (function(x) {
       names(x) <- NULL
       x
     })) %>%
     group_by(id) %>%
     summarize(tot_mono_learned = mean(monosyll)) %>%
     bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
     bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD"))) %>%
     cbind(tot_mono_not_learned = df2 %>%
             mutate(monosyll = monosyll_tokens %>% (function(x) {
               names(x) <- NULL
               x
             })) %>%
             group_by(id) %>%
             summarize(tot_mono_not_learned = mean(monosyll)) %>%
             bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
             bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD"))) %>%
             select(tot_mono_not_learned))
 }
  
  
}

summary_bis <- list(id = summary_bysillables(mot_nouns_bis_learned,
                                             mot_nouns_bis_not_learned, 
                                             type = TRUE))

#### mono types in bis by stage ####
# bisyllabic words list by stage
mot_nouns_bis <- mot_uni_cat %>%
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
  filter(cat == "N") %>%
  filter(syllable == 2)

# monosyllabic words list by stage 
mot_mono <- mot_uni %>%
  select(baby, section, uni_diff) %>%
  apply(1, function(x) {
    tibble(id = x$baby, section = x$section, word = unlist(x$uni_diff))
  }) %>% 
  rbindlist() %>%
  inner_join(., mot_phon, "word") %>%
  mutate(syllable = str_split(phon, "_") %>% 
           sapply(function(x) {
             x %in% vowels %>% 
               sum()
           })) %>%
  filter(syllable == 1) %>%
  arrange(id, phon)

# split into learned not-learned
mot_nouns_bis %<>%
  mutate(learned = unique(id) %>%
           sapply(function(x) {
             mot_nouns_bis %>%
               filter(id == x) %>%
               .$word %>%
               unlist() %>%
               {. %in% (chi_uni_cat %>%
                          filter(baby == x) %>%
                          .$types_root %>%
                          unlist())}
           }) %>%
           unlist())

mot_nouns_bis_learned <- mot_nouns_bis %>%
  filter(learned == TRUE)

chi_learned <- chi_uni_cat %>% 
  apply(1, function(x) {
    tibble(id = x$baby, section = x$section, 
           word = unlist(x$types_root))
  }) %>%
  rbindlist() 

for (ids in unique(mot_nouns_bis_learned$id)) {
  for (w in filter(mot_nouns_bis_learned, id == ids)$word) {
    mot_nouns_bis_learned$section[which(mot_nouns_bis_learned$id == ids & mot_nouns_bis_learned$word == w)] <-
      chi_learned %>%
      filter(id == ids, word == w) %>%
      .$section
  }
}

mot_nouns_bis_not_learned <- mot_nouns_bis %>%
  filter(learned == FALSE)

# assign mono in bis by stage
assign_mono_type <- function(df) {
  df %>%
    mutate(monosyll_type = apply(df %>% 
                              distinct(id, section), 1, function(y) {
                                (df %>% filter(id == y[1], section == as.numeric(y[2])))$phon %>%
                                  sapply(function(x) {
                                    x %>%
                                      str_detect((mot_mono %>% filter(id == y[1], section %in% 0:as.numeric(y[2])))$phon %>% 
                                                   unique()) %>%
                                      sum()
                                  })
                              }) %>%
             unlist())
}

mot_nouns_bis_learned %<>% assign_mono_type()
mot_nouns_bis_not_learned %<>% assign_mono_type()

# mean for each child and total mean
summary_bis %<>% 
  c(list(id_stage = summary_bysillables(mot_nouns_bis_learned, mot_nouns_bis_not_learned, 
                                        type = TRUE)))

#### mono tokens in bis by stage ####
# create mono tokens by section
mot_mono_tokens <- mot_na_baby_section %>%
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
  filter(syllable == 1)

# assign mono
for (i in 1:nrow(mot_nouns_bis_learned)) {
  mot_nouns_bis_learned$monosyll_tokens[i] <- 
    mot_mono_tokens %>%
    filter(baby == mot_nouns_bis_learned$id[i]) %>%
    filter(section %in% 0:mot_nouns_bis_learned$section[i]) %>%
    (function(x) {
      str_detect(mot_nouns_bis_learned$phon[i], x$phon) %>%
        sum()
    })
}

for (i in 1:nrow(mot_nouns_bis_not_learned)) {
  mot_nouns_bis_not_learned$monosyll_tokens[i] <- 
    mot_mono_tokens %>%
    filter(baby == mot_nouns_bis_not_learned$id[i]) %>%
    filter(section %in% 0:mot_nouns_bis_not_learned$section[i]) %>%
    (function(x) {
      str_detect(mot_nouns_bis_not_learned$phon[i], x$phon) %>%
        sum()
    })
}

# add count to summary
summary_bis %<>% 
  c(list(id_stage_tokens = summary_bysillables(mot_nouns_bis_learned, mot_nouns_bis_not_learned, 
                                        type = FALSE)))

#### phonological frames (tokens) in each mono, bi and tri word ####
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
  })

# for each bis, calculate all possible adjacent combs and match in mot_tokens
adj_comb <- function(Data) {
  A <- lapply(2:(length(Data)-1), sequence)
  B <- lapply(rev(vapply(A, length, 1L))-1, function(x) c(0, sequence(x)))
  unlist(lapply(seq_along(A), function(x) {
    lapply(B[[x]], function(y) Data[A[[x]]+y])
  }), recursive = FALSE, use.names = FALSE)
}

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
           sapply(adj_comb)) %>%
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

mot_nouns_bis_learned %<>%
  seq2_5()

mot_nouns_bis_not_learned %<>%
  seq2_5()

mot_tokens %<>%
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
  select(baby.x:phonemes.x, split.y:seq5) %>%
  (function(x) {
    colnames(x) <- str_remove(colnames(x), "[.]{1}.*$")
    x
  })

# unique sequences in seq2 of mot_tokens
mot_tokens %<>%
  mutate(seq2 = sapply(seq2, unique))

# check how many times sequences appear in mot_tokens
assign_seq <- function(df, seq_num, df_ref, df_count = mot_tokens) {
  apply(df, 1, function(rows) {
    sapply((df_ref %>%
              filter(id == rows[1], section == as.numeric(rows[2])) %>%
              .[[seq_num]]),  function(x) {
                (df_count %>%
                   filter(baby == rows[1], section %in% 0:as.numeric(rows[2])) %>%
                   .[[seq_num]] %>%
                   unlist()) %>%
                  str_count(paste("^", x, "$", collapse = "|", sep="")) %>% sum() 
              }) %>% mean(na.rm = TRUE)
    })
}

summary_seq_tokens <- function(df) {
  df %>%
    arrange(id, section) %>%
    group_by(id) %>%
    summarise(mean_seq2 = mean(seq2, na.rm = TRUE),
              mean_seq3 = mean(seq3, na.rm = TRUE),
              mean_seq4 = mean(seq4, na.rm = TRUE),
              mean_seq5 = mean(seq5, na.rm = TRUE)) %>%
    bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(., na.rm = TRUE) else "MEAN"))) %>%
    bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(., na.rm = TRUE) else "SD")))
}

summary_seq_tokens_learned <- mot_nouns_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_learned),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_learned),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_learned),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_learned)) %>%
  summary_seq_tokens()

summary_seq_tokens_not_learned <- mot_nouns_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_not_learned),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_not_learned),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_not_learned),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_not_learned)) %>%
  summary_seq_tokens()

# same calculations but for monosyllabic and trisyllabic words!
# bisyllabic words list by stage
mot_nouns <- mot_uni_cat %>%
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
  filter(cat == "N")

mot_nouns_mono <- mot_nouns %>%
  filter(syllable == 1)

mot_nouns_tri <- mot_nouns %>%
  filter(syllable == 3)

# MONO 
# split into learned not-learned
mot_nouns_mono %<>%
  mutate(learned = unique(id) %>%
           sapply(function(x) {
             mot_nouns_mono %>%
               filter(id == x) %>%
               .$word %>%
               unlist() %>%
               {. %in% (chi_uni_cat %>%
                          filter(baby == x) %>%
                          .$types_root %>%
                          unlist())}
           }) %>%
           unlist())

mot_nouns_mono_learned <- mot_nouns_mono %>%
  filter(learned == TRUE)

chi_learned <- chi_uni_cat %>% 
  apply(1, function(x) {
    tibble(id = x$baby, section = x$section, 
           word = unlist(x$types_root))
  }) %>%
  rbindlist() 

for (ids in unique(mot_nouns_mono_learned$id)) {
  for (w in filter(mot_nouns_mono_learned, id == ids)$word) {
    mot_nouns_mono_learned$section[which(mot_nouns_mono_learned$id == ids & mot_nouns_mono_learned$word == w)] <-
      chi_learned %>%
      filter(id == ids, word == w) %>%
      .$section
  }
}

mot_nouns_mono_not_learned <- mot_nouns_mono %>%
  filter(learned == FALSE)

mot_nouns_mono_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

mot_nouns_mono_not_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

summary_seq_tokens_mono_learned <- mot_nouns_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_learned),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_learned),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_learned),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_learned)) %>%
  summary_seq_tokens()

summary_seq_tokens_mono_not_learned <- mot_nouns_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_not_learned),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_not_learned),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_not_learned),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_not_learned)) %>%
  summary_seq_tokens()

# TRI
# split into learned not-learned
mot_nouns_tri %<>%
  mutate(learned = unique(id) %>%
           sapply(function(x) {
             mot_nouns_tri %>%
               filter(id == x) %>%
               .$word %>%
               unlist() %>%
               {. %in% (chi_uni_cat %>%
                          filter(baby == x) %>%
                          .$types_root %>%
                          unlist())}
           }) %>%
           unlist())

mot_nouns_tri_learned <- mot_nouns_tri %>%
  filter(learned == TRUE)

chi_learned <- chi_uni_cat %>% 
  apply(1, function(x) {
    tibble(id = x$baby, section = x$section, 
           word = unlist(x$types_root))
  }) %>%
  rbindlist() 

for (ids in unique(mot_nouns_tri_learned$id)) {
  for (w in filter(mot_nouns_tri_learned, id == ids)$word) {
    mot_nouns_tri_learned$section[which(mot_nouns_tri_learned$id == ids & mot_nouns_tri_learned$word == w)] <-
      chi_learned %>%
      filter(id == ids, word == w) %>%
      .$section
  }
}

mot_nouns_tri_not_learned <- mot_nouns_tri %>%
  filter(learned == FALSE)

mot_nouns_tri_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

mot_nouns_tri_not_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

summary_seq_tokens_tri_learned <- mot_nouns_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_learned),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_learned),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_learned),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_learned)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_not_learned <- mot_nouns_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_not_learned),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_not_learned),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_not_learned),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_not_learned)) %>%
  summary_seq_tokens()

#### frequency of learned/not learned words up to current stage ####
freq_up_stage <- function(df) {
  sapply(1:nrow(df), function(i) {
    mot_na_baby_section %>%
      filter(word == df$word[i],
             baby == df$id[i],
             section %in% 0:df$section[i]) %>%
      nrow()
  })
}

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

#### (1) learned not learned matched on frequency ####
summary_freq <- function(df, name_var) {
  df %>%
    group_by(id) %>%
    summarize(!!name_var := mean(freq))
}

summary_freq_sets <- summary_freq(mot_nouns_mono_learned, "Mono_learned_freq") %>%
  cbind(summary_freq(mot_nouns_mono_not_learned, "Mono_not_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_bis_learned, "Bis_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_bis_not_learned, "Bis_not_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_tri_learned, "Tri_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_nouns_tri_not_learned, "Tri_not_learned_freq")[,2]) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD")))

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

temp <- tibble(dfs = c("mot_nouns_mono_learned",
               "mot_nouns_bis_learned",
               "mot_nouns_tri_learned"),
       vars = colnames(summary_freq_sets)[-1] %>% .[grepl("not",.)])

matched_dfs <- list()
for (rows in 1:nrow(temp)) {
  matched_dfs[[paste(temp$dfs[rows], "match", sep = "_")]] <-
    unique(get(temp$dfs[rows])$id) %>%
    lapply(function(x) {
      equalize_mean(get(temp$dfs[rows]), x, temp$vars[rows])
    }) %>%
    rbindlist()
} ; rm(rows, temp)
  
# repeat calculation for new dfs
summary_seq_tokens_mono_learned_match <- matched_dfs[["mot_nouns_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_mono_learned_match"]]),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_mono_learned_match"]]),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_mono_learned_match"]]),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_mono_learned_match"]])) %>%
  summary_seq_tokens()

summary_seq_tokens_bis_learned_match <- matched_dfs[["mot_nouns_bis_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_bis_learned_match"]]),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_bis_learned_match"]]),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_bis_learned_match"]]),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_bis_learned_match"]])) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_learned_match <- matched_dfs[["mot_nouns_tri_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_tri_learned_match"]]),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_tri_learned_match"]]),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_tri_learned_match"]]),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_tri_learned_match"]])) %>%
  summary_seq_tokens()

#### phonotactic frames counted in children vocabulary  ####
# children's utterances, adding section and long formatting
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

# sequences 
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

summary_seq_tokens_mono_learned_chi <- mot_nouns_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_mono_not_learned_chi <- mot_nouns_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_not_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_not_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_not_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_not_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_learned_chi <- mot_nouns_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_not_learned_chi <- mot_nouns_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_not_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_not_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_not_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_not_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_learned_chi <- mot_nouns_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_not_learned_chi <- mot_nouns_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_not_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_not_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_not_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_not_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_mono_learned_match_chi <- matched_dfs[["mot_nouns_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_bis_learned_match_chi <- matched_dfs[["mot_nouns_bis_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_learned_match_chi <- matched_dfs[["mot_nouns_tri_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_tokens)) %>%
  summary_seq_tokens()

# use proportions for summaries1
adorn_perc_summary <- function(df) {
  df %<>%
    .[1:12,] %>%
    adorn_percentages() %>%
    bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
    bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD")))
  
}

summary_seq_tokens_bis_learned_match %<>%
  adorn_perc_summary
summary_seq_tokens_bis_learned_match_chi %<>%
  adorn_perc_summary
summary_seq_tokens_learned %<>%
  adorn_perc_summary
summary_seq_tokens_learned_chi %<>%
  adorn_perc_summary
summary_seq_tokens_mono_learned %<>%
  adorn_perc_summary
summary_seq_tokens_mono_learned_chi %<>%
  adorn_perc_summary
summary_seq_tokens_mono_learned_match %<>%
  adorn_perc_summary
summary_seq_tokens_mono_learned_match_chi %<>%
  adorn_perc_summary
summary_seq_tokens_mono_not_learned %<>%
  adorn_perc_summary
summary_seq_tokens_mono_not_learned_chi %<>%
  adorn_perc_summary
summary_seq_tokens_not_learned %<>%
  adorn_perc_summary
summary_seq_tokens_not_learned_chi %<>%
  adorn_perc_summary
summary_seq_tokens_tri_learned %<>%
  adorn_perc_summary
summary_seq_tokens_tri_learned_chi %<>%
  adorn_perc_summary
summary_seq_tokens_tri_learned_match %<>%
  adorn_perc_summary
summary_seq_tokens_tri_learned_match_chi %<>%
  adorn_perc_summary
summary_seq_tokens_tri_not_learned %<>%
  adorn_perc_summary
summary_seq_tokens_tri_not_learned_chi %<>%
  adorn_perc_summary

#### mot_types phonotactic frames in mot_nouns ####
# creare mot_types from mot_tokens
mot_types <- mot_tokens %>%
  arrange(baby, section) %>%
  group_by(baby) %>%
  distinct(phon, .keep_all = TRUE) 

# ripetere i summary for all datasets (matched and not) using mot_types
summary_seq_types_mono_learned <- mot_nouns_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_mono_not_learned <- mot_nouns_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_not_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_not_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_not_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_not_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_learned <- mot_nouns_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_not_learned <- mot_nouns_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_not_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_not_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_not_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_not_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_learned <- mot_nouns_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_not_learned <- mot_nouns_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_not_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_not_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_not_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_not_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_mono_learned_match <- matched_dfs[["mot_nouns_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = mot_types),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = mot_types),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = mot_types),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_bis_learned_match <- matched_dfs[["mot_nouns_bis_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = mot_types),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = mot_types),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = mot_types),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_learned_match <- matched_dfs[["mot_nouns_tri_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = mot_types),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = mot_types),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = mot_types),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = mot_types)) %>%
  summary_seq_tokens()


#### chi_types phonotactic frames in mot_nouns ####
# creare chi_types from chi_tokens
chi_types <- chi_tokens %>%
  arrange(baby, section) %>%
  group_by(baby) %>%
  distinct(phon, .keep_all = TRUE) 

# ripetere i summary for all datasets (matched and not) using chi_types
summary_seq_types_mono_learned_chi <- mot_nouns_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_mono_not_learned_chi <- mot_nouns_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_mono_not_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_nouns_mono_not_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_nouns_mono_not_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_nouns_mono_not_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_learned_chi <- mot_nouns_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_not_learned_chi <- mot_nouns_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_bis_not_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_nouns_bis_not_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_nouns_bis_not_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_nouns_bis_not_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_learned_chi <- mot_nouns_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_not_learned_chi <- mot_nouns_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_nouns_tri_not_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_nouns_tri_not_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_nouns_tri_not_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_nouns_tri_not_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_mono_learned_match_chi <- matched_dfs[["mot_nouns_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_types),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_types),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_types),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_mono_learned_match"]], df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_bis_learned_match_chi <- matched_dfs[["mot_nouns_bis_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_types),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_types),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_types),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_bis_learned_match"]], df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_learned_match_chi <- matched_dfs[["mot_nouns_tri_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_types),
         seq3 = assign_seq(., "seq3", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_types),
         seq4 = assign_seq(., "seq4", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_types),
         seq5 = assign_seq(., "seq5", matched_dfs[["mot_nouns_tri_learned_match"]], df_count = chi_types)) %>%
  summary_seq_tokens()

# to proportions
summary_seq_types_bis_learned_match %<>%
  adorn_perc_summary
summary_seq_types_bis_learned_match_chi %<>%
  adorn_perc_summary
summary_seq_types_learned %<>%
  adorn_perc_summary
summary_seq_types_learned_chi %<>%
  adorn_perc_summary
summary_seq_types_mono_learned %<>%
  adorn_perc_summary
summary_seq_types_mono_learned_chi %<>%
  adorn_perc_summary
summary_seq_types_mono_learned_match %<>%
  adorn_perc_summary
summary_seq_types_mono_learned_match_chi %<>%
  adorn_perc_summary
summary_seq_types_mono_not_learned %<>%
  adorn_perc_summary
summary_seq_types_mono_not_learned_chi %<>%
  adorn_perc_summary
summary_seq_types_not_learned %<>%
  adorn_perc_summary
summary_seq_types_not_learned_chi %<>%
  adorn_perc_summary
summary_seq_types_tri_learned %<>%
  adorn_perc_summary
summary_seq_types_tri_learned_chi %<>%
  adorn_perc_summary
summary_seq_types_tri_learned_match %<>%
  adorn_perc_summary
summary_seq_types_tri_learned_match_chi %<>%
  adorn_perc_summary
summary_seq_types_tri_not_learned %<>%
  adorn_perc_summary
summary_seq_types_tri_not_learned_chi %<>%
  adorn_perc_summary


#### mot_tokens phonotactic frames in mot_verbs ####
# create mot_verbs (learned and not learned)
mot_verbs <- mot_uni_cat %>%
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
  filter(cat == "V")

mot_verbs_mono <- mot_verbs %>%
  filter(syllable == 1)

mot_verbs_bis <- mot_verbs %>%
  filter(syllable == 2)

mot_verbs_tri <- mot_verbs %>%
  filter(syllable == 3)

# mono verbs
mot_verbs_mono %<>%
  mutate(learned = unique(id) %>%
           sapply(function(x) {
             mot_verbs_mono %>%
               filter(id == x) %>%
               .$word %>%
               unlist() %>%
               {. %in% (chi_uni_cat %>%
                          filter(baby == x) %>%
                          .$types_root %>%
                          unlist())}
           }) %>%
           unlist())

mot_verbs_mono_learned <- mot_verbs_mono %>%
  filter(learned == TRUE)

chi_learned <- chi_uni_cat %>% 
  apply(1, function(x) {
    tibble(id = x$baby, section = x$section, 
           word = unlist(x$types_root))
  }) %>%
  rbindlist() 

for (ids in unique(mot_verbs_mono_learned$id)) {
  for (w in filter(mot_verbs_mono_learned, id == ids)$word) {
    mot_verbs_mono_learned$section[which(mot_verbs_mono_learned$id == ids & mot_verbs_mono_learned$word == w)] <-
      chi_learned %>%
      filter(id == ids, word == w) %>%
      .$section
  }
}

mot_verbs_mono_not_learned <- mot_verbs_mono %>%
  filter(learned == FALSE)

mot_verbs_mono_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

mot_verbs_mono_not_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

# bis verbs
mot_verbs_bis %<>%
  mutate(learned = unique(id) %>%
           sapply(function(x) {
             mot_verbs_bis %>%
               filter(id == x) %>%
               .$word %>%
               unlist() %>%
               {. %in% (chi_uni_cat %>%
                          filter(baby == x) %>%
                          .$types_root %>%
                          unlist())}
           }) %>%
           unlist())

mot_verbs_bis_learned <- mot_verbs_bis %>%
  filter(learned == TRUE)

chi_learned <- chi_uni_cat %>% 
  apply(1, function(x) {
    tibble(id = x$baby, section = x$section, 
           word = unlist(x$types_root))
  }) %>%
  rbindlist() 

for (ids in unique(mot_verbs_bis_learned$id)) {
  for (w in filter(mot_verbs_bis_learned, id == ids)$word) {
    mot_verbs_bis_learned$section[which(mot_verbs_bis_learned$id == ids & mot_verbs_bis_learned$word == w)] <-
      chi_learned %>%
      filter(id == ids, word == w) %>%
      .$section
  }
}

mot_verbs_bis_not_learned <- mot_verbs_bis %>%
  filter(learned == FALSE)

mot_verbs_bis_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

mot_verbs_bis_not_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

# tri verbs
mot_verbs_tri %<>%
  mutate(learned = unique(id) %>%
           sapply(function(x) {
             mot_verbs_tri %>%
               filter(id == x) %>%
               .$word %>%
               unlist() %>%
               {. %in% (chi_uni_cat %>%
                          filter(baby == x) %>%
                          .$types_root %>%
                          unlist())}
           }) %>%
           unlist())

mot_verbs_tri_learned <- mot_verbs_tri %>%
  filter(learned == TRUE)

chi_learned <- chi_uni_cat %>% 
  apply(1, function(x) {
    tibble(id = x$baby, section = x$section, 
           word = unlist(x$types_root))
  }) %>%
  rbindlist() 

for (ids in unique(mot_verbs_tri_learned$id)) {
  for (w in filter(mot_verbs_tri_learned, id == ids)$word) {
    mot_verbs_tri_learned$section[which(mot_verbs_tri_learned$id == ids & mot_verbs_tri_learned$word == w)] <-
      chi_learned %>%
      filter(id == ids, word == w) %>%
      .$section
  }
}

mot_verbs_tri_not_learned <- mot_verbs_tri %>%
  filter(learned == FALSE)

mot_verbs_tri_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

mot_verbs_tri_not_learned %<>%
  filter(str_detect(phon, "_")) %>%
  seq2_5()

# summaries verbs
summary_seq_tokens_mono_learned_verbs <- mot_verbs_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_learned, df_count = mot_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_learned, df_count = mot_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_learned, df_count = mot_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_learned, df_count = mot_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_mono_not_learned_verbs <- mot_verbs_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_not_learned, df_count = mot_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_not_learned, df_count = mot_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_not_learned, df_count = mot_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_not_learned, df_count = mot_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_learned_verbs <- mot_verbs_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_learned, df_count = mot_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_learned, df_count = mot_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_learned, df_count = mot_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_learned, df_count = mot_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_not_learned_verbs <- mot_verbs_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_not_learned, df_count = mot_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_not_learned, df_count = mot_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_not_learned, df_count = mot_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_not_learned, df_count = mot_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_learned_verbs <- mot_verbs_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_learned, df_count = mot_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_learned, df_count = mot_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_learned, df_count = mot_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_learned, df_count = mot_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_not_learned_verbs <- mot_verbs_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_not_learned, df_count = mot_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_not_learned, df_count = mot_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_not_learned, df_count = mot_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_not_learned, df_count = mot_tokens)) %>%
  summary_seq_tokens()

# to proportions
summary_seq_tokens_learned_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_not_learned_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_mono_learned_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_mono_not_learned_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_tri_learned_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_tri_not_learned_verbs %<>%
  adorn_perc_summary

#### assign and treat frequency (verbs) ####
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

#### chi_tokens phonotactic frames in mot_verbs ####
# summaries verbs
summary_seq_tokens_mono_learned_chi_verbs <- mot_verbs_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_mono_not_learned_chi_verbs <- mot_verbs_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_not_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_not_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_not_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_not_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_learned_chi_verbs <- mot_verbs_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_not_learned_chi_verbs <- mot_verbs_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_not_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_not_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_not_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_not_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_learned_chi_verbs <- mot_verbs_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_tokens_tri_not_learned_chi_verbs <- mot_verbs_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_not_learned, df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_not_learned, df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_not_learned, df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_not_learned, df_count = chi_tokens)) %>%
  summary_seq_tokens()

# to proportions
summary_seq_tokens_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_not_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_mono_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_mono_not_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_tri_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_tokens_tri_not_learned_chi_verbs %<>%
  adorn_perc_summary

#### mot_types phonotactic frames in mot_verbs ####
# summaries verbs
summary_seq_types_mono_learned_verbs <- mot_verbs_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_mono_not_learned_verbs <- mot_verbs_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_not_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_not_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_not_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_not_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_learned_verbs <- mot_verbs_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_not_learned_verbs <- mot_verbs_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_not_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_not_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_not_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_not_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_learned_verbs <- mot_verbs_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_not_learned_verbs <- mot_verbs_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_not_learned, df_count = mot_types),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_not_learned, df_count = mot_types),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_not_learned, df_count = mot_types),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_not_learned, df_count = mot_types)) %>%
  summary_seq_tokens()

# to proportions
summary_seq_types_learned_verbs %<>%
  adorn_perc_summary
summary_seq_types_not_learned_verbs %<>%
  adorn_perc_summary
summary_seq_types_mono_learned_verbs %<>%
  adorn_perc_summary
summary_seq_types_mono_not_learned_verbs %<>%
  adorn_perc_summary
summary_seq_types_tri_learned_verbs %<>%
  adorn_perc_summary
summary_seq_types_tri_not_learned_verbs %<>%
  adorn_perc_summary

#### chi_types phonotactic frames in mot_verbs ####
# summaries verbs
summary_seq_types_mono_learned_chi_verbs <- mot_verbs_mono_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_mono_not_learned_chi_verbs <- mot_verbs_mono_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_mono_not_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_verbs_mono_not_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_verbs_mono_not_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_verbs_mono_not_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_learned_chi_verbs <- mot_verbs_bis_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_not_learned_chi_verbs <- mot_verbs_bis_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_bis_not_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_verbs_bis_not_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_verbs_bis_not_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_verbs_bis_not_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_learned_chi_verbs <- mot_verbs_tri_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_types_tri_not_learned_chi_verbs <- mot_verbs_tri_not_learned %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", mot_verbs_tri_not_learned, df_count = chi_types),
         seq3 = assign_seq(., "seq3", mot_verbs_tri_not_learned, df_count = chi_types),
         seq4 = assign_seq(., "seq4", mot_verbs_tri_not_learned, df_count = chi_types),
         seq5 = assign_seq(., "seq5", mot_verbs_tri_not_learned, df_count = chi_types)) %>%
  summary_seq_tokens()

# to proportions
summary_seq_types_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_types_not_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_types_mono_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_types_mono_not_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_types_tri_learned_chi_verbs %<>%
  adorn_perc_summary
summary_seq_types_tri_not_learned_chi_verbs %<>%
  adorn_perc_summary

#### matched verbs in input tokens and types ####
summary_freq <- function(df, name_var) {
  df %>%
    group_by(id) %>%
    summarize(!!name_var := mean(freq))
}

summary_freq_sets <- summary_freq(mot_verbs_mono_learned, "Mono_learned_freq") %>%
  cbind(summary_freq(mot_verbs_mono_not_learned, "Mono_not_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_verbs_bis_learned, "Bis_learned_freq")[,2]) %>%
  cbind(summary_freq(mot_verbs_bis_not_learned, "Bis_not_learned_freq")[,2]) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD")))

temp <- tibble(dfs = c("mot_verbs_mono_learned",
                       "mot_verbs_bis_learned"),
               vars = colnames(summary_freq_sets)[-1] %>% .[grepl("not",.)])

matched_dfs_verbs <- list()
for (rows in 1:nrow(temp)) {
  matched_dfs_verbs[[paste(temp$dfs[rows], "match", sep = "_")]] <-
    unique(get(temp$dfs[rows])$id) %>%
    lapply(function(x) {
      equalize_mean(get(temp$dfs[rows]), x, temp$vars[rows])
    }) %>%
    rbindlist()
} ; rm(rows, temp)

# repeat calculation for new dfs
summary_seq_tokens_mono_learned_match_verbs <- matched_dfs_verbs[["mot_verbs_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs_verbs[["mot_verbs_mono_learned_match"]]),
         seq3 = assign_seq(., "seq3", matched_dfs_verbs[["mot_verbs_mono_learned_match"]]),
         seq4 = assign_seq(., "seq4", matched_dfs_verbs[["mot_verbs_mono_learned_match"]]),
         seq5 = assign_seq(., "seq5", matched_dfs_verbs[["mot_verbs_mono_learned_match"]])) %>%
  summary_seq_tokens()

summary_seq_tokens_mono_learned_match_chi_verbs <- matched_dfs_verbs[["mot_verbs_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_tokens),
         seq3 = assign_seq(., "seq3", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_tokens),
         seq4 = assign_seq(., "seq4", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_tokens),
         seq5 = assign_seq(., "seq5", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_tokens)) %>%
  summary_seq_tokens()

summary_seq_types_mono_learned_match_verbs <- matched_dfs_verbs[["mot_verbs_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = mot_types),
         seq3 = assign_seq(., "seq3", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = mot_types),
         seq4 = assign_seq(., "seq4", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = mot_types),
         seq5 = assign_seq(., "seq5", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = mot_types)) %>%
  summary_seq_tokens()


summary_seq_types_mono_learned_match_chi_verbs <- matched_dfs_verbs[["mot_verbs_mono_learned_match"]] %>%
  distinct(id, section) %>%
  mutate(seq2 = assign_seq(., "seq2", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_types),
         seq3 = assign_seq(., "seq3", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_types),
         seq4 = assign_seq(., "seq4", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_types),
         seq5 = assign_seq(., "seq5", matched_dfs_verbs[["mot_verbs_mono_learned_match"]], df_count = chi_types)) %>%
  summary_seq_tokens()

summary_seq_tokens_mono_learned_match_verbs %<>%
  adorn_perc_summary()
summary_seq_tokens_mono_learned_match_chi_verbs %<>%
  adorn_perc_summary()
summary_seq_types_mono_learned_match_verbs %<>%
  adorn_perc_summary()
summary_seq_types_mono_learned_match_chi_verbs %<>%
  adorn_perc_summary()

