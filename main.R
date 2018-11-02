# load libraries
lib <- c("magrittr", "tidyverse", 
         "data.table", "fastmatch",
         "beepr", "mailR")
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

#### Mot bisyllabic nouns learned ####
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
                          cat = (chi_uni_cat %>% filter(baby == y))$cat %>% unlist()) %>%
        filter(cat == "N"))$word
      
      mot_nouns_bis %>%
        filter(id == y) %>%
        filter(word %in% reference)
    }) %>% rbindlist()
  })

mot_nouns_bis_not_learned <- mot_nouns_bis %>%
  (function(x) {
    lapply(unique(x$id), function(y) {
      reference <- (tibble(word = (chi_uni_cat %>% filter(baby == y))$types_root %>% unlist(), 
                           cat = (chi_uni_cat %>% filter(baby == y))$cat %>% unlist()) %>%
                      filter(cat == "N"))$word
      
      mot_nouns_bis %>%
        filter(id == y) %>%
        filter(!word %in% reference)
    }) %>% rbindlist()
  })

# assign to each bis how many mono contains (grouped by id)
mot_nouns_bis_learned %<>%
  mutate(monosyll = sapply(unique(id), function(name) {
  sapply((mot_nouns_bis_learned %>% filter(id == name))$phon, function(x) {
    x %>%
      str_detect((mot_mono %>% filter(id == name))$phon) %>%
      sum()
    })
    }) %>%
    unlist())

mot_nouns_bis_not_learned %<>%
  mutate(monosyll = sapply(unique(id), function(name) {
    sapply((mot_nouns_bis_not_learned %>% filter(id == name))$phon, function(x) {
      x %>%
        str_detect((mot_mono %>% filter(id == name))$phon) %>%
        sum()
    })
  }) %>%
    unlist())

# sum number of mono in (1a) and (1b) and take average by id
summary_bysillables <- mot_nouns_bis_learned %>%
  mutate(monosyll = monosyll %>% (function(x) {
    names(x) <- NULL
    x
  })) %>%
  group_by(id) %>%
  summarize(tot_mono_learned = mean(monosyll)) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
  bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD"))) %>%
  cbind(tot_mono_not_learned = mot_nouns_bis_not_learned %>%
          mutate(monosyll = monosyll %>% (function(x) {
            names(x) <- NULL
            x
          })) %>%
          group_by(id) %>%
          summarize(tot_mono_not_learned = mean(monosyll)) %>%
          bind_rows(summarise_all(., funs(if(is.numeric(.)) mean(.) else "MEAN"))) %>%
          bind_rows(summarise_all(., funs(if(is.numeric(.)) sd(.) else "SD"))) %>%
          select(tot_mono_not_learned))
