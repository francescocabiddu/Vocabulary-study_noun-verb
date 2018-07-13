# load libraries
lib <- c("magrittr", "tidyverse", 
         "data.table", "fastmatch",
         "beepr", "mailR")
lapply(lib, require, character.only = TRUE)
rm(lib)

# load homemade funs
source("homemade_funs.R")

load("~/Dropbox/tirocinio notts/mot_chi study/mot chi_study noun_verb/ws.RData")

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
mot_uni_cat %<>%
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

# redo unique types in types_root
mot_uni_cat %<>%
  (function(x) {
    x$types_root <- setdif_list(x, "types_root")
    x
  })

# assign cat from mot_mor_uni
mot_uni_cat %<>%
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


mot_uni_cat_table <- mot_uni_cat %>% 
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

#### CHILDREN ####