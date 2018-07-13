round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

notify_me <- function(subj = "R: Come back to me", pass, recipient) {
  send.mail(from = "fcabiddur@gmail.com",
            to = recipient,
            subject= subj,
            body = "DO NOT REPLY",
            smtp = list(host.name = "smtp.gmail.com", port = 465, 
                        user.name="fcabiddur@gmail.com", passwd=pass, ssl=TRUE),
            authenticate = TRUE,
            send = TRUE)
}

setdif_list <- function(x, var_name) {
  # consecutive set differences between list elements
  uni_diff <- c()
  for (i in unique(x$baby)) {
    sub <- subset(x, baby == i)
    learnt_voc <- sub[[var_name]][[1]]
    for (i in 2:length(sub[[var_name]])) {
      sub[[var_name]][[i]] <- setdiff(sub[[var_name]][[i]], learnt_voc)
      learnt_voc <- c(learnt_voc,  sub[[var_name]][[i]]) 
    }
    uni_diff <- c(uni_diff, sub[[var_name]])
  }
  
  return(uni_diff)
} 
