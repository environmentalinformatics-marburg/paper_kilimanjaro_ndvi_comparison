mkStats <- function(x, digits = 3) {
  
  # trends
  num_mk <- x[]
  num_mk_abs <- sum(!is.na(num_mk))
  num_mk_rel <- round(num_mk_abs / ncell(x), digits)
  
  # positive trends
  num_mk_abs_pos <- sum(num_mk > 0, na.rm = TRUE)
  num_mk_rel_pos <- round(num_mk_abs_pos / num_mk_abs, digits)
  
  # negative trends
  num_mk_abs_neg <- sum(num_mk < 0, na.rm = TRUE)
  num_mk_rel_neg <- round(num_mk_abs_neg / num_mk_abs, digits)
  
  # return data
  data.frame(trends_abs = num_mk_abs, trends_rel = num_mk_rel, 
             trends_abs_pos = num_mk_abs_pos, trends_rel_pos = num_mk_rel_pos, 
             trends_abs_neg = num_mk_abs_neg, trends_rel_neg = num_mk_rel_neg)
}