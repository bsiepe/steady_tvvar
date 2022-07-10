# Description -------------------------------------------------------------
# This script contains auxiliary functions for multiple scripts


# floorceiling ------------------------------------------------------------
# simple function for showing floor/ceiling effects of variables
# by computing relative frequency
floorceiling <- function(df, var){
  df %>% 
    group_by(id) %>% 
    # count occurence of distinct values
    count({{var}}) %>% 
    # relative frequency of individual values
    mutate(rel_freq = n/sum(n)) %>% 
    # look for the maximum or minimum value
    filter({{var}} == min({{var}}) | {{var}} == max({{var}})) %>% 
    ungroup() %>% 
    # sort by size of relative frequency
    arrange(desc(rel_freq))
}


# search_na_gap -----------------------------------------------------------
### Write a function that searches for consec. number of missings per item
# and is set to TRUE when consec missings > maxgap
# col is the variable name, gap is the highest number of consec. missings allowed
# this is adapted from the "maxgap" option from imputeTS::na_kalman 
# see https://rdrr.io/github/SteffenMoritz/imputeTS/src/R/na_kalman.R
# Set gap to 7 per default according to preregistration
search_na_gap <- function(col, gap = 7){
  # compute run of NAs
  encoding <- rle(is.na(col))
  
  # set runs lower than the maximum gap to FALSE (they can be included)
  # FALSE means that this row does not need to be deleted
  encoding$values[encoding$lengths <= gap] <- FALSE
  
  # Use reverse.cls to indicate for each row if it's FALSE (to be included)
  # or TRUE (to be removed). TRUE means that the length of missingness run > gap
  en <- inverse.rle(encoding)
  return(en)
} 

# niceTable_moredigits ----------------------------------------------------
# This function is used in script 3. Analysis
# it simply extends the already used niceTable function to show more digits
# https://raw.githubusercontent.com/RemPsyc/niceplots/master/niceTableFunction.R
# this is relevant for our results table 
niceTable_moredigits <- function (dataframe, italics = NULL, highlight = FALSE) {
  if(!require(flextable)){install.packages("flextable")}
  if(!require(dplyr)){install.packages("dplyr")}
  library(flextable)
  library(dplyr)
  if("CI_lower" %in% names(dataframe) & "CI_upper" %in% names(dataframe)) {
    dataframe[,c("CI_lower", "CI_upper")] <- lapply(lapply(dataframe[,c("CI_lower", "CI_upper")], as.numeric), round, 2)
    dataframe["95% CI"] <- apply(dataframe[,c("CI_lower", "CI_upper")], 1, function(x) paste0("[", x[1], ", ", x[2], "]"))
    dataframe <- select(dataframe, -c("CI_lower", "CI_upper"))
  }
  if(highlight == TRUE) {
    dataframe %>%
      mutate(signif = ifelse(p < .05, TRUE, FALSE)) -> dataframe
  }
  nice.borders <- list("width" = 0.5, color = "black", style = "solid")
  dataframe %>%
    {if(highlight == TRUE) flextable(., col_keys = names(dataframe)[-length(dataframe)]) 
      else flextable(.)} %>%
    theme_booktabs %>%
    hline_top(part="head", border = nice.borders) %>%
    hline_bottom(part="head", border = nice.borders) %>%
    hline_top(part="body", border = nice.borders) %>%
    hline_bottom(part="body", border = nice.borders) %>%
    fontsize(part = "all", size = 12) %>%
    font(part = "all", fontname = "Times New Roman") %>%
    align(align = "center", part = "all") %>%
    #line_spacing(space = 2, part = "all") %>%
    height(height = 0.55, part = "body") %>%
    height(height = 0.55, part = "head") %>%
    hrule(rule = "exact", part = "all") %>%
    set_table_properties(layout = "autofit") -> table
  if(!missing(italics)) {
    table %>%
      italic(j = italics, part = "header") -> table
  }
  format.p <- function(p, precision = 0.001) {
    digits <- -log(precision, base = 10)
    p <- formatC(p, format = 'f', digits = digits)
    p[p == formatC(0, format = 'f', digits = digits)] <- paste0('< ', precision)
    sub("0", "", p)
  }
  format.r <- function(r, precision = 0.01) {
    digits <- -log(precision, base = 10)
    r <- formatC(r, format = 'f', digits = digits)
    sub("0", "", r)}
  if("p" %in% names(dataframe)) {
    table %>%
      italic(j = "p", part = "header") %>%
      set_formatter(p = function(x) {
        format.p(x)}) -> table
  }  
  if("r" %in% names(dataframe)) {
    table %>%
      italic(j = "r", part = "header") %>%
      set_formatter(r = function(x)
        format.r(x)) -> table
  }
  if("t" %in% names(dataframe)) {
    table %>%
      italic(j = "t", part = "header") %>%
      colformat_double(j = "t", big.mark=",", digits = 2) -> table
  }
  if("SE" %in% names(dataframe)) {
    table %>%
      italic(j = "SE", part = "header") %>%
      colformat_double(j = "SE", big.mark=",", digits = 2) -> table
  }
  if("SD" %in% names(dataframe)) {
    table %>%
      italic(j = "SD", part = "header") %>%
      colformat_double(j = "SD", big.mark=",", digits = 2) -> table
  }
  if("F" %in% names(dataframe)) {
    table %>%
      italic(j = "F", part = "header") %>%
      colformat_double(j = "F", big.mark=",", digits = 2) -> table
  }
  if("df" %in% names(dataframe)) {
    table %>%
      italic(j = "df", part = "header") %>%
      colformat_double(j = "df", big.mark=",", digits = 0) -> table
  }
  if("b" %in% names(dataframe)) {
    table %>%
      italic(j = "b", part = "header") %>%
      colformat_double(j = "b", big.mark=",", digits = 2) -> table
  }
  if("M" %in% names(dataframe)) {
    table %>%
      italic(j = "M", part = "header") %>%
      colformat_double(j = "M", big.mark=",", digits = 2) -> table
  }
  if("B" %in% names(dataframe)) {
    table %>%
      compose(i = 1, j = "B", part = "header",
              value = as_paragraph("ß")) %>%
      colformat_double(j = "B", big.mark=",", digits = 2) -> table
  }
  if("R2" %in% names(dataframe)) {
    table %>%
      compose(i = 1, j = "R2", part = "header",
              value = as_paragraph(as_i("R"), as_sup("2"))) %>%
      set_formatter(R2 = function(x)
        format.r(x)) -> table
  }
  if("sr2" %in% names(dataframe)) {
    table %>%
      compose(i = 1, j = "sr2", part = "header",
              value = as_paragraph(as_i("sr"), as_sup("2"))) %>%
      set_formatter(sr2 = function(x)
        format.r(x)) -> table
  }
  if("np2" %in% names(dataframe)) {
    table %>%
      compose(i = 1, j = "np2", part = "header",
              value = as_paragraph("<U+03B7>", as_sub("p"), as_sup("2"))) %>%
      colformat_double(j = "np2", big.mark=",", digits = 2) -> table
  }
  if("ges" %in% names(dataframe)) {
    table %>%
      compose(i = 1, j = "ges", part = "header",
              value = as_paragraph("<U+03B7>", as_sub("G"), as_sup("2"))) %>%
      colformat_double(j = "ges", big.mark=",", digits = 2) -> table
  }
  if("dR" %in% names(dataframe)) {
    table %>%
      compose(i = 1, j = "dR", part = "header",
              value = as_paragraph(as_i("d"), as_sub("R"))) %>%
      colformat_double(j = "dR", big.mark=",", digits = 2) -> table
  }
  if("d" %in% names(dataframe)) {
    table %>%
      italic(j = "d", part = "header") %>%
      colformat_double(j = "d", big.mark=",", digits = 2) -> table
  }
  if(highlight == TRUE) {
    table %>%
      bold(i = ~ signif == TRUE,
           j = table$col_keys) %>%
      bg(i = ~ signif == TRUE,
         j = table$col_keys,
         bg = "#D9D9D9") -> table
  }
  table %>%
    colformat_double(j = (select(dataframe, where(is.numeric)) %>%
                            select(-matches("^p$|^r$|^t$|^SE$|^SD$|^F$|^df$|
                                    ^b$|^M$|^B$|^R2$|^sr2$|^np2$|^dR$",
                                            ignore.case =F)) %>% names), 
                     # the following line was changed to accomodate more digits
                     big.mark=",", digits = 4) -> table
  table
}

# tvpar_ggplot ------------------------------------------------------------
# Function for plotting single parameters over time
# Idea: It might be helpful to be able to plot individual parameters
# investigate individual networks better
# Call this function tvpar_ggplot
tvpar_ggplot <- function(object, id, dv, iv,
                         ci = FALSE, ci_object){
  # object = mgm::tvmvar list object where results are stored
  # id = numerical id
  # dv = dependent variable
  # iv = independent variable
  # ci = TRUE implies that uncertainty bands are plotted
  # ci_object is the object where bootstrapped results are stored
  
  # check if ci_object was specified
  if(missing(ci_object)) ci_object <- NULL 
  
  
  # Before plotting: get proper sign again
  # Extract parameters
  ests <- object[[id]]$wadj
  # Indicator for negative parameters
  ind_neg <- which(object[[i]]$signs == -1, arr.ind = T)
  # give parameters their correct negative sign
  ests[ind_neg] <- ests[ind_neg] * -1   # TODO: Something here does not work, gives wrong sign!
  # only use relevant variables
  ests <- ests[dv, iv,1,]
  # Plot without CIs
  if(ci == FALSE){
    print(ggplot()+
            geom_line(aes(x = seq(1, length(ests), by = 1),
                          y = ests))+
            theme_light()+
            labs(x = "Estimation Point",
                 y = "Estimate")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
    print(ests)
  }
  
  # Add confidence intervals
  if(ci == TRUE){
    cints <- apply(ci_object[[id]]$bootParameters[dv, iv, 1, , ], 1, 
                   function(x) {quantile(x, probs = c(.05, .95))})
    low <- cints[1,]
    up <- cints[2,]
    plot_data <- as.data.frame(cbind(ests,low,up))
    print(ggplot(plot_data, aes (x = seq(1, nrow(plot_data), by = 1)))+
            geom_line(aes(y = ests))+
            geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.2)+
            theme_light()+
            labs(x = "Estimation Point",
                 y = "Estimate"))
  }
}
##  Example:
# tvpar_ggplot(object = tvvar_res, id = 6,
             # dv = 4, iv = 4,
             # ci = FALSE)

