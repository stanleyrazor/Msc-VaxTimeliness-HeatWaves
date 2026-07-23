
detach("package:autoReg", unload = TRUE)
library(autoReg)

mod_fit2model <- function (fit)
{
  if ("coxph" %in% class(fit)) {
    dataname = as.character(fit$call)[3]
    data = eval(parse(text = dataname))
    f = fit$formula
    y = as.character(f)[2]
    temp1 = str_remove_all(y, "Surv\\(|\\)| ")
    temp1 = unlist(strsplit(temp1, ","))
    timevar = temp1[1]
    statusvar = temp1[2]
    xvars = attr(fit$terms, "term.labels")
    xvars
    timevar
    if (!is.na(statusvar)) {
      if (str_detect(statusvar, "==")) {
        statusvar = unlist(strsplit(statusvar, "=="))[1]
      }
      else if (str_detect(statusvar, "!=")) {
        statusvar = unlist(strsplit(statusvar, "!="))[1]
      }
    }
    add = xvars[str_detect(xvars, "strata\\(|cluster\\(|frailty\\(")]
    if (length(add) > 0) {
      xvars = setdiff(xvars, add)
      add = str_remove_all(add, "strata\\(|cluster\\(|frailty\\(|\\)")
      xvars = c(xvars, add)
    }
    if (is.na(statusvar)) {
      myformula = paste0(timevar, "~", paste0(xvars, collapse = "+"))
    }
    else {
      myformula = paste0(timevar, "~", paste0(c(statusvar,
                                                xvars), collapse = "+"))
    }
    myformula
    fit0 = lm(myformula, data = data, weights = wt)
    modelData = fit0$model
  }
  else if ("glmerMod" %in% class(fit)) {
    modelData = fit@frame
    data = modelData
  }
  else if ("glm" %in% class(fit)) {
    y = as.character(fit$formula)[2]
    y
    if (str_detect(y, "==")) {
      dataname = as.character(fit$call)[length(as.character(fit$call))]
      data = eval(parse(text = dataname))
      f = fit$formula
      y = as.character(f)[2]
      y
      xvars = attr(fit$terms, "term.labels")
      xvars
      if (str_detect(y, "==")) {
        temp = unlist(strsplit(y, "=="))[1]
        temp = str_replace_all(temp, " ", "")
        xvars = c(xvars, temp)
      }
      else if (str_detect(y, "!=")) {
        temp = unlist(strsplit(y, "!="))[1]
        temp = str_replace_all(temp, " ", "")
        xvars = c(xvars, temp)
      }
      add = xvars[str_detect(xvars, "strata\\(|cluster\\(|frailty\\(")]
      if (length(add) > 0) {
        xvars = setdiff(xvars, add)
        add = str_remove_all(add, "strata\\(|cluster\\(|frailty\\(|\\)")
        xvars = c(xvars, add)
      }
      myformula = paste0(y, "~", paste0(xvars, collapse = "+"))
      myformula
      fit0 = lm(myformula, data = data, weights = wt)
      modelData = fit0$model
    }
    else {
      dataname = 'model_data' # as.character(fit$call)[4] #* [length(as.character(fit$call))]
      modelData = eval(parse(text = dataname))
    }
    data = modelData
  }
  else {
    dataname = 'model_data'
    modelData = eval(parse(text = dataname))
    #* originally:
    # dataname = as.character(fit$call)[3]
    # modelData = eval(parse(text = dataname))
  }
  df <- modelData %>% autoReg:::restoreData() %>% autoReg:::restoreData2() %>%
    autoReg:::restoreData3()
  df
}

mod_fit2list <- function (fit)
{
  data = mod_fit2model(fit)
  mode = 1
  if ("glm" %in% attr(fit, "class")) {
    mode = 2
    family = fit$family$family
  }
  xvars = attr(fit$terms, "term.labels")
  xno = length(xvars)
  yvar = as.character(attr(fit$terms, "variables"))[2]
  xvars
  yvar
  data
  fitlist = map(xvars, function(x) {
    myformula = paste0(yvar, "~", x)
    if (mode == 1) {
      fit = lm(as.formula(myformula), data = data, weights = wt)
    }
    else if (mode == 2) {
      fit = glm(as.formula(myformula), family = family,
                data = data, weights = wt)
    }
  })
  class(fitlist) = "fitlist"
  fitlist
}

mod_autoReg_sub <- function (fit, threshold = 2, uni = FALSE, multi = TRUE, final = FALSE,
                             imputed = FALSE, keepstats = FALSE, showstats = TRUE, ...)
{
  xvars = attr(fit$terms, "term.labels")
  yvar = as.character(attr(fit$terms, "variables"))[2]
  xvars
  yvar
  data1 = mod_fit2model(fit)
  data1
  mode = 1
  if ("glm" %in% attr(fit, "class")) {
    mode = 2
    family = fit$family$family
  }

  #* FORCING LM BACK TO MODE 1
  if (fit$family$family == 'gaussian') {
    mode <- 1
    family = fit$family$family
  }

  if (uni == FALSE)
    threshold = 1
  if (length(xvars) <= 1) {
    uni = TRUE
    multi = FALSE
    threshold = 1
  }
  result = getSigVars(fit, threshold = threshold, final = final)
  result
  xvars
  others = setdiff(xvars, names(data1))
  others
  xvars = setdiff(xvars, others)
  if (length(xvars) > 0) {
    formula = paste0(yvar, "~", paste0(xvars, collapse = "+"))
    formula
    if (any(str_detect(names(data1), fixed("I(")))) {
      temp = names(data1)[str_detect(names(data1), fixed("I("))]
      data1 <- data1 %>% select(-all_of(temp))
    }
    df = gaze(x = as.formula(formula), data = data1, show.n = keepstats,
              show.p = FALSE)
    df = as.data.frame(df)
  }
  else {
    df = data.frame(name = "", desc = "", id = "")
    df = df[-1, ]
    df
  }
  df
  others
  if (length(others) > 0) {
    for (i in 1:length(others)) {
      name = others[i]
      desc = "others"
      if (str_detect(name, ":")) {
        desc = "interaction"
        temp = getInteraction(name, data = data1)
        temp$n = NULL
      }
      else if (str_detect(name, fixed("I("))) {
        desc = "interpretation"
        temp = data.frame(name = name, desc = desc,
                          id = name)
      }
      df
      temp
      class(df) = "data.frame"
      df = bind_rows(df, temp)
    }
  }
  df
  if (length(which(duplicated(df$id))) > 0) {
    df = df[-which(duplicated(df$id)), ]
  }
  dflist = list()
  dflist[[1]] = df
  stat = ifelse(mode == 1, "Coefficient", "OR")
  if (uni) {
    fit
    fitlist = mod_fit2list(fit)
    if (keepstats) {
      df1 = map_dfr(fitlist, mod_fit2stats) %>% filter(.data$id !=
                                                     "(Intercept)")
      df1$mode = "univariate"
    }
    else {
      df1 = mod_fit2summary(fitlist)
      names(df1)[2] = paste(stat, "(univariable)")
    }
    dflist[["uni"]] = df1
  }
  if (multi) {
    formula2 = paste0(yvar, "~", paste0(result$sigVars,
                                        collapse = "+"))
    if (mode == 1) {
      # fit2 = lm(as.formula(formula2), data = data1, weights = wt)
      fit2 = fit
    }
    else if (mode == 2) {
      #* fit2 = glm(as.formula(formula2), data = data1, family = family, weights = wt)
      fit2 = fit
    }
    if (keepstats) {
      df2 = mod_fit2stats(fit2)
      df2$mode = "multivariate"
    }
    else {
      df2 = mod_fit2summary(fit2)
      names(df2)[2] = paste(stat, "(multivariable)")
    }
    dflist[["multi"]] = df2
    df2
    if (imputed & (!final)) {
      df4 = imputedReg(fit2, data = data1, ...)
      if (keepstats) {
        if (mode == 2) {
          df4 = df4 %>% select(.data$OR, .data$lower,
                               .data$upper, .data$p.value, .data$id, .data$stats) %>%
            mutate(mode = "imputed") %>% rename(p = .data$p.value)
        }
        else {
          df4 = df4 %>% select(.data$Estimate, .data$lower,
                               .data$upper, .data$p.value, .data$id, .data$stats) %>%
            mutate(mode = "imputed") %>% rename(p = .data$p.value)
        }
      }
      else {
        df4 = df4 %>% select(.data$id, .data$stats)
        temp = paste0(ifelse(mode == 1, "Coefficients ",
                             "OR "), "(imputed)")
        names(df4)[2] = temp
      }
      dflist[["imputed"]] = df4
    }
  }
  if (final & (length(result$finalVars) > 0)) {
    formula3 = paste0(yvar, "~", paste0(result$finalVars,
                                        collapse = "+"))
    if (mode == 1) {
      fit3 = lm(as.formula(formula3), data = data1, weights = wt)
    }
    else if (mode == 2) {
      fit3 = glm(as.formula(formula3), data = data1, family = family, weights = wt)
    }
    if (keepstats) {
      df3 = mod_fit2stats(fit3)
      df3$mode = "final"
    }
    else {
      df3 = mod_fit2summary(fit3)
      names(df3)[2] = paste(stat, "(final)")
    }
    dflist[["final"]] = df3
    if (imputed) {
      df4 = imputedReg(fit3, data = data1, ...)
      if (keepstats) {
        df4 = df4 %>% select(.data$OR, .data$lower,
                             .data$upper, .data$p.value, .data$id, .data$stats) %>%
          mutate(mode = "imputed") %>% rename(p = .data$p.value)
      }
      else {
        df4 = df4 %>% select(.data$id, .data$stats)
        temp = paste0(ifelse(mode == 1, "Coefficients ",
                             "OR "), "(imputed)")
        names(df4)[2] = temp
      }
      dflist[["imputed"]] = df4
    }
  }
  if (keepstats) {
    Final = reduce(dflist[-1], bind_rows)
    Final
  }
  else {
    Final = reduce(dflist, left_join, by = "id")
    Final
  }
  class(Final) = c("autoReg", "data.frame")
  Final[is.na(Final)] = ""
  attr(Final, "model") = ifelse(mode == 1, "lm", "glm")
  if (is.null(attr(data1[[yvar]], "label"))) {
    attr(Final, "yvars") = yvar
  }
  else {
    attr(Final, "yvars") = attr(data1[[yvar]], "label")
  }
  Final
}

mod_modelPlot <- function (fit, widths = NULL, change.pointsize = TRUE, show.OR = TRUE,
                           show.ref = TRUE, bw = TRUE, legend.position = "top", ...)
{
  if (is.null(widths)) {
    if (show.OR) {
      widths = c(1.2, 1, 2.2, 3.5)
    }
    else {
      widths = c(1.2, 1, 2.2, 3.5)
    }
  }
  mode = 1
  xname = "Estimate"
  xlabel = "Coefficient (95% CI)"
  if ("glm" %in% class(fit)) {
    mode = 2
    xname = "OR"
    xlabel = "Odds Ratio (95% CI)"
  }
  else if ("coxph" %in% class(fit)) {
    mode = 3
    xname = "HR"
    xlabel = "Hazard Ratio (95% CI)"
  }

  #* forcing gaussian to be mode 1
  if ("glm" %in% class(fit)) {
    if (fit$family$family == 'gaussian') {
      mode = 1
      xname = "Estimate"
      xlabel = "Coefficient (95% CI)"
    }
  }

  xvars = attr(fit$term, "term.labels") #xxx
  xvars
  data = mod_fit2model(fit)
  others = setdiff(xvars, names(data))
  xvars = setdiff(xvars, others)
  if (length(xvars) > 0) {
    myformula = paste0("~", paste0(xvars, collapse = "+"))
    df1 = gaze(as.formula(myformula), data = data, show.n = TRUE)
    df1$id = str_replace_all(df1$id, "`", "")
    df1$name = str_replace_all(df1$name, "`", "")
    plusminus = "±"
    df1$desc[df1$desc == paste0("Mean ", plusminus, " SD")] = ""
  }
  else {
    df1 = data.frame(name = "", desc = "", N = NA, stats = "",
                     n = 1, id = "")
    df1 = df1[-1, ]
    df1
  }
  del = str_detect(others, "strata\\(|cluster\\(|frailty\\(")
  if (any(del))
    others = others[-which(del)]
  if (length(others) > 0) {
    for (i in 1:length(others)) {
      name = others[i]
      desc = "others"
      if (str_detect(name, ":")) {
        temp = getInteraction(name, data = data)
        temp$N = temp$n
        temp$stats = ""
        temp = temp[c(1, 2, 5, 6, 4, 3)]
      }
      else if (str_detect(name, fixed("I("))) {
        desc = "interpretation"
        id = name
        n = nrow(data)
        N = NA
        temp = data.frame(name = name, desc = desc,
                          N = N, stats = "", n = n, id = id)
      }
      class(df1) = "data.frame"
      tempname = setdiff(names(df1), names(temp))
      for (i in seq_along(tempname)) {
        df1[[tempname[i]]] = NULL
      }
      df1
      df1 = rbind(df1, temp)
    }
  }
  df1$no = 1:nrow(df1)
  df1$stats = NULL
  df1
  fit
  df2 = autoReg(fit, keepstats = TRUE, threshold = 2, ...) #xxx

  if ("survreg" %in% class(fit)) {
    names(df2) <- c('id', names(df2)[-1]) #* forcing the survreg objec to have name id as 1st
  }

  df2 = df2 %>% dplyr::filter(.data$id != "(Intercept)")
  df2
  df = dplyr::full_join(df1, df2, by = "id")
  df
  del = which(is.na(df$stats) & (df$desc == ""))
  if (length(del) > 0)
    df = df[-del, ]
  df$mode[is.na(df$stats)] = "Reference"
  df$stats[is.na(df$stats)] = "Reference"
  if (max(nchar(df$desc), na.rm = TRUE) < 1)
    widths[2] = 0
  if (mode == 1) {
    # originally
    # df$Estimate[is.na(df$Estimate)] = 0
    # xintercept = 0

    # now modified to:
    if ("survreg" %in% class(fit)) {

      # rename the coef column to Estimate, and proceed to impute NA with 0's
      names(df)[which(names(df) == 'ETR')] <- 'Estimate'
      df$Estimate[is.na(df$Estimate)] = 1
      df$Estimate <- as.numeric(df$Estimate)
      xintercept = 1

      df$lower <- as.numeric(df$LB)
      df$upper <- as.numeric(df$UB)
    } else {
      df$Estimate[is.na(df$Estimate)] = 1
      xintercept = 0
    }
  }
  else if (mode == 2) {
    df$OR[is.na(df$OR)] = 1
    xintercept = 1
  }
  else if (mode == 3) {
    df$HR[is.na(df$HR)] = 1
    xintercept = 1
  }
  df
  if (!show.ref)
    df = shorten(df)
  no = which(findDup(df[[1]]) + findDup(df[[2]]) == 2)
  no
  df$name[no] = ""
  if ("desc" %in% names(df)) {
    df$desc[no] = ""
  }
  if ("levels" %in% names(df)) {
    df$levels[no] = ""
  }
  dodge = length(setdiff(unique(df$mode), c(NA, "Reference"))) >
    1
  dodge
  limits = setdiff(unique(df$mode), "Reference")
  p <- ggplot(df, aes_string(x = xname)) + geom_vline(xintercept = xintercept,
                                                      color = "grey30", lty = 2)
  if (dodge) {
    p = p + geom_errorbar(aes(y = reorder(.data$id, .data$no),
                              xmin = .data$lower, xmax = .data$upper, color = mode),
                          position = position_dodge(width = 0.4), width = 0.1)
    if (change.pointsize) {
      p = p + geom_point(aes(y = reorder(.data$id, .data$no),
                             fill = mode, size = .data$n), position = position_dodge(width = 0.4),
                         pch = 22)
    }
    else {
      p = p + geom_point(aes(y = reorder(.data$id, .data$no),
                             fill = mode), position = position_dodge(width = 0.4),
                         pch = 22)
    }
  }
  else {
    p = p + geom_errorbar(aes(y = reorder(.data$id, .data$no),
                              xmin = .data$lower, xmax = .data$upper), width = 0.1)
    if (change.pointsize) {
      p = p + geom_point(aes(y = reorder(.data$id, .data$no),
                             size = .data$n), pch = 15)
    }
    else {
      p = p + geom_point(aes(y = reorder(.data$id, .data$no)),
                         pch = 15)
    }
  }
  p
  p = p + labs(y = "", x = xlabel) + scale_y_discrete(limits = rev)
  if (bw) {
    p = p + scale_color_grey(start = 0, end = 0.2) + scale_fill_grey()
  }
  p = p + guides(size = "none") + theme_bw() + theme(axis.title.x = element_text(),
                                                     axis.title.y = element_blank(), axis.text.y = element_blank(),
                                                     axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                                                     legend.position = ifelse(dodge, "top", "none"))
  tab_base <- ggplot(df, aes(y = reorder(.data$id, .data$no))) +
    scale_y_discrete(limits = rev) + ylab(NULL) + xlab(" ") +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(color = "white"), axis.line = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.y = element_blank(), legend.position = "none",
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_blank())
  tab1 <- tab_base + geom_text(aes(x = 1, label = .data$name)) +
    ggtitle("")
  tab1
  tab2 <- tab_base + geom_text(aes(x = 1, label = .data$desc)) +
    ggtitle("")
  tab2
  if (dodge) {
    tab3 <- tab_base + geom_text(aes(x = 1, label = .data$stats,
                                     group = mode), position = position_dodge(width = 0.4)) +
      guides(color = "none")
  }
  else {
    tab3 <- tab_base + geom_text(aes(x = 1, label = .data$stats))
  }
  result = list(tab1 = tab1, tab2 = tab2, tab3 = tab3, p = p,
                widths = widths, legend.position = legend.position)
  class(result) = "modelPlot"
  result
}


mod_gaze <- function (x, ...)
{
  df = as.data.frame(summary(x)$table)
  df$id = rownames(df)
  df1 = as.data.frame(confint(x))
  df1$id = rownames(df1)
  df = left_join(df, df1, by = "id")

  #*naming the CI
  names(df)[(ncol(df)-1):(ncol(df))] = c("lower", "upper") #* was [6:7]

  df <- df %>% dplyr::select(.data$id, everything())
  df$ETR = exp(df$Value)
  df$LB = exp(df$lower)
  df$UB = exp(df$upper)
  if (x$dist == "weibull") {
    df$HR = exp(-df$Value/x$scale[1])
    df$lower = exp(-df$lower/x$scale[1])
    df$upper = exp(-df$upper/x$scale[1])
    if (length(x$scale) > 1) {
      for (i in 2:length(x$scale)) {
        df$HR[df$id == names(x$scale)[i]] = exp(-df$Value[df$id ==
                                                            names(x$scale)[i]]/x$scale[i])
      }
    }
  }
  else if (x$dist == "exponential") {
    df$HR = exp(-df$Value)
    df$lower = exp(-df$lower)
    df$upper = exp(-df$upper)
  }
  else if (x$dist == "loglogistic") {
    df$OR = exp(-df$Value/x$scale)
    df$lower = exp(-df$lower/x$scale)
    df$upper = exp(-df$upper/x$scale)
  }
  df$temp = df$lower
  df$lower = df$upper
  df$upper = df$temp
  df = df[, c(1:5, 8:10, 11, 6:7)]
  df
  attr(df, "call") = gsub(" ", "", paste0(deparse(x$call),
                                          collapse = ""))
  attr(df, "yvars") = attr(attr(x$terms, "dataClasses"), "names")[1]
  attr(df, "model") = "survreg"
  attr(df, "lik") = fit2lik(x)
  attr(df, "summary") = TRUE
  class(df) = c("autoReg", "data.frame")
  myformat(df)
}



mod_fit2stats <- function (fit, method = "default", digits = 2, mode = 1)
{
  cmode = 1
  if ("glm" %in% attr(fit, "class")) {
    cmode = 2
    family = fit$family$family
  }
  else if ("glmerMod" %in% class(fit)) {
    cmode = 3
  }
  else if ("coxph" %in% class(fit)) {
    cmode = 4
  }
  else if ("survreg" %in% class(fit)) {
    cmode = 5
  }

  #* forcing gaussian to be cmode 1
  if (("glm" %in% attr(fit, "class"))) {
    if (fit$family$family == 'gaussian') {
      cmode <- 1
      family <- fit$family$family
    }
  }

  cmode
  fmt = paste0("%.", digits, "f")
  if (cmode == 5) {
    df = mod_gaze(fit)
    df = df[df$LB != "NA", ]
    addp = function(x) {
      result = c()
      for (i in seq_along(x)) {
        if (substr(x[i], 1, 1) == "<") {
          temp = paste0("p", x[i])
        }
        else {
          temp = paste0("p=", x[i])
        }
        result = c(result, temp)
      }
      result
    }
    df$p1 = addp(df$p)
    if (mode == 1) {
      df$stats = paste0(sprintf(fmt, as.numeric(df$ETR)),
                        " (", sprintf(fmt, as.numeric(df$LB)), "-", sprintf(fmt,
                                                                            as.numeric(df$UB)), ", ", df$p1, ")")
      df
    }
    else {
      df$stats = paste0(sprintf(fmt, as.numeric(df$HR)),
                        " (", sprintf(fmt, as.numeric(df$lower)), "-",
                        sprintf(fmt, as.numeric(df$upper)), ", ", df$p1,
                        ")")
      df
    }
  }
  else if (cmode == 4) {
    result = extractHR(fit)
    names(result)[2:3] = c("lower", "upper")
    result$id = rownames(result)
    result$stats = paste0(sprintf(fmt, result$HR), " (",
                          sprintf(fmt, result$lower), "-", sprintf(fmt, result$upper),
                          ", ", p2character2(result$p), ")")
    df = result
    df
  }
  else if (cmode > 1) {
    result = moonBook::extractOR(fit, method = method, digits = digits)
    names(result)[2:3] = c("lower", "upper")
    result$id = rownames(result)
    result$stats = paste0(sprintf(fmt, result$OR), " (",
                          sprintf(fmt, result$lower), "-", sprintf(fmt, result$upper),
                          ", ", p2character2(result$p), ")")
    df = result
    df
  }
  else if (cmode == 1) {
    result = base::cbind(summary(fit)$coefficients, confint(fit))
    temp = round(result[, c(1, 5, 6)], digits)
    id = rownames(result)
    stats = paste0(sprintf(fmt, temp[, 1]), " (", sprintf(fmt,
                                                          temp[, 2]), " to ", sprintf(fmt, temp[, 3]), ", ",
                   p2character2(result[, 4]), ")")
    df = data.frame(id = id, Estimate = result[, 1], lower = result[,
                                                                    5], upper = result[, 6], stats = stats)
  }
  df
}


mod_fit2summary <- function (fit, mode = 1, ...)
{
  if ("fitlist" %in% class(fit)) {
    df = map_dfr(fit, mod_fit2stats, ...)
    df <- df %>% filter(.data$id != "(Intercept)")
  }
  else {
    df = mod_fit2stats(fit, mode = mode, ...)
  }
  if ("survreg" %in% class(fit)) {
    colnames(df)[1] = "id"
    class(df) = "data.frame"
  }
  df %>% select(.data$id, .data$stats)
}


mod_autoRegsurvreg <- function (x, threshold = 2, uni = FALSE, multi = TRUE, final = FALSE,
                                imputed = FALSE, keepstats = FALSE, mode = 1, ...)
{
  if (uni == FALSE)
    threshold = 1
  fit = x
  data = mod_fit2model(fit)
  temp = as.character(fit$call)[2]
  temp = strsplit(gsub(" ", "", temp), "~")
  y = temp[[1]][1]
  y
  temp1 = str_remove_all(y, "Surv\\(|\\)| ")
  temp1 = unlist(strsplit(temp1, ","))
  timevar = temp1[1]
  statusvar = temp1[2]
  xvars = attr(fit$terms, "term.labels")
  xvars
  xvars
  add = xvars[str_detect(xvars, "strata\\(|frailty\\(")]
  add
  if (str_detect(paste0(deparse(fit$call), collapse = ""),
                 "cluster")) {
    temp = paste0(deparse(fit$call), collapse = "")
    temp = unlist(strsplit(temp, "cluster"))[2]
    temp
    add = c(add, paste0("cluster=", str_remove_all(temp,
                                                   "=|\\)| ")))
    add
  }
  myformula = paste0("~", paste0(xvars, collapse = "+"))
  myformula
  mylist = list()
  mylist[[1]] = gaze(as.formula(myformula), data = data, ...)
  mylist[[1]]
  no = 2
  if (uni) {
    df = mysurvregSimple(fit, threshold = threshold, mode = mode)
    df
    if (keepstats) {
      df = df[c(1:5, 13)]
      df$mode = "univariable"
    }
    else {
      df = df[c(1, 13)]
      if (mode == 1) {
        df = rename(df, `ETR (univariable)` = .data$stats)
      }
      else {
        df = rename(df, `HR (univariable)` = .data$stats)
      }
    }
    mylist[[no]] = df
    no = no + 1
  }
  if (multi) {
    fit = survreg2multi(fit, threshold = threshold)
    if (keepstats) {
      df = mod_fit2stats(fit, mode = mode)
      df$mode = "multivariable"
    }
    else {
      df = mod_fit2summary(fit, mode = mode) #* modified to add id
      if (mode == 1) {
        df = rename(df, `ETR (multivariable)` = .data$stats)
      }
      else {
        df = rename(df, `HR (multivariable)` = .data$stats)
      }
    }
    mylist[[no]] = df
    no = no + 1
  }
  if (final) {
    final1 = survreg2final(fit, threshold = threshold)
    if (keepstats) {
      df = mod_fit2stats(final1, mode = mode)
      df$mode = "final"
    }
    else {
      df = mod_fit2summary(final1, mode = mode)
      if (mode == 1) {
        df = rename(df, `ETR (final)` = .data$stats)
      }
      else {
        df = rename(df, `HR (final)` = .data$stats)
      }
    }
    mylist[[no]] = df
    no = no + 1
  }
  if (imputed) {
    imputed = imputedReg(fit, mode = mode)
    if (keepstats) {
      df = imputed[c(ifelse(mode == 1, "ETR", ifelse(fit$dist ==
                                                       "loglogistic", "OR", "HR")), "lower", "upper",
                     "p.value", "id", "stats")] %>% rename(p = .data$p.value)
      df$mode = "imputed"
    }
    else {
      df = imputed[c("id", "stats")]
      if (mode == 1) {
        df = rename(df, `ETR (imputed)` = .data$stats)
      }
      else if (fit$dist != "loglogistic") {
        df = rename(df, `HR (imputed)` = .data$stats)
      }
      else {
        df = rename(df, `OR (imputed)` = .data$stats)
      }
    }
    mylist[[no]] = df
  }
  if (keepstats) {
    Final = reduce(mylist[-1], bind_rows)
    Final
  }
  else {
    Final = reduce(mylist, left_join, by = "id")
    Final
  }
  class(Final) = c("autoReg", "data.frame")
  Final[is.na(Final)] = ""
  if (length(add) > 0) {
    attr(Final, "add") = add
    if (str_detect(add, "frailty")) {
      if (final) {
        no = which(str_detect(rownames(summary(final1)$coefficients),
                              "frailty"))
        p = data.frame(summary(final1)$coef)$p[no]
        attr(Final, "add") = paste(add, p2character2(p,
                                                     add.p = TRUE), summary(final1)$print2)
      }
      else {
        no = which(str_detect(rownames(summary(fit)$coefficients),
                              "frailty"))
        p = data.frame(summary(fit)$coef)$p[no]
        attr(Final, "add") = paste(add, p2character2(p,
                                                     add.p = TRUE), summary(fit)$print2)
      }
    }
  }
  attr(Final, "yvars") = attr(attr(fit$terms, "dataClasses"),
                              "names")[1]
  attr(Final, "model") = "survreg"
  temp = summary(fit)$logtest
  Final
}



# -------------------------------------------------------------------------

# Unlock the namespace environment
unlockBinding("fit2list", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("fit2list", mod_fit2list, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("fit2list", asNamespace("autoReg"))



# Unlock the namespace environment
unlockBinding("fit2model", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("fit2model", mod_fit2model, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("fit2model", asNamespace("autoReg"))



# Unlock the namespace environment
unlockBinding("autoReg_sub", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("autoReg_sub", mod_autoReg_sub, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("autoReg_sub", asNamespace("autoReg"))



# Unlock the namespace environment
unlockBinding("modelPlot", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("modelPlot", mod_modelPlot, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("modelPlot", asNamespace("autoReg"))



# Unlock the namespace environment
unlockBinding("fit2summary", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("fit2summary", mod_fit2summary, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("fit2summary", asNamespace("autoReg"))




# Unlock the namespace environment
unlockBinding("fit2stats", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("fit2stats", mod_fit2stats, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("fit2stats", asNamespace("autoReg"))




# Unlock the namespace environment
unlockBinding("gaze", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("gaze", mod_gaze, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("gaze", asNamespace("autoReg"))



# Unlock the namespace environment
unlockBinding("autoRegsurvreg", asNamespace("autoReg"))

# Overwrite the internal a3 function with your modified one
assignInNamespace("autoRegsurvreg", mod_autoRegsurvreg, ns = "autoReg")

# Re-lock it to keep it stable
lockBinding("autoRegsurvreg", asNamespace("autoReg"))

