library(readstata13)
library(ggplot2)
library(dplyr)
library(magrittr)
library(boot)
library(reshape) 
library(nnet)
library(tableone)

CIlo <- function(x, p = 0.025) quantile(x, probs = p, na.rm = TRUE) 
CIhi <- function(x, p = 0.975) quantile(x, probs = p, na.rm = TRUE) 

df <- read.dta13("oai_oarsi_weights.dta") 
colnames(df)[match("_20mtim", colnames(df))] <- "time20m"
df <- select(df, age, ed4, race2, sex, hspss, time20m, pyc, bmi, hsmss, comorb) %>% na.omit()
df <-  mutate(df, 
    race2 = factor(race2, labels = c("White", "Black")),
    ed4 = factor(ed4, labels = c("HS or less", "Some College", "College", "Grad School")),
    obcat = cut(bmi, breaks = quantile(df$bmi, probs=seq(0,1,by = 0.2))),
    pyc = factor(pyc),
    depcat = cut(hsmss, breaks = quantile(df$hsmss, probs=seq(0,1,by = 0.2))),
    agevar = age
  )
df$ed4 <- relevel(df$ed4, ref = levels(df$ed4)[4])
ref <- 1
df[df$comorb > 1,]$comorb = 2
df$comorb <- factor(df$comorb)

trace
S <- 300; len <- length(levels(df$ed4))

sesvar <- "ed4"
confounders <- c("agevar", "sex", "race2")
mediators <- c("obcat", "pyc", "depcat") 
outcomes <- c("hspss", "time20m")
res_b <- array(NA, dim=c(S,length(outcomes), length(mediators)+2, len), 
             dimnames = list(1:S, outcomes, c("Unadj", "demog", "rm: obesity", "rm: smk", "rm: depcat"), 
             c("Intercept", levels(df$ed4)[2:len])))
res_b2 <- array(NA, dim=c(S,length(outcomes), length(mediators)+2, 2), 
               dimnames = list(1:S, outcomes, c("Unadj", "demog", "rm: obesity", "rm: smk", "rm: depcat"))) 
res_a <- array(NA, dim=c(S,length(outcomes), len-1, 2), dimnames = list(1:S, outcomes, levels(df$ed4)[2:len], c("Absolute", "Ratio")))
res_c <- array(NA, dim=c(S,length(outcomes), len-1, 2), dimnames = list(1:S, outcomes, levels(df$ed4)[2:len], c("Absolute", "Ratio")))

iptw_mediat_num <- lapply(mediators, function (i) as.formula(paste(i, "~ ed4")))
med_formulas <- list(
  obcat = paste(sesvar,"pyc", paste(confounders, collapse = "+"), sep = "+"),
  pyc = paste(sesvar,"depcat", paste(confounders, collapse = "+"), sep = "+"),
  depcat = paste(sesvar, paste(confounders, collapse = "+"), sep = "+")
)
iptw_mediat_den <- lapply(mediators, function (i) as.formula(paste(i, med_formulas[[i]], sep = "~")))

pb <- txtProgressBar(style=3)

for(i in 1:S){
  df_temp <- sample_frac(df, 1, replace = TRUE)
  invisible({
    ### CALCULATES WEIGHTS FOR TOTAL EFFECT
    wa_nums <- multinom(ed4 ~ 1, data = df_temp, trace = FALSE) %$% predict(object = ., newdata=df_temp, "probs")
    wa_dens <- multinom(ed4 ~ sex + race2 + agevar, data = df_temp, trace = FALSE) %$% predict(object = ., newdata=df_temp, "probs")
    df_temp$wa <- (wa_nums / wa_dens)[cbind(1:nrow(df_temp),as.numeric(df_temp$ed4))]
    
    #### CALCULATES OVERALL MEDIATED EFFECT (GRAPH A)
    df_temp$wa_une <- (wa_nums / wa_dens)[cbind(1:nrow(df_temp),ref)]
    for(j in 1:length(levels(df_temp$ed4))) df_temp[[paste0("wa_exp_", j)]] <- (wa_nums / wa_dens)[,j]
    for(j in outcomes){
      q1 <- with(df_temp[df_temp$ed4 == levels(df_temp$ed4)[ref],], weighted.mean(get(j), wa_une) )

      q2 <- do.call(rbind,lapply( levels(df_temp$ed4)[-ref], 
                                  function (k) with( df_temp[df_temp$ed4 == k,], weighted.mean(get(j), get(paste0("wa_exp_", match(k, levels(df_temp$ed4))))) ) ))

      tmpred <- do.call(cbind, lapply(levels(df_temp$ed4)[-ref], function (k)
        lm(get(j) ~ sex + agevar + race2 + ed4 + pyc + obcat + depcat, data = df_temp) %$% predict(object = ., newdata = df_temp %>% mutate(ed4 = k ))
      ))
      q3 <- do.call(rbind, lapply(levels(df_temp$ed4)[-ref], function (k) 
                    weighted.mean(tmpred[ df_temp$ed4 == levels(df_temp$ed4)[ref], match(k, levels(df_temp$ed4)) - (ref<k)], 
                                  df_temp[df_temp$ed4 == levels(df_temp$ed4)[ref],]$wa_une, na.rm=TRUE ) )) 
      res_a[i,match(j,outcomes),,1] <- q2-q1+q2-q3
      res_a[i,match(j,outcomes),,2] <- log(q2-q3) - log(q2-q1+q2-q3)
    }
    #### CALCULATES EFFECT MEDIATED BY VARIABLES SEPARATELY (GRAPH B)
    for(j in 1:length(mediators)){
      df_temp[[paste0("wm_num_", mediators[j] )]] <- multinom(iptw_mediat_num[[j]], data=df_temp, trace = FALSE) %$% 
        predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df_temp),as.numeric(df_temp[[mediators[j]]]))]
      df_temp[[paste0("wm_den_", mediators[j] )]] <- multinom(iptw_mediat_den[[j]], data=df_temp, trace = FALSE) %$% 
        predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df_temp),as.numeric(df_temp[[mediators[j]]]))]
      df_temp[[paste0("w_",mediators[j])]] <- df_temp[[paste0("wm_num_", mediators[j] )]] / df_temp[[paste0("wm_den_", mediators[j] )]] * df_temp$wa
    }
    for(j in 1:length(outcomes)){
      res_b[i,j,1,] <- lm(data = df_temp, get(outcomes[j]) ~ ed4)$coefficients[1:len]
      res_b[i,j,2,] <- lm(data = df_temp, get(outcomes[j]) ~ ed4, weight = wa)$coefficients[1:len]
      for(k in 1:length(mediators)){
        tmw <- paste0("w_",mediators[k])
        res_b[i,j,k+2,] <- lm(data = df_temp, get(outcomes[j]) ~ ed4 + get(mediators[k]), weight = get(tmw))$coefficients[1:len]
        res_b2[i,j,k+2,1] <- res_b[i,j,2,2] - res_b[i,j,k+2,2] 
        res_b2[i,j,k+2,2] <- log(res_b[i,j,2,2] - res_b[i,j,k+2,2]) - log(res_b[i,j,2,2])
      }
    }
    #### CALCULATES EFFECT MEDIATED BY COMORBIDITIES KEEPING SMK, BMI UNCHANGED (GRAPH C)
    df_temp$wm_num_comorb <- multinom(comorb ~ ed4, data = df_temp, trace = FALSE) %$% 
      predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df_temp),as.numeric(df_temp$comorb))]
    df_temp$wm_den_comorb <- multinom(comorb ~ ed4 + agevar + sex + race2 + pyc + obcat, data = df_temp, trace = FALSE) %$% 
      predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df_temp),as.numeric(df_temp$comorb))]
    df_temp$w_comorb <- df_temp$wm_num_comorb / df_temp$wm_den_comorb * df_temp$wa
    for(j in 1:length(outcomes)){
      tm1 <- lm(data = df_temp, get(outcomes[j]) ~ ed4, weight = wa)$coefficients[2:len]
      tm2 <- lm(data = df_temp, get(outcomes[j]) ~ ed4 + comorb, weight = w_comorb)$coefficients[2:len]
      res_c[i,j,,1] <- tm1 - tm2
      res_c[i,j,,2] <- log(tm1 - tm2) - log(tm1)
    }
  })
  setTxtProgressBar(pb, i/S)
}
resdiff_a <- res_a
restot <- melt(resdiff_a) %>%
  group_by(X2, X3, X4) %>% 
  summarise(mean = mean(value, na.rm = TRUE), 
            l95 = CIlo(value, 0.025), 
            u95 = CIhi(value, 0.975),
            l90 = CIlo(value, 0.05),
            u90 = CIhi(value, 0.95)
  )
for(i in c("mean", "l95", "u95")) restot[restot$X4=="Ratio",][[i]] <- exp(restot[restot$X4=="Ratio",][[i]]) 

resdiff_b <- res_b
for(i in 3:5) resdiff_b[,,i,] <-  resdiff_b[,,2,] - resdiff_b[,,i,]  
resdat <- melt(resdiff_b[,,3:5,]) %>%  #cbind(.,melt(resdiff_b[,,2,])) %>% #[,,3:5,]
  group_by(X2, X3, X4) %>% 
  summarise(mean = mean(value, na.rm = TRUE), 
            l95 = CIlo(value, 0.025), 
            u95 = CIhi(value, 0.975),
            l90 = CIlo(value, 0.05),
            u90 = CIhi(value, 0.95)
  )
for(i in c("mean", "l95", "u95")) resdat[resdat$X4==2,][[i]] <- exp(resdat[resdat$X4==2,][[i]]) 

resdiff_b <- res_b2
resdat <- melt(resdiff_b[,,3:5,]) %>% #[,,3:5,]) %>%  #cbind(.,melt(resdiff_b[,,2,])) %>% #[,,3:5,]
  group_by(X2, X3, X4) %>% 
  summarise(mean = mean(value, na.rm = TRUE), 
            l95 = CIlo(value, 0.025), 
            u95 = CIhi(value, 0.975),
            l90 = CIlo(value, 0.05),
            u90 = CIhi(value, 0.95)
  )
for(i in c("mean", "l95", "u95")) resdat[resdat$X4==2,][[i]] <- exp(resdat[resdat$X4==2,][[i]]) 

resdiff_c <- res_c
rescomb <- melt(resdiff_c) %>%
  group_by(X2, X3, X4) %>% 
  summarise(mean = mean(value, na.rm = TRUE), 
            l95 = CIlo(value, 0.025), 
            u95 = CIhi(value, 0.975),
            l90 = CIlo(value, 0.05),
            u90 = CIhi(value, 0.95)
  )
for(i in c("mean", "l95", "u95")) rescomb[rescomb$X4=="Ratio",][[i]] <- exp(rescomb[rescomb$X4=="Ratio",][[i]]) 

#       q2 <- do.call(rbind,lapply( 2:length(unique(df_temp$ed4)), 
#                                   function (k) with( df_temp[df_temp$ed4 == levels(df_temp$ed4)[k],], weighted.mean(get(j), get(paste0("wa_exp_", k))) ) ))

#       tmpred <- do.call(cbind, lapply(2:length(unique(df_temp$ed4)), function (k)
#         lm(get(j) ~ sex + agevar + race2 + ed4 + pyc + obcat + depcat, data = df_temp) %$% predict(object = ., newdata = df_temp %>% mutate(ed4 = levels(ed4)[k]) )
#       ))

# 
# filter(resdat, X2 == "hspss") %>% ggplot(aes(x = factor(X3), fill = X4)) + 
#   geom_boxplot(aes(ymin=l95, lower=l_75/e_mean, ymax=m_95/e_mean, upper = m_75/e_mean, middle=p_median/e_mean), stat="identity")
# 
# 
# gdat <- melt(data = df, id=c("id", "sex", "race2", "ed3", "age3", "hspss", "wom_dis", "time20m", "obcat", "pyc", "depressed")) 
# colnames(gdat)[match("value", colnames(gdat))] <- "weight_val"
# res_hspss <- group_by(gdat, variable, ed3) %>% summarize(res = weighted.mean(hspss,weight_val)) 
# gdat %>% ggplot(aes(x = factor(variable), y = hspss, fill = ed3, weight = weight_val)) + geom_point()

# 
# for(i in 1:S)
# {
#   for(j in 1:length(outcomes)){
#     res[i,j,1,] <- sample_frac(df, 1, replace = TRUE) %$% lm(data =., get(outcomes[j]) ~ ed3)$coefficients[1:3]
#     res[i,j,2,] <- sample_frac(df, 1, replace = TRUE, weight = wa_int) %$% lm(data =., get(outcomes[j]) ~ ed3 )$coefficients[1:3]
#     for(k in 1:length(mediators)){
#       tmw <- paste0("w_noint_",mediators[k])
#       res[i,j,k+2,] <- sample_frac(df, 1, replace = TRUE, weight = get(tmw)) %$% lm( data =., get(outcomes[j]) ~ ed3 + get(mediators[k]) )$coefficients[1:3]
#     }
#   }
#   setTxtProgressBar(pb, i/S)
# }

#   
#   wm_obcat_num <- multinom(obcat ~ ed4, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df),as.numeric(df$obcat))]
#   wm_obcat_den <- multinom(obcat ~ ed4 + sex + age3 + race2 + pyc + depressed, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df),as.numeric(df$obcat))]
#   df_temp$w_obcat = wm_obcat_num/wm_obcat_den*df_temp$wa
#   
#   wm_pyc_num <- multinom(pyc ~ ed4, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df),as.numeric(df$pyc))]
#   wm_pyc_den <- multinom(pyc ~ ed4 + sex + age3 + race2 + depressed, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df),as.numeric(df$pyc))]
#   df_temp$w_pyc = wm_pyc_num/wm_pyc_den*df_temp$wa
#   
#   wm_depressed_num <- multinom(pyc ~ ed4, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df),as.numeric(df$depressed))]
#   wm_depressed_den <- multinom(pyc ~ ed4 + sex + age3 + race2, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")[cbind(1:nrow(df),as.numeric(df$depressed))]
#   df_temp$w_depressed = wm_depressed_num/wm_depressed_den*df_temp$wa
# 
# for(j in 2:length(unique((df_temp$ed4)))){
#   df_temp$wa_exp <- (wa_nums / wa_dens)[cbind(1:nrow(df_temp),j)]
#   
#   for(j in outcomes) df_temp$q <- lm(hspss ~ sex + age3 + race2 + ed4 + pyc + obcat + depressed, data = df_temp) %$% predict(object = ., newdata = df_temp)
#   res_a
# }
# 
# wa_nums <- multinom(ed4 ~ 1, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")
# wa_dens <- multinom(ed4 ~ sex + race2 + age3, data = df_temp) %$% predict(object = ., newdata=df_temp, "probs")
# df_temp$wa <- (wa_nums / wa_dens)[cbind(1:nrow(df_temp),as.numeric(df_temp$ed4))]
# 
# q1 <- with(df_temp[df_temp$ed4 == levels(df_temp$ed4)[1],], weighted.mean(hspss, wa_une) )
# q2 <- lapply(2:length(unique(df_temp$ed4)), function (k) with(df_temp[df_temp$ed4 == levels(df_temp$ed4[k])], weighted.mean(hspss, wa_exp)) )

# for(j in outcomes){
#   df_temp$q <- lm(hspss ~ sex + age3 + race2 + ed4 + pyc + obcat + depressed, data = df_temp) %$% predict(object = ., newdata = df_temp)
# }
# df_temp$q <- lm(hspss ~ sex + age3 + race2 + ed4 + pyc + obcat + depressed, data = df_temp) %$% predict(object = ., newdata = df_temp)
# 
# df_temp$q1 <- lm(hspss ~ ed4, data = df_temp) %$% predict(object = ., newdata = df_temp %>% mutate(ed4 = levels(ed4)[1]), weight = df_temp$wa_une)
# df_temp$q2 <- lm(hspss ~ ed4, data = df_temp) %$% predict(object = ., newdata = df_temp %>% mutate(ed4 = levels(ed4)[4]), weight = df_temp$wa_exp)
# df_temp$q3 <- lm(hspss ~ ed4 + pyc + obcat + depressed, data = df_temp) %$% 
#   predict(object = ., newdata = df_temp %>% mutate(ed4 = levels(ed4)[1]), weight = df_temp$wa_une)
# 
# 
# for(j in outcomes){
#   df_temp$q <- lm(hspss ~ sex + age3 + race2 + ed4 + pyc + obcat + depressed, data = df_temp) %$% predict(object = ., newdata = df_temp)
# }
# df_temp$q <- lm(hspss ~ sex + age3 + race2 + ed4 + pyc + obcat + depressed, data = df_temp) %$% predict(object = ., newdata = df_temp)
# 
# df_temp$q1 <- lm(hspss ~ ed4, data = df_temp) %$% predict(object = ., newdata = df_temp %>% mutate(ed4 = levels(ed4)[1]), weight = df_temp$wa_une)
# df_temp$q2 <- lm(hspss ~ ed4, data = df_temp) %$% predict(object = ., newdata = df_temp %>% mutate(ed4 = levels(ed4)[4]), weight = df_temp$wa_exp)
# df_temp$q3 <- lm(hspss ~ ed4 + pyc + obcat + depressed, data = df_temp) %$% 
#   predict(object = ., newdata = df_temp %>% mutate(ed4 = levels(ed4)[1]), weight = df_temp$wa_une)