library(foreign)
library(rpart)
library(randomForest)
library(dplyr)
library(tableone)
library(caret)


df <- read.dta("oaiShortAllT_2015.dta") %>% filter(!is.na(race2) & !is.na(edcv) & time == "BL")
tree <- rpart(edcv ~ age2 + sex + race2, 
              data = df[df$time == "BL",], 
              method = "class",
              control = rpart.control(cp = 0.001, minsplit = 2)
              )
printcp(tree)
df$w_a <- NA
for(i in 1:6) df[[paste0("p",i)]] <- NA
df[,c("p1", "p2", "p3", "p4", "p5","p6")] <- predict(tree, type = "prob")
freq <- table(df$edcv)
for(i in 1:6) df[as.numeric(df$edcv) == i,]$w_a <- 
              (freq[i]/nrow(df)) / df[[paste0("p",i)]][as.numeric(df$edcv) == i]

ranfor <- randomForest(ed3 ~ age2 + sex + race2, 
                       data = df, 
                       mtry = 2, 
                       ntree = 5000
                       )
predict(ranfor, type = "prob")
CreateTableOne(data=df, vars= c("edcv", "ed4", "ed3"), strata = c("sex", "race2", "age3"))

ndf <- df[,c("edcv","age3", "sex", "race2")]
model <- train(edcv ~., ndf, method = 'rf')
