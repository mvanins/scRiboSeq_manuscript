#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-05-22

library(optparse)
library(dplyr)
library(tidyr)
library(data.table)
library(feather)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(cowplot)
library(forcats)
library(mlr)
library(ranger)
library(parallelMap)
theme_set(theme_cowplot())

source('/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/scRibo/plotFunctions.R')

load('modelReads.RData')
figurePrefix <- "train100_0529"
set.seed(1234)

train_set <- sample(nrow(modelReads), size = nrow(modelReads)*.8)
test_set <- (1:nrow(modelReads))[-train_set]
trainReads <- modelReads[train_set,]



task <- makeClassifTask(data = trainReads, target = "cd5_stop")

ranger.learner <- makeLearner("classif.ranger",
							  predict.type = "response")

parallelStartSocket(parallel::detectCores())
res <- tuneParams(ranger.learner, 
                  task = task,
                  resampling = makeResampleDesc("CV", iters = 5L), 
                  par.set = makeParamSet(makeIntegerParam("num.trees", lower = 500, upper = 500),
                                         makeIntegerParam("mtry", lower = 13, upper = 27),
                                         makeIntegerParam("min.node.size", lower = 1, upper = 5)),
                  control = makeTuneControlRandom(maxit=100), 
                  measure = acc)
parallelStop()


assessment.task <- makeClassifTask(data = modelReads, target = "cd5_stop")


tuned.ranger.learner <- makeLearner("classif.ranger",
                                    predict.type = "response",
                                    importance = 'permutation',
                                    respect.unordered.factors = 'partition',
                                    mtry = res$x$mtry,
                                    num.trees = res$x$num.trees)


tuned.model <- train(tuned.ranger.learner, assessment.task, subset = train_set)
tuned.predict <- predict(tuned.model, assessment.task, subset = test_set)

performance(tuned.predict, measures=list(acc))

plt_tuned_predict_tt <- tuned.predict$data %>%
  select(truth,response) %>%
  group_by(truth, response) %>%
  summarize(n=n()) %>%
  ggplot(aes(x=truth,y=response,fill=n)) + 
  geom_tile() +
  scale_fill_gradientn(colours=rev(brewer.pal(7,"RdYlBu")),na.value=rev(brewer.pal(7,"RdYlBu"))[1]) +
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme( panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle(paste0("tuned.predict truth table. acc=",performance(tuned.predict, measures=list(acc))[1]))

save_plot(paste0(figurePrefix, "tuned.predict_truthtable.pdf"), plt_tuned_predict_tt, base_height = 5, base_width =5)


modelInfo <- data.frame(importance = tuned.model$learner.model$variable.importance,
						names = as.factor(names(tuned.model$learner.model$variable.importance)))
plt_importance <- modelInfo %>%
  mutate(names = fct_reorder(names,importance)) %>%
  ggplot(aes(x=names,y=importance))+
  geom_bar(stat="identity")+coord_flip()

save_plot(paste0(figurePrefix,"_importance.pdf"),plt_importance,base_height=5,base_width=5)

save(tuned.model,file="fucci.tuned.model_30-45.RData")

outer <- makeResampleDesc("CV", iters = 5)
forestWrapper <- makeTuneWrapper("classif.ranger",
                                 resampling = makeResampleDesc("CV", iters = 5L),
                                 par.set = makeParamSet(makeIntegerParam("num.trees", lower = 500, upper = 500),
                                                        makeIntegerParam("mtry", lower = 13, upper = 27),
                                                        makeIntegerParam("min.node.size", lower = 1, upper = 5)),
                                 control = makeTuneControlRandom(maxit=100))
parallelStartSocket(parallel::detectCores())
cvWithTuning <- resample(forestWrapper, task, resampling = outer)
parallelStop()

save.image("fucci.model.tuning_30-45.RData")



## CV on prediction accuracy
rdesc <- makeResampleDesc("CV", iters = 5)

parallelStartSocket(parallel::detectCores())
r <- resample(tuned.ranger.learner, 
              assessment.task, 
              rdesc,
              measures = list(acc, mmce))
parallelStop()
save.image("fucci.model.tuning_30-45_PredictionCVNumbers.RData")



plt_tuned_predict_tt_pct <- tuned.predict$data %>%
  select(truth,response) %>%
  group_by(truth, response) %>%
  summarize(n=n()) %>%
  ggplot(aes(x=truth,y=response,fill=n)) + 
  geom_tile() +
  scale_fill_gradientn(colours=(brewer.pal(9,"Greys")),na.value=(brewer.pal(9,"Greys"))[1]) +
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme( panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle(paste0("tuned.predict truth table. acc=",performance(tuned.predict, measures=list(acc))[1]))

save_plot(paste0(figurePrefix, "tuned.predict_truthtable_pct.pdf"), plt_tuned_predict_tt_pct, base_height = 5, base_width =5)
