## ----library, echo = T, results = 'hide', warning=FALSE, message=FALSE--------
library(MTPS)

## ----library2, message=FALSE, warning=FALSE, include=FALSE, results='hide'----
library(dplyr)
library(purrr)
library(tidyverse)
library(caret)
library(ggplot2)
library(reshape2)

## ----data---------------------------------------------------------------------
data("HIV")
head(YY)
XX[8:9, 7:10]
dim(YY)
dim(XX)

## ----ctn-split-data, echo = T, results = 'hide'-------------------------------
set.seed(12345)
xmat <- XX %>% as.matrix()
ymat <- YY %>% as.matrix()
nobs <- nrow(xmat)
id <- XX %>% 
  rowMeans() %>% 
  createFolds(k = 5, list = FALSE)
training_id <- seq_len(nobs) %>% 
  sample(size = 0.8 * nobs)

y_train <- ymat[training_id, ]
y_test  <- ymat[-training_id, ]
x_train <- xmat[training_id, ]
x_test  <- xmat[-training_id, ]

## ----ctn-noss, eval=FALSE, results = 'hide', eval=T---------------------------
 # no stacking
 fit_mult <- multiFit(xmat = x_train, ymat = y_train, method = glmnet1, family = "gaussian")


## ----ctn-train, eval=T, echo=T, message=FALSE, warning=FALSE, results='hide'----
# Standard Stacking 
fit_ss <- MTPS(xmat = x_train, ymat = y_train, family = "gaussian",
                            cv = FALSE, residual = FALSE,
                            method_step1 = glmnet1,
                            method_step2 = rpart1)
# Cross-Validation Stacking 
fit_cv <- MTPS(xmat = x_train, ymat = y_train, family = "gaussian",
                            cv = TRUE, residual = FALSE,
                            method_step1 = glmnet1,
                            method_step2 = rpart1)
# Residual Stacking 
fit_rs <- MTPS(xmat = x_train, ymat = y_train, family = "gaussian",
                            cv = FALSE, residual = TRUE,
                            method_step1 = glmnet1,
                            method_step2 = rpart1)
# Cross-Validation Residual Stacking 
fit_cvrs <- MTPS(xmat = x_train, ymat = y_train, family = "gaussian",
                            cv = TRUE, residual = TRUE,
                            method_step1 = glmnet1,
                            method_step2 = rpart1)

## ----echo=F,eval=T------------------------------------------------------------
data("Internal")

## ----ctn-predict, echo = T, results = 'hide', eval=T--------------------------
# no stacking
pred_mult <- predict(fit_mult, x_test)
# Standard Stacking 
pred_ss <- predict(fit_ss, x_test)
# Cross-Validation Stacking 
pred_cv <- predict(fit_cv, x_test)
# Residual Stacking 
pred_rs <- predict(fit_rs, x_test)
# Cross-Validation Residual Stacking 
pred_cvrs <- predict(fit_cvrs, x_test)

## ----ctn-outcome,  echo = T, eval=T, message=FALSE, warning=FALSE-------------
n_test <- nrow(x_test)

ctn_plot_data_matrix <- cbind(
  rbind(pred_mult, pred_ss, pred_cv, pred_rs, pred_cvrs),
  y_test[rep(seq_len(n_test), 5), ]
)
ctn_plot_data <- data.frame(
  method = rep(c("No-Stacking", "SS", "CV", "RS", "CVRS"), each = n_test),
  ctn_plot_data_matrix
)
colnames(ctn_plot_data) <- c("method", paste0("pred.", colnames(y_test)), colnames(y_test))

dm1 <- ctn_plot_data %>%
  select(method, ABC, `3TC`, AZT, D4T, DDI) %>%
  pivot_longer(cols = -method, names_to = "Y", values_to = "yVal")
dm2 <- ctn_plot_data %>%
  select(method, pred.ABC, pred.3TC, pred.AZT, pred.D4T, pred.DDI) %>%
  pivot_longer(cols = -method, names_to = "Y", values_to = "predictVal")

ctn_plot_data <- bind_cols(dm1, select(dm2, -method))
colnames(ctn_plot_data) <- c("method", "Y", "yVal", "predict", "predictVal")
ctn_plot_data <- ctn_plot_data %>%
  mutate(
    method = factor(method, levels = unique(method)),
    yVal = as.numeric(yVal),
    predictVal = as.numeric(predictVal)
  )

ggplot(ctn_plot_data) +
  geom_point(aes(x = predictVal, y = yVal, color = method), size = 0.5, alpha = 0.8) +
  geom_abline(slope = 1, alpha = 0.2) +
  coord_fixed() +
  ylab("Testing Data Outcome") + xlab("Predicted Outcome on Testing Data") +
  scale_x_discrete(breaks = NULL) +
  scale_y_discrete(breaks = NULL) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  facet_grid(Y ~ method)

## ----eval=FALSE, message=FALSE, warning=FALSE, include=FALSE, results='hide'----
#  n_test <- nrow(x_test)
#  
#  ctn_plot_data_matrix <- cbind(
#    rbind(pred_mult, pred_ss, pred_cv, pred_rs, pred_cvrs),
#    y_test[rep(seq_len(n_test), 5), ]
#  )
#  ctn_plot_data <- data.frame(
#    method = rep(c("No-Stacking", "SS", "CV", "RS", "CVRS"), each = n_test),
#    ctn_plot_data_matrix
#  )
#  colnames(ctn_plot_data) <- c("method", paste0("pred.", colnames(y_test)), colnames(y_test))
#  
#  dm1 <- ctn_plot_data %>%
#    select(method, ABC, `3TC`, AZT, D4T, DDI) %>%
#    pivot_longer(cols = -method, names_to = "Y", values_to = "yVal")
#  dm2 <- ctn_plot_data %>%
#    select(method, pred.ABC, pred.3TC, pred.AZT, pred.D4T, pred.DDI) %>%
#    pivot_longer(cols = -method, names_to = "Y", values_to = "predictVal")
#  
#  ctn_plot_data <- bind_cols(dm1, select(dm2, -method))
#  colnames(ctn_plot_data) <- c("method", "Y", "yVal", "predict", "predictVal")
#  ctn_plot_data <- ctn_plot_data %>%
#    mutate(
#      method = factor(method, levels = unique(method)),
#      yVal = as.numeric(yVal),
#      predictVal = as.numeric(predictVal)
#    )
#  
#  ggplot(ctn_plot_data, aes(x = predictVal, y = yVal, color = method)) +
#    geom_point(size = 0.5, alpha = 0.8) +
#    geom_abline(slope = 1, intercept = 0, alpha = 0.2) +
#    coord_fixed() +
#    labs(
#      y = "Testing Data Outcome",
#      x = "Predicted Outcome on Testing Data"
#    ) +
#    scale_x_continuous(breaks = NULL) +
#    scale_y_continuous(breaks = NULL) +
#    theme_bw() +
#    theme(
#      axis.text = element_blank(),
#      strip.placement = "outside",
#      strip.background = element_blank()
#    ) +
#    facet_grid(Y ~ method)
#  

## ----bin-data, echo = T, results = 'hide', eval=T-----------------------------
# https://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/CUTOFFS/drug.cutoffs
# cutoff value to be used to define drug resistent
cutoffs <- c(2,3,3,1.5,1.5) 
ymat_bin <- ymat
xmat_bin <- xmat
for(ii in 1:5) ymat_bin[,ii] <- (10^ymat[,ii] < cutoffs[ii]) * 1
y_train_bin <- ymat_bin[training_id, ]
y_test_bin  <- ymat_bin[-training_id, ]
x_train_bin <- xmat_bin[training_id, ]
x_test_bin  <- xmat_bin[-training_id, ]

## ----bin-train, echo = T, results = 'hide', eval= T---------------------------
fit_prs_std <- MTPS(xmat = x_train_bin, ymat = y_train_bin,
                               family = "binomial",
                               cv = FALSE, residual = TRUE,
                               method_step1 = rpart1,
                               method_step2 = lm1,
                               resid_type = "pearson", resid_std = TRUE) 
pred_prs_std <- predict(fit_prs_std, x_test_bin)

## ----bin-outcome, echo = F----------------------------------------------------
walk(1:ncol(y_test_bin), function(yy) {
  col_name <- colnames(y_test_bin)[yy]
  print(col_name)
  cm <- table((pred_prs_std[, yy] > 0.5) * 1, y_test_bin[, yy])
  print(cm)
})


## ----mix-data, echo = T, results = 'hide', eval=T-----------------------------
ymat.mix <- bind_cols(as_tibble(ymat)[, 1:3], as_tibble(ymat_bin)[, 4:5]) %>% as.matrix()
xmat.mix <- xmat
y_train_mix <- ymat.mix[training_id, ]
y_test_mix <- ymat.mix[-training_id, ]
x_train_mix <- xmat.mix[training_id, ]
x_test_mix <- xmat.mix[-training_id, ]

## ----mix-training, echo = T, results = 'hide', warning=FALSE, message=FALSE, eval=T----
fit_mix_rs <- MTPS(
  xmat = x_train_mix,
  ymat = y_train_mix,
  family = c("gaussian", "gaussian", "gaussian", "binomial", "binomial"),
  cv = FALSE,
  residual = TRUE,
  method_step1 = glmnet_lasso,
  method_step2 = rpart1
)
pred_mix_rs <- predict(fit_mix_rs, x_test_mix)

## ----mix-outcome, eval=T------------------------------------------------------
n_test <- nrow(x_test)

mix_plot_data <- tibble(
  Y = rep(colnames(y_test_mix)[1:3], each = nrow(y_test_mix)),
  predict = c(pred_mix_rs[, 1], pred_mix_rs[, 2], pred_mix_rs[, 3]),
  outcome = c(y_test_mix[, 1], y_test_mix[, 2], y_test_mix[, 3])
) %>%
  mutate(
    predict = as.numeric(predict),
    outcome = as.numeric(outcome)
  )

ggplot(mix_plot_data) +
  geom_point(aes(x = predict, y = outcome, color = Y),
             size = 0.5,
             alpha = 0.8) +
  ylab("Outcome of Testing Data") + xlab("Predicted Outcome of Testing Data") +
  scale_x_discrete(breaks = NULL) +
  scale_y_discrete(breaks = NULL) +
  geom_abline(slope = 1, alpha = 0.2) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.text = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  facet_grid( ~ Y)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  ggplot(mix_plot_data, aes(x = predict, y = outcome, color = Y)) +
#    geom_point(size = 0.5, alpha = 0.8) +
#    geom_abline(slope = 1,
#                intercept = 0,
#                alpha = 0.2) +
#    coord_fixed() +
#    labs(y = "Outcome of Testing Data", x = "Predicted Outcome of Testing Data") +
#    scale_x_continuous(breaks = NULL) +
#    scale_y_continuous(breaks = NULL) +
#    theme_bw() +
#    theme(
#      legend.title = element_blank(),
#      axis.text = element_blank(),
#      strip.placement = "outside",
#      strip.background = element_blank()
#    ) + facet_grid( ~ Y)

## ----echo=T, eval=T-----------------------------------------------------------
walk(4:5, function(yy) {
  col_name <- colnames(y_test_mix)[yy]
  print(col_name)
  cm <- table((pred_mix_rs[, yy] > 0.5) * 1, y_test_mix[, yy])
  print(cm)
})

## ----data2--------------------------------------------------------------------
data("simdat_mtps")
head(ymat)
head(xmat)

dim(ymat)
dim(xmat)

## ----fit, echo = T, results = 'hide', warning=FALSE, message=FALSE, eval=FALSE----
#  # fit Residual Stacking Model for Survival Data
#  fit <- MTPS(
#    xmat_train,
#    ymat_train,
#    xmat_list = xmat_list,
#    family = 'survival',
#    cv = FALSE,
#    residual
#    = TRUE,
#    nfold = 5,
#    method_step1 = surv,
#    method_step2 = lm1,
#    dist1 = "lognormal"
#  )
#  # predict the survival time on test set
#  pre3 <- predict.MTPS(fit, xmat_test)

## ----Simulat survival data, echo = T, results = 'hide', warning=FALSE, message=FALSE, eval=FALSE----
#  set.seed(1)
#  data("simdat_mtps")
#  # prepare training and test set
#  id.train <- sample(1:100, 80)
#  xmat_train <- xmat[id.train, ]
#  xmat_test <- xmat[-id.train, ]
#  ymat_train <- cbind(list(survival::Surv(ymat[id.train, "time01"], ymat[id.train, "status01"])), list(survival::Surv(ymat[id.train, "time02"], ymat[id.train, "status02"])))
#  
#  # Produce the Kaplan-Meier estimator
#  weights <- find_km_weights_mat(ymat[id.train, ], num_outcome = 2)
#  
#  # fit Residual Stacking Model for Survival Data
#  fit <- MTPS(
#    xmat_train,
#    ymat_train,
#    xmat_list = xmat_list,
#    family = 'survival',
#    cv = FALSE,
#    residual
#    = TRUE,
#    nfold = 5,
#    method_step1 = ela1,
#    method_step2 = lm1,
#    dist1 = "lognormal"
#  )
#  
#  # predict the survival time on test set
#  pre4 <- predict.MTPS(fit, xmat_test)

## ----mix-mtd, echo = T, warning=FALSE,message=FALSE, eval=FALSE---------------
#  fit_mixOut <- MTPS(
#    xmat = x_train,
#    ymat = y_train,
#    family = "gaussian",
#    method_step1 =
#      c(glmnet_lasso, glmnet_lasso, glmnet_lasso, lm1, lm1),
#    method_step2 =
#      c(rpart1, glmnet_lasso, glmnet_lasso, glmnet_lasso, glmnet_lasso)
#  )
#  
#  pred <- predict(fit_mixOut, x_test)

## ----method-mod, eval=F, echo=T, eval=T---------------------------------------
glmnet_lasso <- modify_parameter (glmnet1, alpha = 1)
glmnet_ridge <- modify_parameter (glmnet1, alpha = 0)

## ----method-new, eval=F, echo=T,eval=FALSE------------------------------------
#  glm1 <- function(xmat, ymat, family, ...) {
#    tmp0 <- data.frame(yy = ymat, xmat)
#    model <- glm(yy ~ ., data = tmp0, family = family, ...)
#    y.fitted <- fitted(model)
#    predFun <- function(model, xnew) {
#      predict(model, newdata = data.frame(xnew), type = "response")
#    }
#    return(list(
#      model = model,
#      y.fitted = y.fitted,
#      predFun = predFun
#    ))
#  }

## ----learner-compare, echo=T,eval=T-------------------------------------------
nsim <- 20
mse_lasso_lm <- matrix(NA, nrow = nsim, ncol = ncol(ymat)) %>%
  as_tibble(.name_repair = "minimal") %>%
  set_names(colnames(ymat))
mse_ridge_lm <- mse_lasso_lm

for (ii in 1:nsim) {
  set.seed(ii)
  # lasso stacking with lm
  mse_lasso_lm[ii,] <- cv_MTPS(xmat, ymat, family="gaussian",
                            cv = FALSE, residual = TRUE,
                            method_step1 = glmnet_lasso, method_step2 = lm1,
                            resid_std = FALSE)$continuous
  # ridge stacking with lm
  mse_ridge_lm[ii,] <- cv_MTPS(xmat, ymat, family="gaussian",
                            cv = FALSE, residual = TRUE,
                            method_step1 = glmnet_ridge, method_step2 = lm1,
                            resid_std = FALSE)$continuous
}

mse_data <- tibble(
  lasso = rowMeans(mse_lasso_lm, na.rm = TRUE),
  ridge = rowMeans(mse_ridge_lm, na.rm = TRUE)
)

mse_data <- mse_data %>%
  pivot_longer(cols = everything(), names_to = "Learner", values_to = "MSE")

ggplot(mse_data) +
  geom_boxplot(aes(x=Learner, y=MSE, fill=Learner)) +
  ggtitle("Boxplot of average Mean Square Error") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

## ----eval=FALSE, include=FALSE------------------------------------------------
#  
#  ggplot(mse_data, aes(x = Learner, y = MSE, fill = Learner)) +
#    geom_boxplot() +
#    ggtitle("Boxplot of Average Mean Square Error") +
#    theme_bw() +
#    theme(
#      legend.position = "none",
#      plot.title = element_text(hjust = 0.5),
#      axis.title.x = element_blank()
#    ) +
#    scale_fill_manual(values = c("lasso" = "#1f78b4", "ridge" = "#33a02c"))
#  

## ----fomula of MSE and C-index,echo = T, results = 'hide', warning=FALSE, message=FALSE, eval=FALSE----
#  
#  # MSE
#  mse <- function(pre, object){
#    # pre_time: the predicted event time
#    # obj: 2 columns, time & status
#    ind <- which(object$status == 1)
#    event_time <- object$time[ind]
#    pre_event_time <- pre[ind]
#    names(event_time) <- names(pre_event_time)
#    return (mean((log(event_time) - log(pre_event_time))^2))
#  }
#  
#  # C-index
#  c_index <- function(pre, object){
#    total.pairs <- 0
#    c <- 0
#    test <- cbind(pre, object)
#    test <- test[order(test[,2]),]
#    for (i in 1:nrow(test)){
#      if(i ==nrow(test)) break
#      else if(test[i,3]==0) next # only consider uncensored observations
#      else{
#        for (j in (i+1):nrow(test)){
#          total.pairs <- total.pairs+1
#          if(test[j,2] == test[i,2]) {if(test[j,1] > test[i,1] & test[j,1] == test[i,1]) c <- c + 0.5}
#          #"if(test[j,1] > test[i,1] & test[j,1] == test[i,1]) c <- c + 0.5"
#          #if we want add all possible results we can use straightly "c <- c + 0.5"
#  
#          if(test[j,2] > test[i,2]){if(test[j,1] > test[i,1]) {c <- c+1}}
#          if(test[j,2] > test[i,2]){if(test[j,1] == test[i,1]) {c <- c+0.5}}
#        }
#      }
#  
#    }
#    return (c/total.pairs)
#  }
#  
#  

