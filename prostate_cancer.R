library(data.table)

# prostate cancer data set
url <- "https://hastie.su.domains/ElemStatLearn/datasets/prostate.data"
data <- fread(url, drop = 1)

train <- data[train == 'T', 1:9]
test <- data[train == 'F', 1:9]

X_train <- train[, 1:8]
y_train <- train$lpsa
X_test <- test[, 1:8]
y_test <- test$lpsa

set.seed(1000)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear Regression ------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lm_fit <- lm(lpsa ~ ., data = train)

summary(lm_fit)

lm_train_pred <- predict(lm_fit, X_train)
lm_test_pred <- predict(lm_fit, X_test)

lm_train_error <- mean((y_train - lm_train_pred) ^ 2)
lm_test_error <- mean((y_test - lm_test_pred) ^ 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Forward Selection ------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(leaps)

reg_fit <- regsubsets(lpsa ~ ., data = train, method = 'forward')

reg_fit_summary <- summary(reg_fit)
reg_fit_summary

names(reg_fit_summary)
reg_fit_summary$adjr2 # adjusted R-squared
reg_fit_summary$bic # BIC
which.min(reg_fit_summary$bic)

train_mat <- model.matrix(lpsa ~ ., data = train)
test_mat <- model.matrix(lpsa ~ ., data = test)

# regression coefficients, training error, BIC, and AIC for each model

reg_train_error <- vector('double', length(X_train))
reg_train_bic <- vector('double', length(X_train))
reg_train_aic <- vector('double', length(X_train))

for (i in seq_along(X_train)) {
    
    print(paste0("Model: ", i))
    
    coef_i <- coef(reg_fit, i)
    print(coef_i)
    
    reg_train_pred <- train_mat[, names(coef_i)] %*% coef_i
    reg_train_error[i] <- mean((y_train - reg_train_pred) ^ 2)
    
    n <- nrow(X_train)
    reg_train_bic[i] <- n * log(reg_train_error[i]) + log(n) * length(coef_i)
    reg_train_aic[i] <- n * log(reg_train_error[i]) + 2 * length(coef_i)
    
    print(paste0("Training error: ", round(reg_train_error[i], 4)))
    print(paste0("BIC: ", round(reg_train_bic[i], 4)))
    print(paste0("AIC: ", round(reg_train_aic[i], 4)))
    cat("\n")
}

# refit OLS using the best model
# create a function to select BIC or AIC as evaluation parameter for model selection

reg_best_mod <- function(eval_param) {
    
    num_coef <- which.min(eval_param)
    coef(reg_fit, num_coef)
    
    imp_var <- names(coef(reg_fit, num_coef))[-1] # important variables
    imp_var
    all_var <- append(imp_var, 'lpsa')
    all_var
    
    lm_reg_fit <- lm(lpsa ~ ., data = train[, ..all_var])
    
    summary(lm_reg_fit)
    
    lm_reg_test_pred <- predict(lm_reg_fit, test)
    lm_reg_test_error <- mean((y_test - lm_reg_test_pred) ^ 2)
    
    return(list(num_coef = num_coef, 
                imp_var = imp_var, 
                lm_reg_fit = lm_reg_fit,
                lm_reg_test_error = lm_reg_test_error))
}

## best model based on BIC
reg_best_mod_bic <- reg_best_mod(reg_train_bic)
reg_best_mod_bic$num_coef
reg_best_mod_bic$lm_reg_fit$coefficients
reg_best_mod_bic$lm_reg_test_error

## best model based on AIC
reg_best_mod_aic <- reg_best_mod(reg_train_aic)
reg_best_mod_aic$num_coef
reg_best_mod_aic$lm_reg_fit$coefficients
reg_best_mod_aic$lm_reg_test_error


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LASSO Regression -------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(lars)

lasso_fit <- lars(as.matrix(X_train), y_train)

plot(lasso_fit)
summary(lasso_fit)

# cross-validation
lasso_s <- seq(0, 1, length = 100)
lasso_cv <- cv.lars(as.matrix(X_train), y_train, K = 5, index = lasso_s)

names(lasso_cv)
lasso_mcv <- which.min(lasso_cv$cv)

# minimum CV rule
lasso_mcv_best_s <- lasso_s[lasso_mcv]

lasso_mcv_pred <- predict(lasso_fit, s = lasso_mcv_best_s, type = 'coef', mode = 'frac')
lasso_mcv_pred$s
lasso_mcv_pred$coefficients

lasso_mcv_test_pred <- predict(lasso_fit, newx = X_test, s = lasso_mcv_best_s, type = 'fit', mode = 'frac')
lasso_mcv_test_error <- mean((y_test - lasso_mcv_test_pred$fit) ^ 2)

# one-standard rule
bound <- lasso_cv$cv[lasso_mcv] + lasso_cv$cv.error[lasso_mcv]
lasso_st_best_s <- lasso_s[min(which(lasso_cv$cv < bound))]

lasso_st_pred <- predict(lasso_fit, s = lasso_st_best_s, type = 'coef', mode = 'frac')
lasso_st_pred$s
lasso_st_pred$coefficients

lasso_st_test_pred <- predict(lasso_fit, newx = X_test, s = lasso_st_best_s, type = 'fit', mode = 'frac')
lasso_st_test_error <- mean((y_test - lasso_st_test_pred$fit) ^ 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adaptive LASSO Regression ----------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compute weights
lm_coef <- lm_fit$coefficients

alasso_tune_param <- 1
alasso_w <- 1 / (abs(lm_coef) ^ alasso_tune_param)

# modify X
X_train_star <- X_train * (1 / alasso_w)

# fit adaptive LASSO
alasso_fit <- lars(as.matrix(X_train_star), y_train)

names(alasso_fit)
plot(alasso_fit)

# modify beta
alasso_fit$beta <- alasso_fit$beta * (1 / alasso_w)

# cross-validation
alasso_s <- seq(0, 1, length = 100)
alasso_cv <- cv.lars(as.matrix(X_train_star), y_train, K = 5, index = alasso_s)

names(alasso_cv)
alasso_mcv <- which.min(alasso_cv$cv)

# minimum CV rule
alasso_mcv_best_s <- alasso_s[alasso_mcv]

alasso_mcv_pred <- predict(alasso_fit, s = alasso_mcv_best_s, type = 'coef', mode = 'frac')
alasso_mcv_pred$s
alasso_mcv_pred$coefficients

X_test_star <- X_test * (1 / alasso_w)

alasso_mcv_test_pred <- predict(alasso_fit, newx = X_test_star, s = alasso_mcv_best_s, type = 'fit', mode = 'frac')
alasso_mcv_test_error <- mean((y_test - alasso_mcv_test_pred$fit) ^ 2)

# one-standard rule
bound <- alasso_cv$cv[alasso_mcv] + alasso_cv$cv.error[alasso_mcv]
alasso_st_best_s <- alasso_s[min(which(alasso_cv$cv < bound))]

alasso_st_pred <- predict(alasso_fit, s = alasso_st_best_s, type = 'coef', mode = 'frac')
alasso_st_pred$s
alasso_st_pred$coefficients

alasso_st_test_pred <- predict(alasso_fit, newx = X_test_star, s = alasso_st_best_s, type = 'fit', mode = 'frac')
alasso_st_test_error <- mean((y_test - alasso_st_test_pred$fit) ^ 2)
