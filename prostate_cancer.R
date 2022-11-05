library(data.table)

# prostate cancer data set
url <- "https://hastie.su.domains/ElemStatLearn/datasets/prostate.data"
data <- fread(url)
data$V1 <- NULL

train <- data[train == 'T', 1:9]
test <- data[train == 'F', 1:9]

X_train <- train[, 1:8]
y_train <- train$lpsa
X_test <- test[, 1:8]
y_test <- test$lpsa


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
    
    imp_var <- names(coef(reg_fit, num_coef)) |> tail(-1) # important variables
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

## best mod based on BIC
reg_best_mod_bic <- best_mod(reg_train_bic)
reg_best_mod_bic$num_coef
reg_best_mod_bic$lm_reg_test_error

## best mod based on AIC
reg_best_mod_aic <- best_mod(reg_train_aic)
reg_best_mod_aic$num_coef
reg_best_mod_aic$lm_reg_test_error
