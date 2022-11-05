library(data.table)

# prostate cancer data set
url <- "https://hastie.su.domains/ElemStatLearn/datasets/prostate.data"
data <- fread(url)

data$V1 <- NULL

train <- data[train == 'T', 1:9]
test <- data[train == 'F', 1:9]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Linear Regression ------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lm_fit <- lm(lpsa ~ ., data = train)

summary(lm_fit)

lm_train_pred <- predict(lm_fit, train)
lm_test_pred <- predict(lm_fit, test)

lm_train_error <- mean((train$lpsa - lm_train_pred) ^ 2)
lm_test_error <- mean((test$lpsa - lm_test_pred) ^ 2)
