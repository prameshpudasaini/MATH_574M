library(data.table)
library(e1071)
library(tree)

# data accessed from https://hastie.su.domains/ElemStatLearn/data.html

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SVM --------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

zip.train <- read.csv("ignore/zip.train.gz", header = FALSE, sep = " ")
zip.train <- as.data.table(zip.train)[V1 %in% c(2, 3), ]
zip.train[, V1 := factor(as.factor(V1), levels = c(2, 3))]

zip.test <- read.csv("ignore/zip.test.gz", header = FALSE, sep = " ")
zip.test <- as.data.table(zip.test)[V1 %in% c(2, 3), ]
zip.test[, V1 := factor(as.factor(V1), levels = c(2, 3))]

zip.train <- zip.train[, 1:257]

# Linear SVM

lsvm.fit <- svm(V1 ~ ., data = zip.train, kernel = 'linear')
lsvm.pred <- predict(lsvm.fit, newdata = zip.test[, -1])

lsvm.error <- round(mean(lsvm.pred != zip.test$V1), 6L)
lsvm.error

# Gaussian Kernel SVM

cost_gaussian_kernel <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)

for (i in seq_along(cost_gaussian_kernel)) {
    gsvm.fit <- svm(V1 ~ ., data = zip.train, kernel = 'radial', cost = cost_gaussian_kernel[i])
    gsvm.pred <- predict(gsvm.fit, newdata = zip.test[, -1])
    
    gsvm.test.error <- round(mean(gsvm.pred != zip.test$V1), 6L)
    
    print(paste0("Test error for cost '", cost_gaussian_kernel[i], "' is: ", gsvm.test.error))
}

# Polynomial Kernel SVM

degree_polynomial_kernel <- c(1, 2, 3)

for (i in seq_along(degree_polynomial_kernel)) {
    psvm.fit <- svm(V1 ~ ., data = zip.train, kernel = 'polynomial', degree = degree_polynomial_kernel[i])
    psvm.pred <- predict(psvm.fit, newdata = zip.test[, -1])
    
    psvm.test.error <- round(mean(psvm.pred != zip.test$V1), 6L)
    
    print(paste0("Test error for degree '", degree_polynomial_kernel[i], "' is: ", psvm.test.error))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Classification Tree ----------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

zip.train <- read.csv("ignore/zip.train.gz", header = FALSE, sep = " ")
zip.train <- as.data.table(zip.train)[V1 %in% c(1, 2, 3), ]
zip.train[, V1 := factor(as.factor(V1), levels = c(1, 2, 3))]

zip.test <- read.csv("ignore/zip.test.gz", header = FALSE, sep = " ")
zip.test <- as.data.table(zip.test)[V1 %in% c(1, 2, 3), ]
zip.test[, V1 := factor(as.factor(V1), levels = c(1, 2, 3))]

zip.train <- zip.train[, 1:257]

tree.fit <- tree(V1 ~ ., data = zip.train, split = 'deviance')

tree.train.pred <- predict(tree.fit, type = 'class')
tree.train.error <- round(mean(tree.train.pred != zip.train$V1), 6L)

tree.test.pred <- predict(tree.fit, newdata = zip.test[, -1], type = 'class')
tree.test.error <- round(mean(tree.test.pred != zip.test$V1), 6L)

tree.train.error
tree.test.error
