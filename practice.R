library(tidyverse)
library(Rtsne)
library(GGally)
library(caret)
file <- "data-raw/scrna_data.csv"
rna <- read_csv(file)

file <- "data-raw/scrna_meta.csv"
meta <- read_csv(file) %>% select(louvain)


rna$louvain <- meta$louvain




pca <- rna %>% 
  prcomp(scale. = TRUE)

dat <-  data.frame(pca$x)
ggplot(dat, aes(x = PC1, y = PC2)) +
  geom_point()

dat <-  data.frame(pca$x, type = meta$louvain)
ggplot(dat, aes(x = PC1, y = PC2, colour = type)) +
  geom_point()


tsne <- rna %>% 
  Rtsne(perplexity = 30,
        check_duplicates = FALSE)

dat <- data.frame(tsne$Y)

dat %>% ggplot(aes(x = X1, y = X2)) +
  geom_point()

dat <- data.frame(tsne$Y, type = meta$louvain)

dat %>% ggplot(aes(x = X1, y = X2, colour = type)) +
  geom_point()

#####
lda <- rna %>% 
  select(-louvain) %>% 
  MASS::lda(grouping = rna$louvain)

plda <- rna %>% 
  select(-louvain)  %>%
  predict(object = lda)
table(predicted = plda$class, observed = rna$louvain)

lda_labelled <- data.frame(plda$x,
                           louvain = rna$louvain)

ggplot(lda_labelled, aes(x = LD1,
                         y = LD4,
                         colour = louvain)) +
  geom_point()


lda_labelled %>% 
  select(LD1:LD5, louvain) %>% 
  ggpairs(aes(color = louvain))



file <- "../LGG.csv"
lgg <- read_csv(file)


pca <- lgg %>% 
  select(-patientID, -Mutacion) %>% 
  prcomp(scale. = TRUE)
dat <-  data.frame(pca$x, mut = factor(lgg$Mutacion))
ggplot(dat, aes(x = PC1, y = PC2, colour = mut)) +
  geom_point()



tsne <- Rtsne(lgg, 
              perplexity = 20,
              check_duplicates = FALSE,
              pca = FALSE)

dat <- data.frame(tsne$Y,  mut = factor(lgg$Mutacion))

dat %>% ggplot(aes(x = X1, y = X2, colour = mut)) +
  geom_point()

dat %>% 
  select(PC1:PC5, mut) %>% 
  ggpairs(aes(color = mut))



file <- "data-raw/Liver_GSE14520_U133A.csv"
liver <- read_csv(file)

pca <- liver %>%
  select(-samples, -type) %>%
  prcomp(scale. = TRUE,
         rank. = 100)
dat <- data.frame(pca$x, type = liver$type)
ggplot(dat, aes(x = PC1, y = PC2, colour = type)) +
  geom_point()

tsne <- dat %>% 
  select(-type) %>% 
  Rtsne(perplexity = 30,
        check_duplicates = FALSE)

dat2 <- data.frame(tsne$Y,  type = dat$type)
dat2 %>% ggplot(aes(x = X1, y = X2, colour = type)) +
  geom_point()

# doesn't woek well for liver
ldaout <- liver %>% 
  select(-samples, -type) %>% 
  MASS::lda(grouping = liver$type)

pldaout <- liver %>% 
  select(-samples, -type)  %>%
  predict(object = ldaout)
table(predicted = pldaout$class, observed = liver$type)

ldaout_labelled <- data.frame(pldaout$x,
                              type = liver$type)

ggplot(ldaout_labelled, aes(x = LD1, y = LD2, colour = type)) +
  geom_point()
  

file <- "data-raw/Ovary_GSE12470.csv"
ovary <- read_csv(file)

pca <- ovary %>%
  select(-samples, -type) %>%
  prcomp(scale. = TRUE)

dat <- data.frame(pca$x, type = ovary$type)
ggplot(dat, aes(x = PC1, y = PC2, colour = type)) +
  geom_point()

dat %>% 
  select(PC1:PC8, type) %>% 
  ggpairs(aes(color = type))


tsne <- ovary %>% 
  select(-samples, -type) %>% 
  Rtsne(perplexity = 10,
        check_duplicates = FALSE)
dat2 <- data.frame(tsne$Y,  type = ovary$type)
dat2 %>% ggplot(aes(x = X1, y = X2, colour = type)) +
  geom_point()


ldaout <- ovary %>% 
  select(-samples, -type) %>% 
  MASS::lda(grouping = ovary$type)

pldaout <- ovary %>% 
  select(-samples, -type)  %>%
  predict(object = ldaout)
table(predicted = pldaout$class, observed = ovary$type)

ldaout_labelled <- data.frame(pldaout$x,
                              type = ovary$type)

ggplot(ldaout_labelled, aes(x = LD1, y = LD2, colour = type)) +
  geom_point()

#### BREAST
file <- "data-raw/Breast_GSE70947.csv"
breast <- read_csv(file)

pca <- breast %>%
  select(-samples, -type) %>%
  prcomp(scale. = TRUE)
dat <- data.frame(pca$x, type = breast$type)
ggplot(dat, aes(x = PC1, y = PC2, colour = type)) +
  geom_point()

dat %>% 
  select(PC1:PC10, type) %>% 
  ggpairs(aes(color = type))


# tsne <- dat %>% 
#   select(-type) %>% 
#   Rtsne(perplexity = 30,
#         check_duplicates = FALSE)

# dat2 <- data.frame(tsne$Y,  type = dat$type)
# dat2 %>% ggplot(aes(x = X1, y = X2, colour = type)) +
#   geom_point()

# w/o training and testing
lda <- breast  %>% 
  select(-samples, -type) %>% 
  MASS::lda(grouping = breast$type)

plda <- breast %>% 
  select(-samples, -type) %>%
  predict(object = lda)
# Examining the confusion matrix:
confusionMatrix(plda$class, factor(breast$type))



ids <- createDataPartition(y = breast$type,
                           p = 0.75,
                           list = FALSE)
train <- breast %>% slice(ids)
test <- breast %>% slice(-ids)

lda <- train %>% 
  select(-samples, -type) %>% 
  MASS::lda(grouping = train$type)

# predict on the training set
plda_train <- train %>% 
  select(-samples, -type) %>%
  predict(object = lda)
# Examining the confusion matrix:
confusionMatrix(plda_train$class, factor(train$type))

# predict on the test set
plda_test <- test %>% 
  select(-samples, -type) %>%
  predict(object = lda)

# Examining the confusion matrix:
confusionMatrix(plda_test$class, factor(test$type))

#  Extract the scores from the training set with the cell names:
lda_labelled_train <- data.frame(plda_train$x,
                                 cell = train$cell)
# Extract the scores from the training set with the cell names:
lda_labelled_test <- data.frame(plda_test$x,
                                cell = test$cell)


# Create a scatter plot for the training data:
lda_labelled_train %>% 
  ggplot(aes(x = LD1, y = LD2, color = type)) +
  geom_point()

# Select the first 5 LDs and pipe in to ggpairs():
lda_labelled_train %>% 
  select(LD1:LD5, cell) %>% 
  ggpairs(aes(color = cell))


# Now consider the test set.
# Create a scatter plot for the test data:
lda_labelled_test %>% 
  ggplot(aes(x = LD1, y = LD2, color = cell)) +
  geom_point()

# Select the first 5 LDs and pipe in to ggpairs():
lda_labelled_test %>% 
  select(LD1:LD5, cell) %>% 
  ggpairs(aes(color = cell))
##############################
## PCA LDA demo
n <- 50
m <- 0.8
min <- 2
max <- 6
x <- runif(n, min, max)
y <-  x*m + rnorm(n)
gp <- rep(c("a","b"), each = n/2)
df <- data.frame(x, y, gp)
df %>% ggplot(aes(x = x, y = y, colour = gp)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 8)) +
  scale_y_continuous(limits = c(0,8))

df$gp[df$y < 3] <- "a"
df$gp[df$y > 4.5] <- "b"
df %>% ggplot(aes(x = x, y = y, colour = gp)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 8)) +
  scale_y_continuous(limits = c(0,8))

write.table(df, "data-raw/df.txt", row.names = FALSE, quote = FALSE)

pca <- df %>% 
  select(x, y) %>%
  prcomp(scale. = TRUE)
pca$rotation
dat <- data.frame(pca$x, gp = df$gp)
ggplot(dat, aes(x = PC1, y = PC2, colour = gp)) +
  geom_point()


lda <- df %>% 
  select(x, y) %>%
  MASS::lda(grouping = df$gp)
lda$scaling

plda <- df %>% 
  select(x, y) %>%
  predict(object = lda)
dat <- data.frame(plda$x, type = df$gp)
ggplot(dat, aes(x = LD1, y = 0,
                colour = gp)) +
  geom_point()





