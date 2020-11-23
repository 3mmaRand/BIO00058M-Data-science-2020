library(tidyverse)
library(GGally)
library(caret)

# Import the data and the cell labels:
# the scRNASeq  
file <- "data-raw/scrna_data.csv"
rna <- read_csv(file)

# the cell labels. louvain is the name of the clustering method used to identify
# cell groups
file <- "data-raw/scrna_meta.csv"
meta <- read_csv(file) %>% select(louvain)

# what cell types are there
table(meta$louvain)

# Add the cell labels to the data:
rna$cell <- meta$louvain

# Split the dataset in to training and testing sets
ids <- createDataPartition(y = rna$cell,
                           p = 0.75,
                           list = FALSE)
str(ids)

2638 * 0.75

# Create the training set:
train <- rna %>% slice(ids)

# Create the testing set:
test <- rna %>% slice(-ids)

# Perform the LDA on the training data
lda <- train %>% 
  select(-cell) %>%
  MASS::lda(grouping = train$cell)

# model performance on training set
# predict from training
plda_train <- train %>% 
  select(-cell) %>%
  predict(object = lda)

# confusion matrix
confusionMatrix(plda_train$class,factor(train$cell))

# cell types in training
table(train$cell)

865 / 1981

# stats by class
# Sensitivity measures the proportion of positives that are correctly identified (true positives)
# Specificity measures the proportion of negatives that are correctly identified (true negatives)


# model performance on test set
# predict from test
plda_test <- test %>% 
  select(-cell) %>%
  predict(object = lda)

confusionMatrix(plda_test$class, factor(test$cell))

# will differ from yours because the ids are a random selection.

# Plot the training data and then the test data.
# Extract the scores from the training set with the cell names:

lda_labelled_train <- data.frame(plda_train$x,
                                 cell = train$cell)
# Extract the scores from the training set with the cell names:
lda_labelled_test <- data.frame(plda_test$x,
                                cell = test$cell)

# Create a scatter plot for the training data:
lda_labelled_train %>% 
  ggplot(aes(x = LD1, y = LD2, color = cell)) +
  geom_point()

# Based on this plot, you might be surprised by the accuracy of the model predictions on the training set - there seems to be a lot of overlap.
# 
# However, you are only looking at LD1 and LD2. There are many dimensions in this dataset and the separation of groups might not be obvious from the first to LD.
# 
# pairwise plotting of several LD
lda_labelled_train %>% 
  select(LD1:LD5, cell) %>% 
  ggpairs(aes(color = cell))
# You can see how LD1 really separates Megakaryocytes from the other cell types but that other LD are needed to distinguish all the cell types.

# Create a scatter plot for the test data:

lda_labelled_test %>% 
  ggplot(aes(x = LD1, y = LD2, color = cell)) +
  geom_point()

# There's a lot of overlap here. Perhaps we will better see the difference by examining additional LDs. However, remember that the predictions were less good on the test set so we would expect it to be difficult to distinguish all cells.

lda_labelled_test %>% 
  select(LD1:LD5, cell) %>% 
  ggpairs(aes(color = cell))

library(MASS)
