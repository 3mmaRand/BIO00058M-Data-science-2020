library(tidyverse)

file <- "../scrna_data.csv"
rna <- read_csv(file)


file <- "../scrna_meta.csv"
meta <- read_csv(file)




pca <- rna %>% 
  prcomp(scale. = TRUE)

dat <-  data.frame(pca$x)
ggplot(dat, aes(x = PC1, y = PC2)) +
  geom_point()


dat <-  data.frame(pca$x, type = meta$louvain)
ggplot(dat, aes(x = PC1, y = PC2, colour = type)) +
  geom_point()

library(Rtsne)
tsne <- rna %>% 
  Rtsne(perplexity = 30,
        check_duplicates = FALSE)

dat <- data.frame(tsne$Y, type = meta$louvain)

dat %>% ggplot(aes(x = X1, y = X2, colour = type)) +
  geom_point()





file <- "../LGG.csv"
lgg <- read_csv(file)


pca <- lgg %>% 
  select(-patientID, -Mutacion) %>% 
  prcomp(scale. = TRUE)
dat <-  data.frame(pca$x, mut = factor(lgg$Mutacion))
ggplot(dat, aes(x = PC1, y = PC2, colour = mut)) +
  geom_point()



tsne <- Rtsne(lgg, 
              perplexity = 30,
              check_duplicates = FALSE)

dat <- data.frame(tsne$Y,  mut = factor(lgg$Mutacion))

dat %>% ggplot(aes(x = X1, y = X2, colour = mut)) +
  geom_point()


file <- "../Breast_GSE45827.csv"
breast <- read_csv(file)

pca <- breast %>%
  select(-samples, -type) %>%
  prcomp(scale. = TRUE,
         rank. = 50)
dat <- data.frame(pca$x)
ggplot(dat, aes(x = PC1, y = PC2)) +
  geom_point()
screeplot(pca)

dat <- data.frame(pca$x, type = breast$type)
ggplot(dat, aes(x = PC1, y = PC2, colour = type)) +
  geom_point()

tsne <- dat %>% 
  select(-type) %>% 
  Rtsne(perplexity = 30,
        check_duplicates = FALSE)



dat2 <- data.frame(tsne$Y,  type = dat$type)
dat2 %>% ggplot(aes(x = X1, y = X2, colour = type)) +
  geom_point()

