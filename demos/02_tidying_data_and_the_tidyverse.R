# You will now work through an example of some real data from The Genever Group. The arrangement and format of these data are typical of many protein and gene expression datasets so the processing is representative of that needed in a variety of situations.
# The data are mass spectrometry data of the soluble protein fraction from five immortalised mesenchymal stromal cell (MSC) lines.
# 
# The data are normalised protein abundances. Each row is a protein.
# 

# My project and working directory is "C:/Users/er13/Desktop/BIO00058M-Data-science-2020"
# I have a copy of Y101_Y102_Y201_Y202_Y101-5.csv in data-raw

# Open the file in excel

# Data description
# The cells lines are Y101, Y102, Y201, Y202 and Y101.5 and there are three replicates for each cell line arranged in columns. Also in the file are columns for:
#   
#   the protein accession
# the number of peptides used to identify the protein
# the number of unique peptides used to identify the protein
# a measure of confidence in that identification
# the maximum fold change between the mean abundances of two cell lines (i.e., highest mean / lowest mean)
# a p value for a comparison test between the highest mean and lowest mean
# a q value (corrected p value)
# a measure of the power of that test
# the cell line with the highest mean
# the cell line with the lowest mean
# the protein mass
# whether at least two peptides were detected for a protein.

# load tidyverse
library(tidyverse)

# import
# define file name
filesol <- "data-raw/Y101_Y102_Y201_Y202_Y101-5.csv"
# skip first two lines
sol <- read_csv(filesol, skip = 2) %>% 
  janitor::clean_names()


# 👀 the :: notation gives you access to a package's functions without first using the library() command.
# 
# This is useful when you want to use a single function from a package, or you need to specify which package when a function name is used in two loaded packages.
# 

# Filtering rows
# This dataset includes bovine serum proteins from the medium on which the cells were grown which need to be filtered out.
# We also filter out proteins for which fewer than 2 peptides were detected since we can not be confident about their identity. This is common practice for such proteomic data.

sol <- sol %>% 
  filter(str_detect(description,
                    "OS=Homo sapiens")) %>% 
  filter(x1pep == "x")

# 👀 str_detect(string, pattern) returns a logical vector according to whether 'pattern' is found in 'string'.
# 
# 👀 Notice that we have applied filter() twice using the pipe.

# 
# Processing cells contents
# It would be useful to extract the genename from the description and put it in a column.
# 
# One entry from the description column looks like this:
sol$description[1]
  
# The genename is after GN=. We need to extract the part of the string with the genename and put it in a new column.
# 

# Processing cells contents
# A way to problem-solve your way through this is work with one value carrying out one operation at a time until you've worked out what to do before implementing on an entire column.

# The pipe makes it especially easy to break problems down like this. 


# One step at a time on one value
# Extract the first value of the description to work with:
one_description <- sol$description[1]

# Extract the part of the string after GN= using a regex:
str_extract(one_description,"GN=[^\\s]+")

# [ ] means some characters
# ^ means 'not' when inside [ ]
# \s means white space
# the \ before is an escape character to indicate that the next character, \ should not be taken literally (because it's part of \s)
# + means one or more

# So GN=[^\\s]+ means GN= followed by one or more characters that are not whitespace. This means the pattern stops matching at the first white space after "GN=".

# We're close. Now we will drop the GN= part by replacing it with nothing:
# 
# Add replacing GN= with an empty string, "", to the pipeline:
str_extract(one_description, "GN=[^\\s]+") %>% 
  str_replace("GN=", "")

# Creating a new column
# Now we know how to get the result for one value, we need to apply the same process to the whole column
# 
# 
# mutate() is the dplyr function that adds new variables and preserves existing ones. It takes name = value pairs of expressions where:
# name is the name for the new variable and
# value is the value it takes. This is usually an expression.
# 
# Add a variable genename which contains the processed string from the description variable:
sol <- sol %>%
  mutate(genename =  str_extract(description,"GN=[^\\s]+") %>% 
           str_replace("GN=", ""))

