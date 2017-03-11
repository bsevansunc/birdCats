library('tidyverse')
library('marked')
library('stringr')

# Read in data

bird <- read_csv('data/ch.csv')


# Temporary bird2 frame for tests
bird2 <- bird %>%
  select(imp,tCats,sex,ch)

# Remove active/inactive designation
bird2$ch <- str_replace_all(bird2$ch, '[.]', 0)

# Set sex as a factor with levels 'M', 'F', and 'U'
bird2$sex <- as.factor(bird2$sex)

# Model only runs with 65 unique capture histories or fewer 
# (with the data subset to these variables, 73 rows collapses to 65 histories)
bird2 <- bird2[1:73,]

# Process the data
bird.proc <- process.data(bird2)

# Design parameters
design.Phi <- list(static = c('imp', 'tCats', 'sex'))
design.p <- list(static = c('sex'))
design.parameters <- list(Phi = design.Phi, p = design.p)

bird.ddl <- make.design.data(bird.proc, parameters = design.parameters)

# Model parameters
Phi.test <- list(formula = ~imp+tCats+sex+imp*tCats+sex*tCats+imp^2)
p.dot <- list(formula = ~1)

model.test <- crm(bird.proc, ddl, model.parameters = list(Phi = Phi.test, p = p.dot))
