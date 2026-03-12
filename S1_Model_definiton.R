
# Set the base directory using your favorite method
# setwd("...")


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# LOAD PACKAGES (BEGINNING)
##################################################################################################
library(Hmsc)
##################################################################################################
# LOAD PACKAGES (END)
##################################################################################################


##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################


##################################################################################################
# READ AND EXPLORE THE DATA (BEGINNING)
##################################################################################################
# Load the data
load(file=file.path(dataDir,"allDataGOOD.RData"))


dim(Zoo_data)
head(Zoo_data)
View(Zoo_data)
Y.con = Zoo_data


env_trans <- env_transformed
dim(env_trans)
head(env_trans)
View(env_trans)


trait <- subset(trait, select = -Capture)
dim(trait)
head(trait)
rownames(trait) <- colnames(Y.con)
View(trait)


phylo_tree
plot(phylo_tree,cex=0.5)


# We are ready to set up the HMSC model
##################################################################################################
# READ AND EXPLORE THE DATA (END)
##################################################################################################

##################################################################################################
# SET UP THE MODEL (BEGINNING)
##################################################################################################


XFormula = ~ conductivity+ pH + water_temp + water_depth + colour + chla+ lake_area + Tavg + precip + human + forest


TrFormula = ~ Body_size+Non_predators+Predators+Verification+Primary_filtration+Secondary_filtration+Suction+Crawling+Attachment


# We next define the Hmsc model.

m = Hmsc(Y=Y.con,
         TrData = trait, TrFormula = TrFormula,
         phyloTree = phylo_tree,
         XData = env_trans,  XFormula=XFormula,
         distr="lognormal poisson", studyDesign = study_design, ranLevels= list(Index_random_geo=rL))

# We included environmental covariates through XData and XFormula

# We included trait covariates by TrData and TrFormula

# We included the phylogenetic relationships by phyloTree

# We assumed the lognormal poisson distribution as appropriate for
# abundance data

# The data have specific structure and we include random
# effects, studyDesign & ranLevels

# It is always a good idea to look at the model object.

m


##################################################################################################
# SET UP THE MODEL (END)
##################################################################################################


##################################################################################################
# In the general case, we could have multiple models. We combine them into a list and given them names.
# COMBINING AND SAVING MODELS (START)
models = list(m)
names(models) = c("zoo.pa.model")
save(models, file = file.path(modelDir, "unfitted_models.RData"))
##################################################################################################
# COMBINING AND SAVING MODELS (END)
##################################################################################################


##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (START)
##################################################################################################
for(i in 1:length(models)){
  print(i)
  sampleMcmc(models[[i]],samples=2)
}
##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (END)
##################################################################################################
