## How many species are required to maintain ecosystem functioning at large spatial scales?

This repository contains data and code for modelling the response of dry biomass productivity to species loss in a macroalgae system. It is associated with the following publication:

+ *coming soon (hopefully)*

## data

All raw data along with relevant metadata can be found in the *data* folder.

## code

The code to reproduce the analysis can be found in the *code* folder. The code will directly call the relevant data from the *data* folder. All figures and tables reported in the analysis will be outputted into a folder called *figures-tables* that will be created in the project directory. The scripts are numbered 01 to 08 and should be run in that order.

## output

The folder called *output* contains various outputs that are generated by the scripts during the analysis. These are then called by other scripts.

## renv

This project uses the *renv* R-package for package management. The .lock file contains all relevant information on the packages and the versions of those packages that were used in this project. To reproduce this analysis, users should install the renv R package:

+ install.packages("renv")

Then run the following code in the console:

+ renv::restore()

This will create a local copy of all relevant package versions that were used to perform these analyses.

## DESCRIPTION and LICENSE files

The description and licence files contain a variety of information regarding the licence under which the code was published, the contact details of the creator (me) and some notes on the package versions used to perform the analysis.

