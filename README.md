# Phil's SAE Code

This is the code shared by Phil for the Fae-Harriett estimator.   

The docs directory has the paper for this code, and a bunch of Markdown files (mostly generated with LLMs) to make upderstanding the code a little easier. 

The data directory (ignored by git due to size) holds the data.  The Rscripts folder has the R code, and sql folder holds the SQL query.  


## Code
This is a short (one line) description of the scripts
- The first scripts, 1, 2, and 2a extract data and get it in the correct format for further processing.
- The scripts 3 and 3a do some more post-processing (e.g. binning canopy eight) and output county level estimates of these statistics
- Script 4 glues the data created so far together.
- Then Script 5 fits the models and graphs some results based on those.  
- Script 6 refits the best model with more iterations, and then does some evaluations with that.
- The remaining files have some utility functions.  
    - `adj_YL_FisherScore_eblupFH.R` calculates the model fits (it's long)
    - `adj_YL_FisherScore_mseFH.R` calculates error estimates
    - `reproject_spatial_raster.R` has a utility function for reprojecting rasters (just imagine!)
- The SQL query aggregates over counties and calculates area-based estimates from the tree-level data.

## Blocks

I think that the best way to get started with writing blocks for this code is simply to take the entire scripts and dump them into blocks.  Then we can sub-divide them.  There will be some file management issues since the initial code uses absolute paths and our system handles them differently.  We can also convert the command line arguments to block inputs. Other than these things, I think converting into blocks wholesale should be (relatively) easy.  

It's also possible that there are so many inputs we need to update our (global) input handling.  The current system is best when blocks output one file of each type, and there might be multiple here in Phil's code.  This is a good thing to look at overall.  
