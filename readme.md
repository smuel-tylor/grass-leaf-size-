# Developmental and biophysical determinants of grass leaf size worldwide

The two folders shared here present code used by *Baird et al* 
- *Kew*
  -- Multiple regression analysis of a 1752 species dataset of leaf dimensions and average climate across species ranges, including validation based on plotting of residuals
  > !! The size of the dataset and number of models makes the Kew code very intensive to run


- *LA27*
  -- Demonstration of the approaches used for ANOVA and regression analyses that assessed leaf trait scaling and photosynthetic type comparisons, using the 27 species dataset collected at UCLA, and a matched phylogeny.
	
The code was built in R 3.6.3, and uses packages
- ape
- phytools
- PHYLOGR
- nlme
- nortest

These are called through ```library()```, so will need to be pre-installed for the code to run. Similarly, the code tries to make your life slightly easier by using ```setwd(choose.dir())```, hopefully that's enough...