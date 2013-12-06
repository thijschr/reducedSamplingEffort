### Source file with functions
# Source it in R: source(functionsSamplingEffort.R)

# Function for selecting random rows/samples within a site ---------------

randomRows <- function(data, n) 
    ## data: dataframe with site (1st col), sample (2nd), and species abundances
    ## n: number of samples to sample within a site
{
    return(data[sample(nrow(data), n), ])
}

