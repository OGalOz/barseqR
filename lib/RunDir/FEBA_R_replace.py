


"""
    In R,
        data frame row slices have the following syntax:
            df

    In R, L is a 'vector', where the member values are TRUE if the row has am, 
        and FALSE if otherwise.
     L = df$am == 0 

     i.e. '$' operator is used to access columns of data frame.
        df[2:3,]    # select 2nd and 3rd row


""" 

def FEBA_Fit(expsUsed, all_tbl, genes,
            genesUsed=NULL, strainsUsed=NULL, genesUsed12=NULL,
            minT0Strain=3, minT0Gene=30,
            minT0GeneSide=minT0Gene/2,
            minGenesPerScaffold=10,
            pred=CrudeOp(genes),
            okDay=TRUE, # OK to use Time0 from another day on the same lane, if necessary?
            okLane=TRUE, # OK to compare to Time0 from another lane, if necessary?
            metacol=1:7, # for all
            # names of experiments to ignore; experiments with Drop=TRUE are also ignored
            ignore=NULL,
            # ignore those below this threshold, unless ignore is set
            minSampleReads = getenv_numeric_or_default("FEBA_MIN_SAMPLE_READS", 200*1000),
            debug=FALSE, computeCofit=TRUE,
            dir="."):

    
