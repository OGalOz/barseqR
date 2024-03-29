#!/usr/bin/Rscript
# A command-line wrapper for running FEBA.R
# Can be invoked from the command-line to analyze the fitness experiments and save the results
# (i.e., FEBA_Fit() and FEBA_Save_Tables()) - the two main functions from FEBA.R
#   used to do computations.

usage = paste("",
        	"Usage: RunFEBA.R orgname data_directory FEBA_dir",
		"   Compute fitness values from genes, exps, pool, and all.poolcount",
		"   in the data directory.",
		"   If successful, creates a mini web site in the data_directory,",
		"   including tables of per-gene and per-strain fitness values,",
		"   plots and quality metrics to check if the experiments,",
		"   and an R image with everything from the complex fit data structure.",
		"   (RunFeba.R is normally invoked via BarSeqR.pl)",
		"", sep="\n");

RunFEBA = function(args = commandArgs(trailingOnly=TRUE)) {
        # Note: commandArgs is a built-in function that gets command line args
        # Note: data_dir must contain files all.poolcount, genes, exps, pool
        # Note: FEBAdir must contain /lib/FEBA.R, /lib/desc_short_rules
        # 
        # When running from Kbase App barseqR:
        # dir: 'workdir/outdir'
        # FEBAdir:
        #

	if (length(args) != 3) stop(usage);
	org = args[1];
	dir = args[2]
	FEBAdir = args[3];


	allfile = paste(dir,"/all.poolcount",sep="");
	genesfile = paste(dir,"/genes",sep="");
	expsfile = paste(dir,"/exps",sep="");
	poolfile = paste(dir,"/pool",sep="");
	FEBA_R = paste(FEBAdir,"/lib/FEBA.R",sep="");
        kb_test_exps_fp = paste(dir,"/KBtestexps",sep="");
        kb_test_all_fp = paste(dir,"/KBtest_all",sep="");

	# mode 4 means read permission
	if (file.access(allfile, mode=4) != 0) stop("Cannot read all file: ",allfile);
	if (file.access(genesfile, mode=4) != 0) stop("Cannot read genes file: ",genesfile);
	if (file.access(expsfile, mode=4) != 0) stop("Cannot read exps file: ",expsfile);
	if (file.access(poolfile, mode=4) != 0) stop("Cannot read pool file: ",poolfile);
	if (file.access(FEBAdir, mode=4) != 0) stop("Cannot access ",FEBAdir);
	if (file.access(FEBA_R, mode=4) != 0) stop("Cannot find ",FEBA_R);

        # This function, 'source', allows us to get functions from the FEBA.R file.
	source(FEBA_R);

        # Below functions convert all the tsv files into data frames
        # as.is being set to "TRUE" means convert character variables to factors
        # ("strings to integers?")
	rules = read.table(paste(FEBAdir,"/lib/desc_short_rules",sep=""),as.is=T);

	genes = read.delim(genesfile,quote="",as.is=T);
        for (n in c("scaffoldId","locusId","sysName","desc"))
          if(!n %in% names(genes)) stop("genes table must include field ", n);
        # below check.names checks if names are syntactically valid variable names, 
        #   so setting it to False (F) means we don't check.
	all = read.delim(allfile,as.is=T,check.names=F);
	exps = read.delim(expsfile,as.is=T);
        # Below we get all unique locus Ids from the all.poolcount file
	d = unique(all$locusId[all$locusId != "" & !is.na(all$locusId)]);
        # Then we check if for each unique locusId in d, it also exists in genes?
	if (!all(d %in% genes$locusId)) stop("Unknown genes in ",allfile);
	cat(sprintf("Read %d genes, %d exps, and data for %d barcodes\n",
	            nrow(genes), nrow(exps), nrow(all) ));

	# fix up names of all to be shorter
	SetNames = unique(exps$SetName);
	SetNames2 = ShortSetNames(SetNames);
	expNamesNew = paste(SetNames2[match(exps$SetName, SetNames)], exps$Index, sep="");

        # What is happening below?
	exps$num = 1:nrow(exps);
	exps$name = paste(exps$SetName,exps$Index,sep=".");
        # applyRules function in file lib/FEBA.R
	exps$short = applyRules(rules, exps$Description);

        # Below creates a vector 1, 2, 3 ..., 7
	metacol = 1:7;
        
        # names(data_frame) returns a string vector with the names of the dataframe
        # name_list is the column names of 'all', except the metadata
	name_list = names(all)[-metacol];
	if(length(name_list) != nrow(exps)) stop("Number of data columns in  ",allfile," does not match number of rows in ",expsfile);
	if(any(name_list != exps$name)) stop("Column names in  ",allfile," do not match names from ",expsfile);

	# remove trailing spaces from Group, Condition_1, Condition_2 in exps
	for(n in c("Group","Condition_1","Condition_2")) {
	    if(!is.null(exps[[n]])) exps[[n]] = sub(" +$", "", exps[[n]]);
	}

	names(all)[-metacol] = expNamesNew;
	exps$name = expNamesNew;

        # load strains to use, if this data exists
        strainsUsed = NULL;
        genesUsed = NULL;
        genesUsed12 = NULL;
        fStrainsUsed = paste(dir,"/strainusage.barcodes",sep="");
        fGenesUsed = paste(dir,"/strainusage.genes",sep="");
        fGenesUsed12 = paste(dir,"/strainusage.genes12",sep="");
        # Uses the strain usage files if they exist
	if (file.access(fStrainsUsed, mode=4) == 0) {
        	stopifnot(file.access(fGenesUsed,mode=4)==0);
        	stopifnot(file.access(fGenesUsed12,mode=4)==0);
		barcodesUsed = scan(fStrainsUsed,"");
                strainsUsed = all$barcode %in% barcodesUsed;
		genesUsed = scan(fGenesUsed,"");
		genesUsed12 = scan(fGenesUsed12,"");
                cat(sprintf("Loaded %d strains and %d genes to include in the analysis\n",
                	    length(barcodesUsed), length(genesUsed)));
	}

	options(width=100);


        # DEBUG: Writing Table out to file.
        write.table(exps, kb_test_exps_fp, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

        write.table(all, kb_test_all_fp, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

        
        # 'exps' and 'all' and 'genes' come from read.delim on files. All are dataframes
        # strainsUsed, genesUsed, genesUsed12 are all optional
	fit = FEBA_Fit(exps, all, genes, dir=dir,
		strainsUsed=strainsUsed, genesUsed=genesUsed, 
                genesUsed12=genesUsed12,
                minSampleReads=0);
	# FEBA_Save_Tables(fit, genes, org, expsU=exps, dir=dir, FEBAdir=FEBAdir);
	write(date(), paste(dir,"/.FEBA.success", sep="")); # report success to BarSeqR.pl
}

ShortSetNames = function(sets) {
    # sets is a list of set names (str)
    # We shorten each to fit into format defined below
	simple = grepl("(set|test)[0-9A-Z]+[0-9A-Z0-9]*$", sets);
	sets[simple] = sub("^.*(set|test)", "\\1", sets[simple]);
	nleft = sum(!simple);
	candidates = strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ","")[[1]];
	candidates = paste("set",candidates,sep="");
	candidates = setdiff(candidates, sets[simple]);
	if (nleft > length(candidates)) stop("Too many wierd set names");
	oldComplex = sets[!simple];
	sets[!simple] = candidates[1:nleft];
	if(nleft > 0) for (i in 1:nleft) {
	   write(sprintf("Set %s simplified to %s", oldComplex[i], sets[!simple][i]), stderr());
	}
	if (length(sets) != length(unique(sets))) stop("Non-unique sets");
	return(sets);
}

# Actually do the work
if(!interactive()) {
	RunFEBA();
	quit();
}
