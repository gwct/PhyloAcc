# Parameter file
In the parameter file, each parameter is specified in a line with the parameter name at the beginning followed by paremeter value. The parameters are: 
* **Input and output**: 
  * *PHYTREE_FILE*: the path of Phylogeny (.mod)  
  * *SEG_FILE*: the path of bed file for genomic regions
  * *ALIGN_FILE*: the path of multiple alignment file (.fasta)
  * *RESULT_PREFIX*: the prefix for output files
  * *id_path* (optional): only compute elements in this file. (The element is tagged by its order in the input bed file starting from 0). If not specified, the program will compute all elements in the input file.  
 
* **Specify species on the phylogeny**:
  * *TARGETSPECIES*: species of interest. E.g. species potentially lost conservation or with convergent phenotype changes.
  * *OUTGROUP*: outgroup species of the phylogeny. These species are not considered to be accelerated in our model. 
  * *CONSERVE*: species assumed to be mostly conserved. The algorithm will filter out elements "missing" in more than *CONSERVE_PROP* of the conserved species. 
  * *CONSERVE_PROP*: Default 0.8.
  * *GAPCHAR*: The character for alignment gaps (default: *). 

* **(Hyper)Parameters for the model**:
  * *PREP_GRATE*: the prior transition probabily from neutral to conserved state (default: 0.8).
  * *PREP_LRATE*: the prior transition probabily from conserved to accelerated state (default: 0.1).
  * *CONSERVE_PRIOR_A*: the shape parameter for the gamma prior of conserved rate (default: 10).
  * *CONSERVE_PRIOR_B*: the scale parameter for the gamma prior of conserved rate (default: 0.2).
  * *ACCE_PRIOR_A*: the shape parameter for the gamma prior of acceleraed rate (default: 5).
  * *ACCE_PRIOR_B*: the scale parameter for the gamma prior of accelerated rate (default: 0.04).

* **Control for MCMC**: 
  * *BURNIN*: number of initial iterations to discard before equilibrium of the chain (default: 100). Should set to be larger.
  * *MCMC*: number of MCMC iterations (default: 400). Should set to be larger. 
  * *SEED*: seed for random sampling (default: 5)
  * *RATIO0*: initial value for the conserved rate (default: 0.5). 
  * *RATIO1*: initial value for the neutral rate (default: 1). 
  * *CHAIN*: (default: 1)
