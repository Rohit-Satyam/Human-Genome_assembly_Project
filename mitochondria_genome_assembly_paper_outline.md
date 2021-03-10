# Mitochondrial Genome Assembly Paper Outline

10 checkpoints from [Essential guidelines for computational method benchmarking](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8)

### 1. Define the purpose and scope of the benchmark.
> Benchmarking tools for accurate mitochondrial geneome assembly from Short read WGS data.

### 2. Include all relevant methods.
> Compilation of list of assembly suite involved in Mamalian mitochondrial genome assembly/ Organelle assembly tools. Here is a comprehensive list of tools collwcted so far
>MitoZ
Norgal
mitoMaker
SMART
MITObim
NOVOPlasty
GetOrganelle
MEANGS
IOGA
Org.ASM
>TODO: is there another one we are missing?

### 3. Select (or design) representative datasets.
> We plan to use simulated dataset and real datasets (Taking Mother Father and Child-TRIO of two different sequencing read lengths). Alternately, datasets can also be taken from 1000 genomes project (NA11920, HG01112, NA18941, HG00096,
HG00273, NA18548, NA18510) to account for differences in sequencing depth. Also, adding an ancient Genome Sequence would be better from the [2014 study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4269527/). TODO: select exact list of chloros TODO: produce simulated datasets
>The choice of simulator will be based on the following literature: [Source1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5224698/) or [Source2](https://www.nature.com/articles/nrg.2016.57).

### 4. Choose appropriate parameter values and software versions.
>Documenting the source and installation procedure of each tool. Will try to make docker images for each tool if possible. Try keeping default parameters as possible 
>TODO: select default parameters for each tool

### 5. Evaluate and rank methods according to key quantitative performance metrics.
> Quantitative Metrices: Using the Whole Genome Assembly Quality Assessment tools. Use sanger mitochondrial genome assembly from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/251831106) and CHM13 v0.1 for accuracy assessment.
> Qualitative metrics (success, failure, incomplete, ...) TODO: write script to gather these metrics from output

### 6. Evaluate secondary measures including runtimes and computational requirements,user-friendliness, code quality, and documentation quality.
> we have a script to track all performence metrics (time consumption and peak memory). Will try to have three replicate runs for each thread chosen. 
>TODO: separate performence benchmarking runs

### 7. Interpret results and provide guidelines or recommendations from both user and method developer perspectives.
>How can assembling mitochondrial genomes fast can leverage clinical inferences about associated diseases. 
>How a chimeric Genome assembly (graph based mitochondrial genome) can be instrumental, 
>caveats on the basis of present benchmarking observed and room for improvement (at algorithmic levels)

### 8. Publish and distribute results in an accessible format.
GitHub, zenodo, DockerHub, biorXiv, BMC TODO

### 9. Design the benchmark to enable future extensions.
> Strategy 1 (Effective When Taking both 1000 genome project dataset and account for sequencing depth): Convert BAMs obtained from GIAB to FASTQ and mapping them to revised Cambridge mitochondrial reference sequence (rCRS) which has not changed since 1999 [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/251831106) and filtering out reads. 
> Strategy2 (Effective when taking only GIAB dataset of and do not wish to account for sequencing depth): Directly using samtools retrive Mitochondrial sequences. 
TODO: documentation on GitHub on how to reproduce the benchmarking (incl. extension)

### 10. Follow reproducible research best practices, in particular by making all code and data publicly available.
>Already covered with all previous points
