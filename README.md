## Argentine stem weevil GBS

An honours project analysing genotyping by sequencing data from the Argentine stem weevil (*Listronotus bonariensis*).

Running the first four snakefiles in order will allow the methods for this honours project to be reproduced. The remaining snakefiles are for a subsequent extension of the project. In order to reproduce the methods, use:

``snakemake -s <snakefile.name> --use-singularity --cores=<number> ``

Snakemake (version 5.1.2) and Singularity (version 2.5.1) are required.

The singularity containers used are hosted on singularity hub at:

https://www.singularity-hub.org/collections/1290
and
https://singularity-hub.org/collections/996

Plotting scripts are in src/plots and are not run automatically.
