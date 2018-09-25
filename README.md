## Argentine stem weevil GBS

An honours project analysing genotyping by sequencing data from the Argentine stem weevil (*Listronotus bonariensis*).

Running the four snakefiles in order will allow the methods to be reproduced. In order to do so, use:

``snakemake -s <snakefile.name> --use-singularity --cores=<number> ``

Snakemake (version 5.1.2) and Singularity (version 2.5.1) are required.

The singularity containers used are hosted on singularity hub at:

https://www.singularity-hub.org/collections/1290
and
https://singularity-hub.org/collections/996

Plotting scripts are in src/plots.
