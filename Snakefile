#!/usr/bin/env python3

import os
import pandas
import pathlib
import pickle
import re





def resolve_path(x):
    return(str(pathlib.Path(x).resolve()))

#########
# RULES #
#########

expand('output/022_fastqc/{individual}_fastqc.zip',
               individual=all_indivs)












