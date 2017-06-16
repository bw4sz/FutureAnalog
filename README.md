Future geographic patterns of novel and disappearing assemblages across three dimensions of diversity: a case study with Ecuadorian hummingbirds
============

__Authors:__ Laura J. Graham, Ben Weinstein, Sarah R. Supp, Catherine Graham

Code for [Future geographic patterns of novel and disappearing assemblages across three dimensions of diversity: a case study with Ecuadorian hummingbirds](http://onlinelibrary.wiley.com/doi/10.1111/ddi.12587/full) DOI: 10.1111/ddi.12587

License: This code is available under a BSD 2-Clause License.

Copyright (c) 2017. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
in the documentation and/or other materials provided with the distribution. 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact information Laura J. Graham's email: laurajaneegraham@gmail.com, Ben Weinstein's email: benweinstein2010@gmail.com

Contents
=================================

Data
--------------
`Inputdata` folder contains hummingbird localities for Ecuador and Colombia, the shapefile for the country perimeter of Ecuador
morphology, and phylogeny. 

Figures
--------------
Output results and figures are stored in the folder, with the exception of the large data and .Rdata files which are 
held on dropbox (please contact authors)

Requirements
---------------
R version 2.XX or greater and 
__Packages:__

`vegan`, `picante`, `analogue`, `doSNOW`, `ape`, `cluster`, `RColorBrewer`, `raster`, `ggplot2`, `phylobase`, `rgdal`, `tidyr`, `stringr`, `dplyr`, `biomod2`, `rasterVis`, `grid`, `devtools`, `broom`, `gridExtra`, `proxy`, `geometry`, `rcdd`, `cowplot`, `hypervolume`

## Workflow
Top level script `runscript.R` runs all code required for analysis:

- Installs and loads all required packages
- Sets resolution of analysis
- Sets up output folder structure
- Calls `runSDM()` function from `runSDM.R` (depends `fnSDM.R`)
- Calls `runProjections()` function from `runProjections.R`
- Calls `runBetaDiv()` function from `FutureAnalog.R` (depends `AlphaMappingFunctions.R`)
- Calls `runAnalogAnalysis` function from `FutureAnalog.R` for each threshold value (depends `FutureAnalogFunctions.R`)
- Calls `runFAPlots.R`
