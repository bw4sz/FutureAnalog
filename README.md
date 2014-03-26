Hummingbird Taxonomic, Phylogenetic and Trait Analogs under Climate Change
============

__Authors:__ Ben Weinstein, Sarah R. Supp, Anusha Shankar, Marisa Lim, Catherine Graham

Code for computing analog and non-analog  hummingbird communities under future climate change scenarios.

License: This code is available under a BSD 2-Clause License.

Copyright (c) 2013. All rights reserved.

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

Contact information Ben Weinstein's email: benweinstein2010@gmail.com

Contents
=================================

`Inputdata` folder contains hummingbird localities for Ecuador and Colombia, the shapefile for the country perimeter of Ecuador
morphology, and phylogeny. 

Alpha Diversity
---------------
* `AlphaMapping.R` : Computes Alpha mapping of phylogenetic, taxonomic, and trait metrics
calls 
* `AlphaMappingFunctions.R` to calculate grid based metrics of alpha diversity on grid cells
* `SDM.R` which performs the ensemble niche models using biomod
-The input localities are in the Inputdata folder, but the env layers need to be held locally, they are > 30GB
* `BenHolttraitDiversity.R` to calculate MNND

Beta Diversity
------------
* `FutureAnalog.R` computes betametrics between current and future assemblages, needs to be done carefully, 
consider how parallelization works on your operating system. 
* `GDM_hbird.R` computes beta diversity using a generalized dissimilarity modeling (GDM) framework.

Figures
--------------
Output results and figures are stored in the folder, with the exception of the large data and .Rdata files which are 
held on dropbox (please contact authors)

Requirements
---------------
R version 2.XX or greater and 
__Packages:__

* analogue 
* ape
* biomod2
* doSNOW
* ecodist
* FD
* ggplot2
* maptools
* MASS
* parallel
* picante
* raster
* RColorBrewer
* reshape
* rgdal
* stringr
* vegan
