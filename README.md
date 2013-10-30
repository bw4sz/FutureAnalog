Hummingbird Taxonomic, Phylogenetic and Trait Analogs under Climate Change
============

Authors: Ben Weinstein, Sarah R. Supp, Anusha Shankar

Code for computing non-analog communities under future climate change.

License: This code is available under a BSD 2-Clause License.

Copyright (c) 2013. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact information Sarah Supp's email: sarah@weecology.org and sarah.supp@stonybrook.edu


Readme
=================================

_AlphaMapping.R_ : Computes Alpha mapping of phylogenetic, taxonomic, and trait metrics
calls _SDM.R_ which performs the ensemble niche models using biomod
-The input localities are in the Inputdata folder, but the env layers need to help locally, there are > 30GB

_FutureAnalog.R_ computes betametrics between current and future assemblages, needs to be done carefully, consider how parallelization works on your operating system. 
