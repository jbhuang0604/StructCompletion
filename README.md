StructCompletion
================

Image Completion using Planar Structure Guidance

This is a Matlab re-implementation of the paper.

Jia-Bin Huang, Sing Bing Kang, Narendra Ahuja, and Johannes Kopf, 
Image Completion using Planar Structure Guidance,
ACM Transactions on Graphics (Proceedings of SIGGRAPH 2014), 33(4), 2014

It is provided for educational/research purpose only. Please cite our paper if you found the software useful for your work.

To run the code, please use the main function sc_complete.m.
 
===
Example
===

imgFileName = '005_input_hole.png';

imgCompletion = sc_complete(imgFileName);

===
Content
===

There are three main directories:

- cache: pre-computed data for vanishing point detection

- data: the folder for input images

- external: external codes

1. MeanShift: the Mean-shift clustering algorithm
http://www.mathworks.com/matlabcentral/fileexchange/10161-mean-shift-clustering

2. vlfeat-0.9.19: for feature detection and matching
http://www.vlfeat.org/
Please download the vlfeat library and unzip it to the 'external' folder

3. mirt2D_mexinterp: Fast 2D linear interpolation
http://www.mathworks.com/matlabcentral/fileexchange/24183-2d-interpolation/content/mirt2D_mexinterp/mirt2D_mexinterp.m

===
Contact 
===

Jia-Bin Huang

jbhuang1@illinois.edu

www.jiabinhuang.com
