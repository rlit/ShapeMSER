ShapeMSER v1.0 (code rev 258)

- This code packeage implements the algorithm of MSER on deformable shapes
- Original paper in - http://dx.doi.org/10.1016/j.cag.2011.03.011
- Perform "SVN checkout" of all subfolders - currently there are 2 main folders: "code" and "data". 
- The code was tested on matlab 2010b (version 7.11.0) on a windows platform (should contain mex files for both 32 & 64 bit architecture).
- Before running anything please run "InitPath.m" in the code folder to add the necceassary folders to the matlab path.
- If you use this package for your publication, please cite this work using the "CiteME.bib" file.

to make sure the framework is working
- (without ground-truth) run "test_centaur4.m" 
- (with ground-truth)    run "shrec2010_bench.m" 

NOTE:
- This package does not contain the datasets necessary to run the benchmark - SHREC 2010 and SHREC 2011, nor the ground-truth correspondances.
- data may be acquired from http://tosca.cs.technion.ac.il/book/shrec.html (or by google-ing e.g. "SHREC 2010")
- for the full ground-truth please contact the author
