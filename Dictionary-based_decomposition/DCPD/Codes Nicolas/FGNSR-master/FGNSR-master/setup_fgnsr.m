function setup_fgnsr
% Build fast projection mex file from matlab and setup path

mex -largeArrayDims -outdir matlab -O ckernel/parproj_priv.c ckernel/proj.c ckernel/util.c
addpath(fullfile(pwd, 'matlab'));

example_midpoints
disp('FGNSR setup done, issue ''savepath'' in order to save your path.');
end
