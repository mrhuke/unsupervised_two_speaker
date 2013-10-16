Unsupervised Two Speaker Separation
===================================

The programs corresponds to the algorithm described in "An unsupervised approach to cochannel speech separation" by 
K. Hu and D. L. Wang (IEEE Trans. Audio, Speech, and Lang. Process., 2013). This is an unsupervised
algorithm for two-speaker separation.


Requirements:

- The input mixture (the wav file) needs to have a sampling frequency of 16 kHz

- The algorithms is developed in MATLAB under Linux

- The main MATLAB program will call external executables compiled in Linux version 2.6.32-358.2.1.el6.x86_64 in a REHL 6.4 distribution. In other systems, you have to compile the tandem algorithm (source code in folder "tandem") and segmentation algorithm (source code in folder "segment") and generate your own executables.

- The tandem algorithm used here is not the most unpdated version. It does not generate and group T-segments. The newest version can be found at github.com/mrhuke/tandem

Run an example:

1. Under a Linux system, go to the "run" folder
2. start MATLAB
3. mixture = load('mixture');
4. mask = twoSpk_unsupervised(mixture);


A description of the main steps performed:

1. Run a tandem algorithm (Hu & Wang'11) to generate simultaneous streams (SS)
2. Order SS by time
3. Extract GFCCs for each SS
4. Group voiced SS by beam search
5. Generate unvoiced speech segments by onset/offset based segmentation
6. Group unvoiced-voiced and unvoiced-unvoiced segments
