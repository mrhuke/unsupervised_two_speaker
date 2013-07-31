Unsupervised Two Speaker Separation
===================================

The programs corresponds to the algorithm described in "An unsupervised approach to cochannel speech separation" by 
K. Hu and D. L. Wang (IEEE Trans. Audio, Speech, and Lang. Process., 2013). This is an unsupervised
algorithm for two-speaker separation.

- System requirements:
It is recommended to run the algorithms under a Linux system in MATLAB. The main program will call external executable files which are already compiled in Linux version 2.6.32-358.2.1.el6.x86_64 in a REHL 6.4 distribution. Otherwise you have to compile the tandem algorithm (source code in folder "tandem") and segmentation algorithm (source code in folder "segment") and generate your own executables.

- Run an example:
In a Linux system, go to "run" folder and execute seqGrp.m to start the program. It requires a 16-kHz time-domain cochannel speech as the input.
To run a sample program, do the following:
- Under a Linux system, go to the "run" folder
- start MATLAB
- mixture = load('mixture');
- mask = twoSpk_unsupervised(mixture);


- A description of the main steps performed:

1. Run a tandem algorithm (Hu & Wang'11) to estimate simultaneous streams
2. Rank simultaneous streams by time
3. Extract GFCC features (Shao'07) for each simultaneous stream
4. Perform beam search to group voiced simultaneous streams
5. Onset/offset segmentation to generate unvoiced segments
6. Group unvoiced-voiced and unvoiced-unvoiced segments
