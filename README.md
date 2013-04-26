unsupervised_two_speaker
========================

Unsupervised Two Speaker Separation

This program implements the algorithm in "An unsupervised approach to cochannel speech separation" by 
K. Hu and D. L. Wang (submited to IEEE Trans. Audio, Speech, and Lang. Process.). This is an unsupervised
algorithm for two-speaker separation.

In a Linux system, go to "run" folder and execute seqGrp.m to start the program. It requires a 16-kHz time-domain cochannel speech as the input.
To run a sample program, do the following:
- Under a Linux system, go to the "run" folder
- start MATLAB
- mixture = load('mixture');
- mask = twoSpk_unsupervised(mixture);


Following is a text description of the main steps in seqGrp.m:

1. Run a tandem algorithm (Hu & Wang'11) to estimate simultaneous streams

2. Rank simultaneous streams by time

3. Extract GFCC features (Shao'07) for each simultaneous stream

4. Perform beam search to group voiced simultaneous streams

5. Onset/offset segmentation to generate unvoiced segments

6. Group unvoiced-voiced and unvoiced-unvoiced segments
