# FourierVisualizer
A program to create cool visuals to help understand the fourier transform

Inspired by SmarterEveryDay's awesome video (https://www.youtube.com/watch?v=ds0cmAV-Yek), I have attempted to recreate 
the same effect in the video. 

The program is fairly simple. It takes in a path as a .ppm image file which it then extracts the path and takes
the fourier transform of the points. Then it visualizes the reverse fourier transform by drawing out every circle
that each complex exponential represents.

To use the program you should build with the provided Makefile and do ./ftest <image_file>. The image file is in a rather 
specific format, the trace should be 1px thick to give optimal results though it may work (not guaranteed) otherwise. I
plan to make a version using FFTW as it is much faster, however the bottleneck of this program is really in the reverse
transform where it would be far more difficult to use an FFT algorithm as it draws the results out iteratively. It also
uses OpenCV to save the video file, so make sure you install that with your favorite package manager before building & running.

Beyond that, I hope you enjoy!
