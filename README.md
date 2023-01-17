# AUDDIT

Version 2
Creators: Andrew Clark and Zach Kalmanson
San Miguel Lab, NCSU College of Chemical and Biomolecular Engineering 

Installation instructions: download the AUDDIT installer .exe and run the application. No additional files are necessary to install it. 

Using AUDDIT 
---------------
1. Input: enter the pixel size of the camera used for taking the images in microns. Select the folder you want to read the images from and the folder you want to save the results to in the respective file directory and save directory options. All images used must be .tif or .tiff file types. Each subfolder of images in the folder you select will be treated as a experimental group. 
2. Press run algorithm to begin the analysis. A window will pop up when the analysis is complete. 
3. You can look at the analysis of the images in the folder you selected as the save directory. If there are any images that were incorrectly analyzed by the algorithm, you can remove those data points in this step. 
4. Check all of the boxes of properties you would like to observe. After pressing "Print Results", figures will be generated and sent to the save directory folder. 

Additional Tips
---------------
To ensure the algorithm can accurately analyze neuron images, it is important that certain practices are folowed. 
All images in each experimental group should be taken with the same imaging set up and at the same magnification. 
Try to ensure that the head region of the worm is in the center of the frame. 
Each image should have only one worm in the frame. Additional worms or objects in the frame will result in inaccurate analysis. 
When looking at projections of z-stack images of the CEP dendrites, it is common for the worm to be oriented such that one or more of the dendrites are overlapping with each other. Any images with overlapping dendrites are not able to be analyzed and should be discarded. 
Please see the supplemental of the associated paper for examples of usable and unusuable images. 
