# Multiple Sequence Alignment

The object of this python code is multiply align three sequences using a 3-D Manhattan Cube with each axis representing a sequence.
<p align="center">
<img src="https://github.com/Sulstice/Bioinformatics/blob/master/Screen%20Shot%202018-01-06%20at%203.27.28%20PM.png" width="310">
 </p>

The edges of the cube are 7 and thus can be represented mathematically like so

Geometric View             |  Mathematical View
:-------------------------:|:-------------------------:
![](https://github.com/Sulstice/Bioinformatics/blob/master/Screen%20Shot%202018-01-06%20at%203.28.33%20PM.png)  |  ![](https://github.com/Sulstice/Bioinformatics/blob/master/Screen%20Shot%202018-01-06%20at%203.28.54%20PM.png)

The run time of the algorithm is **O(n^3)**, which did not scale well with the language of python. Initially, the matrix data was stored in lists of lists which proved to be slow with the algorithm by processing one triple base alignment every 2-3 seconds. In short sequences this can prove to be viable since it is easy code to write and learn from. With sequences upwards of 1000 it would take 24 hours. 

<p align="center">
<img src="https://github.com/Sulstice/Bioinformatics/blob/master/Screen%20Shot%202018-01-06%20at%203.34.02%20PM.png" width="310">
 </p>
 
 To increase processing time, in this code is implemented a time-saving recursion algorithm to approach the center of the cube from two different entry points and ultimately align the sequences together. 
 
 Although this saved time, due to the algorithm's current nature it still did not resolve the scaling issue and python parallel processing can be more than troublesome to implement. 
