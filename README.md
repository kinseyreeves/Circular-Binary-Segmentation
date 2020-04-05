# Circular-Binary-Segmentation


Python implementation of circular binary segmentation in Olshen et als Paper. Segments copy number changes into different intervals along the dataset.
Recursively cuts segments until no segment can find a Z value greater than the threshold.

First, the CBS algorithm splices segments and then recursively searches them for further segments with a sufficient Z value. If it does not find one, it will output this whole segment into a staging array. Secondly, from this array segments are discarded if their average absolute log ratio value is less than 0.1. Finally, we must check for contiguous segments in the output as some segments are circularised and therefore spliced out from others.  

Paper : Olshen AB, Circular binary segmentation for the analysis of array-based DNA copy number data. Department of Epidemiology and Biostatistics

Usage:

`python cbs.py <input file> <z-threshold> <outputfile>`

