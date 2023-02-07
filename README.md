# hilbertgenome-py
Create Hilbert curve based data structures for visualization of genomic data.

The human genome contains roughly 3 billion basepairs. There are many research datasets that cover vast stretches of the genome with useful information. Visualizing these datasets is a challenge given the scale of potentially billions of data points.

This library helps process BED files into an efficient file format for use in web-based data visualizations that allow interactively zooming out to see the aggregations across the entire genome down to data associated with individual basepairs.

## Usage
The `hilbertgenome` python module...

```py
# TODO: example code
```

TODO: Jupyter Notebook link



## Examples
TODO: altius app

TODO: simple observable notebook demo


## Why Hilbert curves?
Hilbert curves have a couple of nice properties for allowing us to visualize the genome. 

The first is that by laying out 1D points on a 2D hilbert curve, the points close together in 1D will be relatively close together in 2D (and vice versa). 
This means we can take a very long list of points (like the basepairs in the reference human genome) and lay them out in a 2D plane. 
Placing the genome on a 2D plane allows us to treat it as a map, which we can then design with cartographic and data visualization principles.  

It also means that given some zoomed in region we can efficiently query a 1D range of numbers within the 2D bounding box. This allows us to store a very simple datastructure that can be quickly and cheaply queried.

Another nice property is that each order of the Hilbert curve is 4x more points than the order before, and the extra points are easily addressable. So we can compute 2D aggregations at various zoom levels quickly and efficiently as well.

## Hilbert indexing explained
TODO: a few images showing mappings

TODO: observable notebook with visuals