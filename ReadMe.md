# PCF - Point Cloud Fit

```txt
     ______  ______         ______  __  _______   
    /  _  / / ____/  ___   / ____/ / / /__  __/   
   / ____/ / /___   /__/  / ___/  / /    / /      
  /_/     /_____/        /_/     /_/    /_/       
```

PCF - Point Cloud Fit, is a robust automatic detection scheme for geometric primitives such as planes, cylinders and cuboids.

In this scheme, cylinders are first detected in iteration of energy-based geometric model fitting (by [GCO](http://vision.csd.uwo.ca/code/)) and cylinder parameter estimation (by [GTEngine](https://www.geometrictools.com/Downloads/Downloads.html)). Planes are then detected by Hough transform, and further described as bounded patches with their minimum bounding rectangles. Cuboids are finally detected with pair-wise geometry relations from the detected patches. After successively detection of cylinders, planar patches and cuboids, a mid-level geometry representation can be delivered.

> The detection results can be observed by rendering in 3D model view tool [PlyWin](https://github.com/weiquanmao/PlyWin).

## Projects

+ Obj2Ply : sampling the mesh model data (`*.obj`) into point cloud model data (`*.ply`).

+ PCE : the robust automatic detection scheme. This scheme takes the point cloud model `*.ply` as input and finally output the `*.strcut` file in which the detected components are recorded.

+ Struct2Ply : sampling the the detected components (`*.strcut`) into point cloud model data (`*.ply`).


The accuracy and completeness estimations for 3D reconstruction can be utilized here to quantitatively evaluate the accuracy of the detection results. Project **Struct2Ply** can offer you a convenience to convert the detection results into point cloud models, and the accuracy, completeness and F-score can then be estimated with [PCE](https://github.com/weiquanmao/PCE).

## Dependencies

+ [Qt](www.qt.io/download/)
+ [VCG Library](http://vcg.isti.cnr.it/vcglib/) (Embedded)
+ [GCO](http://vision.csd.uwo.ca/code/) (Embedded)
+ [GTEngine](https://www.geometrictools.com/Downloads/Downloads.html) (Embedded)

## Examples

Detection results of the test data:

![Results](https://github.com/weiquanmao/PCF/blob/master/TestData/Result.jpg)

## Publications

This program and data in `TestData` have been used in the following publications:

```
@article{weiSensors18,
  author  = {Quanmao Wei and Zhiguo Jiang and Haopeng Zhang},
  title   = {Robust Spacecraft Component Detection in Point Clouds},
  journal = {Sensors},
  volume  = {18},
  year    = {2018},
  number  = {4},
  article number = {933},
  url     = {http://www.mdpi.com/1424-8220/18/4/933},
  issn    = {1424-8220},
  doi     = {10.3390/s18040933}
}
```
```
@inproceedings{weiIGTA17,
  author    = {Quanmao Wei and Zhiguo Jiang and Haopeng Zhang and Shanlan Nie},
  title     = {Spacecraft Component Detection in Point Clouds},
  booktitle = {Advances in Image and Graphics Technologies},
  year      = {2017},
  publisher = {Springer Singapore},
  address   = {Singapore},
  pages     = {210--218},
  isbn      = {978-981-10-7389-2}
}
```

---
By [WeiQM](https://weiquanmao.github.io) at D409.IPC.BUAA