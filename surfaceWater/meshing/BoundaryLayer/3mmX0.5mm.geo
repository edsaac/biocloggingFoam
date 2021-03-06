/*
Purpose: Mesh to be used in a LES/RANS flow simulation
	of surface water over dunes following dimensions in
	Fox et al. 2018 (WRR)
*/

// Options for export
Mesh.SaveAll = 0;
Mesh.MshFileVersion = 2.2;

convertToMeter = 1.0/100;

//#### Dimensons ####################

Hd = 1.50 * convertToMeter;
Yd = 9.75 * convertToMeter;
Lee = 5.00 * convertToMeter;
Stoss = 10.00 * convertToMeter;
Width = 29.00 * convertToMeter;

//### Mesh Size #####################
meshSizeX = 0.30 * convertToMeter;
meshSizeY = 0.30 * convertToMeter;


nElementsX = Ceil((Lee+Stoss)/meshSizeX);
nElementsY = Ceil(Yd/meshSizeY);

// This growth ratio gives a del0 = 0.05mm
progY = 1.08953673902451; 

// Dummy value
meshSize  = 1.0;

Point(1) = {0,0,0,meshSize};
Point(2) = {Stoss/2,Hd/2,0,meshSize};
Point(3) = {Stoss/2+Lee,-Hd/2,0,meshSize};
Point(4) = {Stoss+Lee,0,0,meshSize};

Point(5) = {Stoss+Lee,Yd,0,meshSize};
Point(6) = {Stoss/2+Lee,Yd,0,meshSize};
Point(7) = {Stoss/2,Yd,0,meshSize};
Point(8) = {0,Yd,0,meshSize};

Spline(1) = {1,2,3,4};
Line(2) = {4,5};
Line(3) = {5,6,7,8};
Line(4) = {8,1};

Transfinite Line{1} = nElementsX;
Transfinite Line{3} = nElementsX;

Transfinite Line{2} = nElementsY Using Progression progY;
Transfinite Line{-4} = nElementsY Using Progression progY;

Line Loop(1) = {1 ... 4};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1,4,5,8};
Recombine Surface {1};
extrusionName[] = Extrude {0, 0, Width} {
    Surface{1};
    Layers{1};
    Recombine;
  };

Physical Surface("back") = {1};
Physical Surface("front") = extrusionName[0];
Physical Volume("bed") = extrusionName[1];
Physical Surface("bottom") = extrusionName[2];
Physical Surface("right") = extrusionName[3];
Physical Surface("top") = extrusionName[4];
Physical Surface("left") = extrusionName[5];
