// #pragma once

// #include <iostream>
// #include "Eigen/Eigen"

// using namespace std;
// using namespace Eigen;

// namespace PolygonalLibrary{

//     struct PolygonalMesh{

//         unsigned int NumberCell0Ds = 0; // numero delle Cell0Ds
//         vector<unsigned int> IdCell0Ds = {}; // Cell0Ds id, size 1 x NumberCell0Ds
//         vector<Vector3d> CoordinatesCell0Ds = {}; // coordinate Cell0Ds, size 3 x NumberCell0Ds (x,y)

//         unsigned int NumberCell1Ds = 0; // numero delle Cell1Ds
//         vector<unsigned int> IdCell1Ds = {}; // Cell1Ds id, size 1 x NumberCell1Ds
//         vector<Vector2i> VerticesCell1Ds = {}; // indici vertici Cell1Ds, size 2 x NumberCell1Ds (Origin,End) --> (DaId, All'Id)

//         unsigned int NumberCell2Ds = 0; // numero delle Cell2Ds
//         vector<unsigned int> NumVertices = {};
//         vector<unsigned int> NumEdges = {};
//         vector<unsigned int> IdCell2Ds = {}; // Cell2Ds id, size 1 x NumberCell2Ds
//         vector<vector<unsigned int>> VerticesCell2Ds = {}; // indici vertici Cell2Ds, size 1 x NumberCell2DsVertices[NumberCell2Ds]
//         vector<vector<unsigned int>> EdgesCell2Ds = {}; // indici Cell1Ds Cell2Ds, size 1 x NumberCell2DsEdges[NumberCell2Ds]
//     };
// }
