#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolygonalLibrary{

    struct PolygonalMesh{

        unsigned int NumberCell0Ds = 0; // numero delle Cell0Ds
        vector<unsigned int> IdCell0Ds = {}; // Cell0Ds id, size 1 x NumberCell0Ds
        vector<Vector3d> CoordinatesCell0Ds = {}; // coordinate Cell0Ds, size 2 x NumberCell0Ds (x,y)
        map<unsigned int, list<unsigned int>> MarkersCell0Ds = {}; // Cell0Ds markers, size 1 x NumberCell0Ds (marker)

        unsigned int NumberCell1Ds = 0; // numero delle Cell1Ds
        vector<unsigned int> IdCell1Ds = {}; // Cell1Ds id, size 1 x NumberCell1Ds
        vector<Vector2i> VerticesCell1Ds = {}; // indici vertici Cell1Ds, size 2 x NumberCell1Ds (Origin,End) --> (DaId, All'Id)
        map<unsigned int, list<unsigned int>> MarkersCell1Ds = {}; // proprietà Cell1Ds, size 1 x NumberCell1Ds (marker)

        unsigned int NumberCell2Ds = 0; // numero delle Cell2Ds
        unsigned int NumVertices = 0;
        unsigned int NumEdges = 0;
        vector<unsigned int> IdCell2Ds = {}; // Cell2Ds id, size 1 x NumberCell2Ds
        vector<vector<unsigned int>> VerticesCell2Ds = {}; // indici vertici Cell2Ds, size 1 x NumberCell2DsVertices[NumberCell2Ds]
        vector<vector<unsigned int>> EdgesCell2Ds = {}; // indici Cell1Ds Cell2Ds, size 1 x NumberCell2DsEdges[NumberCell2Ds]
        map<unsigned int, list<unsigned int>> MarkersCell2Ds = {}; // proprietà Cell2Ds, size 1 x NumberCell2Ds (marker)
    };
}
