#pragma once

#include <iostream>
#include "DFN.hpp"

using namespace std;

namespace DFNLibrary{

bool stampaDatiSulFileFrattureParaview(const string &percorsoFileFrattureParaview, Frattura &Fratt);

bool stampaDatiSulFileTracceParaview(const string &percorsoFileTracceParaview, Traccia &Trac);

}
