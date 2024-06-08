#pragma once

#include <iostream>
#include "DFN.hpp"

using namespace std;

namespace DFNLibrary{

bool stampaDatiSulFileFrattureParaview(const string &percorsoFileFrattureParaview, DFN &Fract);

bool stampaDatiSulFileTracceParaview(const string &percorsoFileTracceParaview, DFN &Fract);

}
