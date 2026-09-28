#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# author: Daniel Kloeser

import numpy as np
import os
from pathlib import Path

def getTrack(filename):
    track_file = os.path.join(str(Path(__file__).parent), filename)
    array=np.loadtxt(track_file)
    sref=array[:,0]
    xref=array[:,1]
    yref=array[:,2]
    psiref=array[:,3]
    kapparef=array[:,4]
    return sref,xref,yref,psiref,kapparef
