%%
%%	This file is part of qpOASES.
%%
%%	qpOASES -- An Implementation of the Online Active Set Strategy.
%%	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
%%	Christian Kirches et al. All rights reserved.
%%
%%	qpOASES is free software; you can redistribute it and/or
%%	modify it under the terms of the GNU Lesser General Public
%%	License as published by the Free Software Foundation; either
%%	version 2.1 of the License, or (at your option) any later version.
%%
%%	qpOASES is distributed in the hope that it will be useful,
%%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%	See the GNU Lesser General Public License for more details.
%%
%%	You should have received a copy of the GNU Lesser General Public
%%	License along with qpOASES; if not, write to the Free Software
%%	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%%



%%
%%	Filename:  interfaces/simulink/load_example_QProblemB.m
%%	Author:    Hans Joachim Ferreau (thanks to Aude Perrin)
%%	Version:   3.1embedded
%%	Date:      2007-2015
%%



clear all;


%% setup QP data
simulationTime  =  [0;0.1];

H = [ 1.0,0.0; ...
      0.0,0.5 ];

g.time = simulationTime;
data1 = [ 1.5,1.0 ];
data2 = [ 1.0,1.5 ];
g.signals.values = [data1; data2];
g.signals.dimensions = numel(data1);

lb.time = simulationTime;
data1 = [ 0.5,-2.0 ];
data2 = [ 0.0,-1.0 ];
lb.signals.values = [data1; data2];
lb.signals.dimensions = numel(data1);

ub.time = simulationTime;
data1 = [ 5.0, 2.0 ];
data2 = [ 5.0,-0.5 ];
ub.signals.values = [data1; data2];
ub.signals.dimensions = numel(data1);

clear simulationTime data1 data2


%% open corresponding simulink example
open( 'example_QProblemB.mdl' );



%%
%%	end of file
%%
