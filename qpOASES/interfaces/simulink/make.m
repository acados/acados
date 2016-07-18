function [] = make( varargin )
%MAKE Compiles the Simulink interface of qpOASES_e.
%
%Type  make            to compile all interfaces that 
%                      have been modified,
%type  make clean      to delete all compiled interfaces,
%type  make clean all  to first delete and then compile 
%                      all interfaces,
%type  make 'name'     to compile only the interface with
%                      the given name (if it has been modified),
%type  make 'opt'      to compile all interfaces using the 
%                      given compiler options
%
%Copyright (C) 2013-2015 by Hans Joachim Ferreau, Andreas Potschka,
%Christian Kirches et al. All rights reserved.

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
%%	Filename:  interfaces/simulink/make.m
%%	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
%%	Version:   3.1embedded
%%	Date:      2007-2015
%%

       
    %% consistency check
    if ( exist( [pwd, '/make.m'],'file' ) == 0 )
        error( ['ERROR (',mfilename '.m): Run this make script directly within the directory', ...
                '<qpOASES-inst-dir>/interfaces/simulink, please.'] );
    end


    if ( nargin > 2 )
        error( ['ERROR (',mfilename '.m): At most two make arguments supported!'] );
    else
        [ doClean,fcnNames,userFlags ] = analyseMakeArguments( nargin,varargin );
    end

    
    %% define compiler settings
    QPOASESPATH = '../../';
    
    DEBUGFLAGS = ' ';
    %DEBUGFLAGS = ' -g CXXDEBUGFLAGS=''$CXXDEBUGFLAGS -Wall -pedantic -Wshadow'' ';

    IFLAGS = [ '-I. -I',QPOASESPATH,'include',' -I',QPOASESPATH,'src',' ' ];
    CFLAGS = [ IFLAGS, DEBUGFLAGS, '-largeArrayDims -D__MATLAB__ -D__SINGLE_OBJECT__ -Dinline="" -Dsnprintf="_snprintf"',' ' ];
    defaultFlags = '-O -D__NO_COPYRIGHT__ '; %% -D__NO_COPYRIGHT__ -D__SUPPRESSANYOUTPUT__ -D__MANY_CONSTRAINTS__ 

    if ( ispc == 0 )
        CFLAGS  = [ CFLAGS, '-DLINUX ',' ' ]; 
    else
        CFLAGS  = [ CFLAGS, '-DWIN32 ',' ' ];
    end

    if ( isempty(userFlags) > 0 )
        CFLAGS = [ CFLAGS, defaultFlags,' ' ];
    else
        CFLAGS = [ CFLAGS, userFlags,' ' ];
    end

    mexExt = eval('mexext');
    
  
	%% ensure copyright notice is displayed
    if ~isempty( strfind( CFLAGS,'-D__NO_COPYRIGHT__' ) )
        printCopyrightNotice( );
    end
    
  
    %% clean if desired
    if ( doClean > 0 )
        
        eval( 'delete *.o;' );
        eval( ['delete *.',mexExt,'*;'] );
        disp( [ 'INFO (',mfilename '.m): Cleaned all compiled files.'] );
        pause( 0.2 );
        
    end
    
    
    if ( ~isempty(userFlags) )
        disp( [ 'INFO (',mfilename '.m): Compiling all files with user-defined compiler flags (''',userFlags,''')...'] );
    end
    
    
    %% call mex compiler
    for ii=1:length(fcnNames)
        
        cmd = [ 'mex -output ', fcnNames{ii}, ' ', CFLAGS, [fcnNames{ii},'.c'] ];
        
        if ( exist( [fcnNames{ii},'.',mexExt],'file' ) == 0 )
            
            eval( cmd );
            disp( [ 'INFO (',mfilename '.m): ', fcnNames{ii},'.',mexExt, ' successfully created.'] );
            
        else
            
            % check modification time of source/Make files and compiled mex file
            cppFile = dir( [pwd,'/',fcnNames{ii},'.c'] );
            cppFileTimestamp = getTimestamp( cppFile );
                       
            makeFile = dir( [pwd,'/make.m'] );
            makeFileTimestamp = getTimestamp( makeFile );
            
            mexFile = dir( [pwd,'/',fcnNames{ii},'.',mexExt] );
            if ( isempty(mexFile) == 0 )
                mexFileTimestamp = getTimestamp( mexFile );
            else
                mexFileTimestamp = 0;
            end
            
            if ( ( cppFileTimestamp   >= mexFileTimestamp ) || ...
                 ( makeFileTimestamp  >= mexFileTimestamp ) )
                eval( cmd );
                disp( [ 'INFO (',mfilename '.m): ', fcnNames{ii},'.',mexExt, ' successfully created.'] );
            else            
                disp( [ 'INFO (',mfilename '.m): ', fcnNames{ii},'.',mexExt, ' already exists.'] );
            end
            
        end
        
    end

    %% add qpOASES directory to path
    path( path,pwd );

end


function [ doClean,fcnNames,userIFlags ] = analyseMakeArguments( nArgs,args )

    doClean = 0;
    fcnNames = [];
    userIFlags = [];

    switch ( nArgs )
        
        case 1
            if ( strcmp( args{1},'all' ) > 0 )
                fcnNames = { 'qpOASES_e_QProblemB','qpOASES_e_QProblem' };
            elseif ( strcmp( args{1},'qpOASES_e_QProblemB' ) > 0 )
                fcnNames = { 'qpOASES_e_QProblemB' };
            elseif ( strcmp( args{1},'qpOASES_e_QProblem' ) > 0 )
                fcnNames = { 'qpOASES_e_QProblem' };
            elseif ( strcmp( args{1},'clean' ) > 0 )
                doClean = 1;
            elseif ( strcmp( args{1}(1),'-' ) > 0 )
                % make clean all with user-specified compiler flags
                userIFlags = args{1};
                doClean = 1;
                fcnNames = { 'qpOASES_e_QProblemB','qpOASES_e_QProblem' };
            else
                error( ['ERROR (',mfilename '.m): Invalid first argument (''',args{1},''')!'] );
            end

        case 2
            if ( strcmp( args{1},'clean' ) > 0 )
                doClean = 1;
            else
                error( ['ERROR (',mfilename '.m): First argument must be ''clean'' if two arguments are provided!'] );
            end
            
            if ( strcmp( args{2},'all' ) > 0 )
                fcnNames = { 'qpOASES_e_QProblemB','qpOASES_e_QProblem' };
            elseif ( strcmp( args{2},'qpOASES_e_QProblemB' ) > 0 )
                fcnNames = { 'qpOASES_e_QProblemB' };
            elseif ( strcmp( args{2},'qpOASES_e_QProblem' ) > 0 )
                fcnNames = { 'qpOASES_e_QProblem' };
            else
                error( ['ERROR (',mfilename '.m): Invalid second argument (''',args{2},''')!'] );
            end
            
        otherwise
            doClean = 0;
            fcnNames = { 'qpOASES_e_QProblemB','qpOASES_e_QProblem' };
            userIFlags = [];            
    end
    
end


function [ timestamp ] = getTimestamp( dateString )

    try 
        timestamp = dateString.datenum;
    catch
        timestamp = Inf;
    end

end


function [ ] = printCopyrightNotice( )

    disp( ' ' );
    disp( 'qpOASES -- An Implementation of the Online Active Set Strategy.' );
    disp( 'Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,' );
    disp( 'Christian Kirches et al. All rights reserved.' );
    disp( ' ' );
    disp( 'qpOASES is distributed under the terms of the' );
    disp( 'GNU Lesser General Public License 2.1 in the hope that it will be' );
    disp( 'useful, but WITHOUT ANY WARRANTY; without even the implied warranty' );
    disp( 'of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.' );
    disp( 'See the GNU Lesser General Public License for more details.' );
    disp( ' ' );
    disp( ' ' );

end


%%
%%	end of file
%%
