GKW
======

GKW is a gyro-kinetic simulation code for the study of turbulence
in magnetised plasmas developed from the original linear code LINART.

  *Copyright (C) 2003, 2004, 2005*
    A.G. Peeters, D. Strintzi

  *Copyright (C) 2007, 2008, 2009, 2010, 2011*
    A.G. Peeters, Y. Camenen, F.J. Casson, W.A. Hornsby, A.P. Snodin,
    D. Strintzi, G. Szepesi
    
  *Copyright (C) 2012, 2013, 2014, 2015*
    A.G. Peeters, Y. Camenen, F.J. Casson, W.A. Hornsby, A.P. Snodin,
    D. Strintzi, G. Szepesi, R. Buchholz, S. Grosshauser, P. Manas, 
    P. Migliano, M. Siccinio, T. Sung,  D. Zarzoso  

This directory contains the GKW source code, documentation and
associated tools. This file contains a brief overview of the source code
project and package contents. Some files which constitute parts of the
code may be the works of (or derivative works of) other authors; notes
to this effect are included at the top of any such files. Please see the
end of this file for license conditions.

A detailed description of the code, how to build it, and how to run it
can be found in the manual, located in the doc/ directory, and in
the associated paper, http://dx.doi.org/10.1016/j.cpc.2009.07.001

If you use GKW (or some results obtained from it) in any publication, we
politely request that you cite paper **[1]** below, in which the code is
comprehensively described. You may also wish to cite the original LINART
paper **[2]** and/or the first paper in which GKW was used **[3]**.

 **[1]** A.G. Peeters Y. Camenen, F.J. Casson, W.A. Hornsby, A.P. Snodin,
     D. Strintzi, and G. Szepesi,
     Computer Physics Communications, 180, 2650 (2009)
     http://dx.doi.org/10.1016/j.cpc.2009.07.001

 **[2]** A.G. Peeters, D. Strintzi, Phys. Plasmas, 11, 3748 (2004)
     http://dx.doi.org/10.1063/1.1762876

 **[3]** A.G. Peeters, C. Angioni, D. Strintzi,
     Phys. Rev. Lett. 98, 265003 (2007) 
     http://dx.doi.org/10.1103/PhysRevLett.98.265003

Versions of the GKW code corresponding directly to the article **[1]** can
be downloaded from http://cpc.cs.qub.ac.uk/summaries/AEES_v1_0.html.
Alternatively, the code can be obtained from https://bitbucket.org/gkw/gkw;
versions found there may contain new features (and also new bugs).

Although you are under no obligation to do so, it is most helpful if
you report any bugs you find in the code (and perhaps submit any fixes
if you have them). Please do this via the issue tracker at
http://bitbucket.org/gkw/gkw/issues.
The issue tracker is also used for more general discussions concerning the code.

Participation in the project is encouraged and contributions are
welcome; please contact a project owner at http://bitbucket.org/gkw/
if you are interested.


 CONTENTS SUMMARY
==================

**README.md**  
(this file)

**LICENSE**  
license details

**GNUmakefile**  
top level makefile which builds the code in src/

**src/**  
    All Fortran source code files.

**config/**  
    Host specific makefiles, see the GNUmakefile for details.

**doc/**  
    Contains the code manual and other pieces of documentation.
    
**doc/input/**  
    Some sample input files for the code.
    
**doc/benchmarks/**  
    A directory containing a collection of code inputs for specific
    benchmarks; see any README files there for detail.

**scripts/**  
    Various scripts, mainly connected with running and testing the code.

**libs/ **  
    Third party libraries distributed with GKW, used to extend the code 
    functionality; these must be built separately if required.

**matlab/**  
    A collection of matlab files used by some developers that may
    be of use to others who have access to the matlab program.

**tests/**  
    Short test input and reference files for use when making changes to
    the code.


 LICENSE
=========

(see the LICENSE file)

This file is part of GKW.

GKW is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GKW is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GKW.  If not, see <http://www.gnu.org/licenses/>.