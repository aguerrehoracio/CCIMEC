/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::flameletModelA

Description
    

SourceFiles
    flameletModelA.C

\*---------------------------------------------------------------------------*/

#ifndef flameletModelA_H
#define flameletModelA_H

#include "flameletModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletModelA Declaration
\*---------------------------------------------------------------------------*/

class flameletModelA
:
    public flameletModel
{
private:

    //-Dictionary
	dictionary dictionary_;

    //Constants
	scalar alphaSigma_;
    scalar GammaK_;
    scalar betaSigma_;

    //Fields
    volScalarField& rho_;    
	volVectorField& U_;
	volScalarField& Su_;
	volScalarField& Sigma_;
	volScalarField& b_;

    //TurbulenceModels
	compressible::turbulenceModel& turbulence_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        flameletModelA(const flameletModelA&);

        //- Disallow default bitwise assignment
        void operator=(const flameletModelA&);


public:

    //- Runtime type information
    TypeName("flameletModelA");


    // Constructors

        //- Construct for the given flameletModel
        flameletModelA
        (
            volScalarField& rho,
            volVectorField& U,
			volScalarField& Su,
			volScalarField& Sigma,
			volScalarField& b,
			compressible::turbulenceModel& turbulence,
			IOdictionary& mdData
        );


    //- Destructor
    virtual ~flameletModelA();


    // Member Functions
    
     //Functions necessary on model 
	 virtual void update();

    virtual void name();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
