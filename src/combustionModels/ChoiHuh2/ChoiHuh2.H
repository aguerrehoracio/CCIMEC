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
    Foam::flameletModel::ChoiHuh2

Description
    Implementation of the flamelet model proposed by Choi and Huh 
    described in
    \verbatim
        Choi, Ch. R., and Huh, K. Y.,
        "Development of a coherent flamelet model for a spark-ignited
	turbulent premixed flame in a closed vessel",
	Combustion and Flame, 1998,
        114:pp. 336 - 348.
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        ChoiHuh2Coeffs
        {
            alphaSigma    10.0;
            betaSigma     0.1;
	    ltc           0.00126;
        }
    \endverbatim

SourceFiles
    ChoiHuh2.C

\*---------------------------------------------------------------------------*/

#ifndef ChoiHuh2_H
#define ChoiHuh2_H

#include "flameletModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletModelA Declaration
\*---------------------------------------------------------------------------*/

class ChoiHuh2
:
    public flameletModel
{
private:

    //-Dictionary
    dictionary dictionary_;

    //Constants
    scalar alphaSigma_;
    scalar betaSigma_;
    scalar ltc_;

    //Fields
    volScalarField& rho_;    
    volVectorField& U_;
    volScalarField& Su_;
    volScalarField& Sigma_;
    volScalarField& b_;

    //Thermo
    psiuReactionThermo& thermo_;

    //TurbulenceModels
    compressible::turbulenceModel& turbulence_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        ChoiHuh2(const ChoiHuh2&);

        //- Disallow default bitwise assignment
        void operator=(const ChoiHuh2&);

public:

    //- Runtime type information
    TypeName("ChoiHuh2");

    // Constructors

        //- Construct for the given flameletModel
        ChoiHuh2
        (
            volScalarField& rho,
            volVectorField& U,
	    volScalarField& Su,
	    volScalarField& Sigma,
	    volScalarField& b,
	    psiuReactionThermo& thermo,
	    compressible::turbulenceModel& turbulence,
	    IOdictionary& mdData
        );

    //- Destructor
    virtual ~ChoiHuh2();

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
