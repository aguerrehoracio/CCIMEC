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
    Foam::flameletModel::ChoiHuh1

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
        ChoiHuh1Coeffs
        {
            alphaSigma    1.2;
            betaSigma     0.1;
        }
    \endverbatim

SourceFiles
    ChoiHuh1.C

\*---------------------------------------------------------------------------*/

#ifndef ChoiHuh1_H
#define ChoiHuh1_H

#include "flameletModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletModelA Declaration
\*---------------------------------------------------------------------------*/

class ChoiHuh1
:
    public flameletModel
{
private:

    //-Dictionary
    dictionary dictionary_;

    //Constants
    scalar alphaSigma_;
    scalar betaSigma_;

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
        ChoiHuh1(const ChoiHuh1&);

        //- Disallow default bitwise assignment
        void operator=(const ChoiHuh1&);

public:

    //- Runtime type information
    TypeName("ChoiHuh1");

    // Constructors

        //- Construct for the given flameletModel
        ChoiHuh1
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
    virtual ~ChoiHuh1();

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
