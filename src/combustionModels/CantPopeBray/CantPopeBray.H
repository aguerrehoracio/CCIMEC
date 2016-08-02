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
    Foam::flameletModel::CantPopeBray

Description
    Implementation of the flamelet model proposed by Cant, Pope and
    Bray described in
    \verbatim
        Cant, R. S., Pope, S. B., and Bray, K. N. C.,
        "Modelling of flamelet surface-to-volume ratio in turbulent 
	premixed combustion",
        Twenty-Third Symposium (International) on Combustion/The
	Combustion Institute, 1990,
        pp. 809 - 815.
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        CantPopeBrayCoeffs
        {
            alphaSigma    0.28;
            betaSigma     1.0;
            aCoeff        10.0;
        }
    \endverbatim

    alphaSigma = C_A coefficient from the paper.
    aCoeff = A coefficient from the paper.

SourceFiles
    CantPopeBray.C

\*---------------------------------------------------------------------------*/

#ifndef CantPopeBray_H
#define CantPopeBray_H

#include "flameletModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletModelA Declaration
\*---------------------------------------------------------------------------*/

class CantPopeBray
:
    public flameletModel
{
private:

    //-Dictionary
    dictionary dictionary_;

    //Constants
    scalar alphaSigma_;
    scalar betaSigma_;
    scalar aCoeff_;

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
        CantPopeBray(const CantPopeBray&);

        //- Disallow default bitwise assignment
        void operator=(const CantPopeBray&);

public:

    //- Runtime type information
    TypeName("CantPopeBray");

    // Constructors

        //- Construct for the given flameletModel
        CantPopeBray
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
    virtual ~CantPopeBray();

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
