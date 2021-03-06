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
    Foam::flameletModel::coherentFlameModel2b

Description
    Implementation of the flamelet model proposed in
    \verbatim
        Duclos, J. M., Veynante, D., and Poinsot, T.,
        "A comparison of flamelet models for premixed turbulent combustion",
	Combustion and Flame, 1993,
        95:pp. 101 - 117.
    \endverbatim

    The name of the model was adopted from
    \verbatim
        Poinsot, T., and Veynante, D.,
        "Theoretical and numerical combustion", Third Edition, 2012,
	Table 5.3 (pp 241).
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        coherentFlameModel2bCoeffs
        {
            alphaSigma    2.1;
            betaSigma     1.0;
            CSigma        0.5;
        }
    \endverbatim

SourceFiles
    coherentFlameModel2b.C

\*---------------------------------------------------------------------------*/

#ifndef coherentFlameModel2b_H
#define coherentFlameModel2b_H

#include "flameletModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletModelA Declaration
\*---------------------------------------------------------------------------*/

class coherentFlameModel2b
:
    public flameletModel
{
private:

    //-Dictionary
    dictionary dictionary_;

    //Constants
    scalar alphaSigma_;
    scalar betaSigma_;
    scalar CSigma_;
    scalar Clt_;
    scalar quenchingCoeff_;
    bool fittedGammaK_;

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
        coherentFlameModel2b(const coherentFlameModel2b&);

        //- Disallow default bitwise assignment
        void operator=(const coherentFlameModel2b&);

public:

    //- Runtime type information
    TypeName("coherentFlameModel2b");

    // Constructors

        //- Construct for the given flameletModel
        coherentFlameModel2b
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
    virtual ~coherentFlameModel2b();

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
