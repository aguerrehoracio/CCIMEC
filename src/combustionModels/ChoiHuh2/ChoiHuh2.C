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

\*---------------------------------------------------------------------------*/

#include "ChoiHuh2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ChoiHuh2, 0);

    addToRunTimeSelectionTable
    (
        flameletModel,
        ChoiHuh2,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChoiHuh2::ChoiHuh2
(   
	volScalarField& rho,
	volVectorField& U,
	volScalarField& Su,
	volScalarField& Sigma,
	volScalarField& b,
	psiuReactionThermo& thermo,
	compressible::turbulenceModel& turbulence,
	IOdictionary& mdData
)
:
    flameletModel
    (
        rho
    ),
    dictionary_(mdData.subDict("ChoiHuh2Coeffs")),
    alphaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "alphaSigma",
	    10.0
	)
    ),
    betaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "betaSigma",
	    0.1
	)
    ),
    ltc_
    (
        dictionary_.lookupOrDefault
	(
	    "l_tc",
	    0.00126
	)
    ),
    rho_(rho),	
    U_(U),
    Su_(Su),
    Sigma_(Sigma),
    b_(b),
    thermo_(thermo),
    turbulence_(turbulence)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ChoiHuh2::~ChoiHuh2()
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::ChoiHuh2::update()
{
  Info<<"Updating Sigma Source terms"<<endl;

  ProdRateForSigma_ = rho_ * alphaSigma_ * sqrt(2.0 / 3.0 * turbulence_.k()) / 
    dimensionedScalar("ltc", dimLength, ltc_ + SMALL);

  DestrRateForSigma_ = rho_ * betaSigma_ * Su_ * 
    Sigma_ / (b_ * (1.0 - b_) + SMALL);

}

void Foam::ChoiHuh2::name
()
{
  Info << endl << "I am CH2 model" << endl;
}


// ************************************************************************* //
