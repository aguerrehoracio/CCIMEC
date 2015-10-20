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

#include "ChoiHuh1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ChoiHuh1, 0);

    addToRunTimeSelectionTable
    (
        flameletModel,
        ChoiHuh1,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChoiHuh1::ChoiHuh1
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
    dictionary_(mdData.subDict("ChoiHuh1Coeffs")),
    alphaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "alphaSigma",
	    1.2
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
    rho_(rho),	
    U_(U),
    Su_(Su),
    Sigma_(Sigma),
    b_(b),
    thermo_(thermo),
    turbulence_(turbulence)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ChoiHuh1::~ChoiHuh1()
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::ChoiHuh1::update()
{
  Info<<"Updating Sigma Source terms"<<endl;

  volScalarField nu = (thermo_.muu() / thermo_.rhou() * b_ +
  		       thermo_.mub() / thermo_.rhob() * (1.0 - b_));

  ProdRateForSigma_ = rho_ * alphaSigma_ * 
    sqrt(turbulence_.epsilon() / nu / 15.0);

  DestrRateForSigma_ = rho_ * betaSigma_ * Su_ * 
    Sigma_ / (b_ * (1.0 - b_) + SMALL);

}

void Foam::ChoiHuh1::name
()
{
  Info << endl << "I am CH1 model" << endl;
}


// ************************************************************************* //
