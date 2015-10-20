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

#include "ChengDiringer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ChengDiringer, 0);

    addToRunTimeSelectionTable
    (
        flameletModel,
        ChengDiringer,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChengDiringer::ChengDiringer
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
    dictionary_(mdData.subDict("ChengDiringerCoeffs")),
    alphaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "alphaSigma",
	    8.5
	)
    ),
    betaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "betaSigma",
	    10.0
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

Foam::ChengDiringer::~ChengDiringer()
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::ChengDiringer::update()
{
  Info<<"Updating Sigma Source terms"<<endl;

  //By the time, no quenching due to excesive strain is incorporated
  //(see the paper by Cheng & Diringer)

  ProdRateForSigma_ = rho_ * alphaSigma_ * turbulence_.epsilon() / 
    (turbulence_.k() + dimensionedScalar("tol", pow(dimVelocity,2), SMALL));

  DestrRateForSigma_ = rho_ * betaSigma_ * Su_ * Sigma_/(b_ + SMALL);

}

void Foam::ChengDiringer::name
()
{
  Info << endl << "I am CD model" << endl;
}


// ************************************************************************* //
