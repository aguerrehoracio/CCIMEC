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

#include "flameletModelA.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flameletModelA, 0);

    addToRunTimeSelectionTable
    (
        flameletModel,
        flameletModelA,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletModelA::flameletModelA
(   
	volScalarField& rho,
    volVectorField& U,
	volScalarField& Su,
	volScalarField& Sigma,
	volScalarField& b,
	compressible::turbulenceModel& turbulence,
	IOdictionary& mdData
)
:
    flameletModel
    (
        rho
    ),
    dictionary_(mdData.subDict("flameletModelACoeffs")),
    alphaSigma_(readScalar(dictionary_.lookup("alphaSigma"))),
    GammaK_(readScalar(dictionary_.lookup("GammaK"))),
    betaSigma_(readScalar(dictionary_.lookup("betaSigma"))),
	rho_(rho),	
	U_(U),
	Su_(Su),
	Sigma_(Sigma),
	b_(b),
	turbulence_(turbulence)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameletModelA::~flameletModelA()
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::flameletModelA::update()
{
	Info<<"Updating Sigma Source terms"<<endl;

    volScalarField P1 = rho_ * alphaSigma_ * GammaK_ * turbulence_.epsilon()/turbulence_.k();

    volScalarField P2 = rho_ * (2.0 / 3.0) * fvc::div(U_);

    ProdRateForSigma_ = P1 + P2;
    DestrRateForSigma_ = rho_ * rho_ * betaSigma_ * Su_ * (Sigma_)/(b_ + SMALL);

}

void Foam::flameletModelA::name
()
{
    scalar forInfo = readScalar(dictionary_.lookup("GammaK"));

	Info << endl << "I am model A with a GammaK equal to " << forInfo << endl;
}


// ************************************************************************* //
