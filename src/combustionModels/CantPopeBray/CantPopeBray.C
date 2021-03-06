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

#include "CantPopeBray.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CantPopeBray, 0);

    addToRunTimeSelectionTable
    (
        flameletModel,
        CantPopeBray,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CantPopeBray::CantPopeBray
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
    dictionary_(mdData.subDict("CantPopeBrayCoeffs")),
    alphaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "alphaSigma",
	    0.28
	)
    ),
    betaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "betaSigma",
	    1.0
	)
    ),
    aCoeff_
    (
        dictionary_.lookupOrDefault
	(
	    "aCoeff",
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

Foam::CantPopeBray::~CantPopeBray()
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::CantPopeBray::update()
{
  Info<<"Updating Sigma Source terms"<<endl;

  volVectorField M(fvc::grad(b_));
  volScalarField mgb_ = mag(M);

  dimensionedScalar dSigma = 1.0e-3*
    (b_* (scalar(1.0) - b_) * mgb_)().weightedAverage(rho_.mesh().V())
    /((b_ * (scalar(1.0) - b_))().weightedAverage(rho_.mesh().V()) + SMALL)
    + dimensionedScalar("dSig", Sigma_.dimensions(), SMALL);

  M /= (max(Sigma_, mgb_) + dSigma);

  volScalarField magM = mag(M);

  // M /= Sigma_ + dimensionedScalar("tol", dimless/dimLength, 1e-2);
  // M /= (magM + dimensionedScalar("tol", dimless/dimLength, SMALL));

  volScalarField orientationFactor = scalar(1.0) - (M & M);
  orientationFactor.max(0.0);
  orientationFactor.min(1.0);
  Info<< "min(alpha) = " << min(orientationFactor).value() << endl;

  volTensorField A_ = I * (1.0 - orientationFactor / 3.0) - (M * M);
  volTensorField gradU(fvc::grad(U_));

  volScalarField R_ = b_ * turbulence_.epsilon() /
    (Sigma_ * Su_ * turbulence_.k() + 
     dimensionedScalar("tol", sqr(dimVelocity)/dimTime, SMALL));

  volScalarField nu = (thermo_.muu() / thermo_.rhou() * b_ +
		       thermo_.mub() / thermo_.rhob() * (1.0 - b_));

  volScalarField P1 = rho_ * alphaSigma_ * orientationFactor *
    sqrt(turbulence_.epsilon() / nu);
  
  volScalarField P2 = rho_ * (A_ && gradU);

  ProdRateForSigma_ = P1 + P2;

  DestrRateForSigma_ = rho_ * betaSigma_ * orientationFactor / 3.0 *
    (2.0 + exp(-aCoeff_*R_)) * Su_ * Sigma_ / (b_ + SMALL);
}

void Foam::CantPopeBray::name
()
{
  Info << endl << "I am CPB model" << endl;
}


// ************************************************************************* //
