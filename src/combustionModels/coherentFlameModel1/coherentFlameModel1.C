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

#include "coherentFlameModel1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coherentFlameModel1, 0);

    addToRunTimeSelectionTable
    (
        flameletModel,
        coherentFlameModel1,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coherentFlameModel1::coherentFlameModel1
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
    dictionary_(mdData.subDict("coherentFlameModel1Coeffs")),
    alphaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "alphaSigma",
	    1.7
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
    CSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "CSigma",
	    0.5
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

Foam::coherentFlameModel1::~coherentFlameModel1()
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::coherentFlameModel1::update()
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

  volTensorField A_ = I * (scalar(1.0) - orientationFactor / 3.0) - (M * M);
  volTensorField gradU(fvc::grad(U_));

  volScalarField P1 = rho_ * alphaSigma_ * turbulence_.epsilon() / 
    (turbulence_.k() + dimensionedScalar("ktol", pow(dimVelocity,2), 1.0e-3));
  
  volScalarField P2 = rho_ * (A_ && gradU);

  ProdRateForSigma_ = P1 + P2;

  DestrRateForSigma_ = rho_ * betaSigma_ * 
    (Su_ + CSigma_ * sqrt(turbulence_.k())) * Sigma_/(b_ + SMALL);

}

void Foam::coherentFlameModel1::name
()
{
  Info << endl << "I am CFM1 model" << endl;
}


// ************************************************************************* //
