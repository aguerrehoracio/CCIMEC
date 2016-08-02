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

#include "coherentFlameModel2a.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coherentFlameModel2a, 0);

    addToRunTimeSelectionTable
    (
        flameletModel,
        coherentFlameModel2a,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coherentFlameModel2a::coherentFlameModel2a
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
    dictionary_(mdData.subDict("coherentFlameModel2aCoeffs")),
    alphaSigma_
    (
        dictionary_.lookupOrDefault
	(
	    "alphaSigma",
	    2.1
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
    Clt_
    (
        dictionary_.lookupOrDefault
	(
	    "Coeff_lt",
	    1.0
	 )
    ),
    quenchingCoeff_
    (
        dictionary_.lookupOrDefault
	(
	    "quenchingCoeff",
	    1.0
	)
    ),
    fittedGammaK_
    (
        dictionary_.lookupOrDefault
	(
	    "fittedGammaK",
	    false
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

Foam::coherentFlameModel2a::~coherentFlameModel2a()
{}

// * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //
void Foam::coherentFlameModel2a::update()
{
  Info<<"Updating Sigma Source terms"<<endl;

  volVectorField M(fvc::grad(b_));
  volScalarField mgb_ = mag(M);

  dimensionedScalar dSigma = 1.0e-3*
    (b_* (scalar(1.0) - b_) * mgb_)().weightedAverage(rho_.mesh().V())
    /((b_ * (scalar(1.0) - b_))().weightedAverage(rho_.mesh().V()) + SMALL)
    + dimensionedScalar("dSig", Sigma_.dimensions(), SMALL);

  M /= (max(Sigma_, mgb_) + dSigma);

  // volScalarField magM = mag(M);
  // M /= Sigma_ + dimensionedScalar("tol", dimless/dimLength, 1e-2);
  // M /= (magM + dimensionedScalar("tol", dimless/dimLength, SMALL));
  
  volScalarField orientationFactor = scalar(1.0) - (M & M);
  orientationFactor.max(0.0);
  orientationFactor.min(1.0);

  volTensorField A_ = I * (1.0 - orientationFactor / 3.0) - (M * M);
  volTensorField gradU(fvc::grad(U_));

  volScalarField up(sqrt((2.0/3.0)*turbulence_.k()));
  
  //Efficiency function from the ITNFS model
  volScalarField GammaK("GammaK", rho_/rho_);

  scalar T1 = min(thermo_.Tu()).value();
  scalar T2 = max(thermo_.Tb()).value();

  volScalarField deltaL(2.0 * thermo_.alpha() / thermo_.rhou() / 
			Su_ * pow(T2/T1, 0.7));

  volScalarField lt(Clt_ * pow(up, 3) / (turbulence_.epsilon() + 
					 dimensionedScalar("tol", 
                 		         pow(dimVelocity,2)/dimTime, SMALL)));

  volScalarField lRatio(lt / deltaL);

  volScalarField uRatio(up / Su_);

  volScalarField s = max
    (
       log10 (lRatio),
       scalar(-0.4+SMALL)
    );

  if(fittedGammaK_){
    GammaK = 0.75 * exp(-1.2 / pow(uRatio, 0.3)) * pow(lRatio, 2.0/3.0);
  } else {
    volScalarField sigma1 = 2.0 / 3.0 * 
      (1.0 - 0.5 * exp(-pow(uRatio, 1.0/3.0)));
    volScalarField r = -exp(-(s+0.4)) / (s+0.4) + 
      (1.0 - exp(-(s+0.4))) * (sigma1 * s - 0.11);

    GammaK = pow(10.0, r);
  }

  if(quenchingCoeff_ != 0.0){
    volScalarField g = (0.7 + 1.0 / s) * exp(-s) + 
      (1.0 - exp(-s)) * (1.0 + 0.36 * s);
    volScalarField x = (log10 (uRatio) - g) / s / 0.04;
    volScalarField Pq = 0.5 * (1.0 + tanh (sign(x) * x * x));
    
    GammaK -= quenchingCoeff_ * 1.5 * lRatio / uRatio * log (1.0 / (1.0 - Pq));
  }

  volScalarField P1 = rho_ * alphaSigma_ * GammaK * turbulence_.epsilon() / 
    (turbulence_.k() + dimensionedScalar("tol", pow(dimVelocity,2), SMALL));
  
  volScalarField P2 = rho_ * (A_ && gradU);

  ProdRateForSigma_ = P1 + P2;

  DestrRateForSigma_ = rho_ * betaSigma_ * 
    (Su_ + CSigma_ * sqrt(3.0/2.0) * up) * Sigma_/(b_ + SMALL);

}

void Foam::coherentFlameModel2a::name
()
{
  Info << endl << "I am CFM2a model" << endl;
}


// ************************************************************************* //
