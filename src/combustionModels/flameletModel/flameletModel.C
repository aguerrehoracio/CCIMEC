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

#include "flameletModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(flameletModel, 0);
    defineRunTimeSelectionTable(flameletModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletModel::flameletModel
(
	volScalarField& rho
)
:
	ProdRateForSigma_
    (
        IOobject
        (
            "ProdRateForSigma",
            rho.mesh().time().timeName(),
            rho.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho.mesh(),
        dimensionSet(1, -3, -1, 0, 0)
    ),

    DestrRateForSigma_
    (
        IOobject
        (
            "DestrRateForSigma",
            rho.mesh().time().timeName(),
            rho.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho.mesh(),
        dimensionSet(1, -3, -1, 0, 0)
    ),

    I0_
    (
        IOobject
        (
            "I0",
            rho.mesh().time().timeName(),
            rho.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho.mesh(),
        dimensionSet(0, 0, 0, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::flameletModel> Foam::flameletModel::New
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
{
    const word flameletType(mdData.lookup("flameletModel"));

    Info<< "Selecting flamelet diffusion: " << flameletType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(flameletType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "flameletModel::New("
            "const IOdictionary& dict)"
        )   << "Unknown diffusion type "
            << flameletType << nl << nl
            << "Valid diffusion types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<flameletModel>
      (cstrIter()(rho, U, Su, Sigma, b, thermo, turbulence, mdData));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameletModel::~flameletModel()
{}



// ************************************************************************* //
