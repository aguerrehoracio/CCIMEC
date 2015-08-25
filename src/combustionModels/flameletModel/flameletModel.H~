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
    Foam::flameletModel

Description
    Abstract base class for cell-centre mesh motion diffusivity.

SourceFiles
    flameletModel.C

\*---------------------------------------------------------------------------*/

#ifndef flameletModel_H
#define flameletModel_H

#include "surfaceFieldsFwd.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "IOdictionary.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletModel Declaration
\*---------------------------------------------------------------------------*/

class flameletModel
{
protected:

    // Private data
	volScalarField ProdRateForSigma_;
	volScalarField DestrRateForSigma_;
	volScalarField I0_;

public:

    //- Runtime type information
    TypeName("flameletModel");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            flameletModel,
            dictionary,
            (
				volScalarField& rho,	
				volVectorField& U,
				volScalarField& Su,
				volScalarField& Sigma,
				volScalarField& b,
				compressible::turbulenceModel& turbulence,
				IOdictionary& mdData
            ),
            (rho, U, Su, Sigma, b, turbulence, mdData)
        );


    // Selectors

        //- Select null constructed
        static autoPtr<flameletModel> New
        (
			volScalarField& rho,
            volVectorField& U,
			volScalarField& Su,
			volScalarField& Sigma,
			volScalarField& b,
			compressible::turbulenceModel& turbulence,
			IOdictionary& mdData
        );


    // Constructors

        //- Construct for the given flameletModel
        flameletModel
        (
            volScalarField& rho
        );


    //- Destructor
    virtual ~flameletModel();


    // Member Functions

     virtual void name()=0;

     //Functions necessary on model 
	 virtual void update()=0;

	 const volScalarField& SigmaProd() const
	 {
	 	return ProdRateForSigma_;
	 }

	 const volScalarField& SigmaDestr() const
	 {
	 	return DestrRateForSigma_;
	 }

	 const volScalarField& I0() const
	 {
	 	return I0_;
	 }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
