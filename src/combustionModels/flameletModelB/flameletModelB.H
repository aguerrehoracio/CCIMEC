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
    Foam::flameletModelB

Description


SourceFiles
    flameletModelB.C

\*---------------------------------------------------------------------------*/

#ifndef flameletModelB_H
#define flameletModelB_H

#include "flameletModel.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flameletModelB Declaration
\*---------------------------------------------------------------------------*/

class flameletModelB
:
    public flameletModel
{

protected:

    // Protected data

        surfaceScalarField faceDiffusivity_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        flameletModelB(const flameletModelB&);

        //- Disallow default bitwise assignment
        void operator=(const flameletModelB&);


public:

    //- Runtime type information
    TypeName("flameletModelB");


    // Constructors

        //- Construct for the given flamelet model
        flameletModelB(Istream& mdData);


    //- Destructor
    virtual ~flameletModelB();


    // Member Functions

    virtual void name();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
