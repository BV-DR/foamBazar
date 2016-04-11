/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::sixDoFRigidBodyMotionSolver

Description
	Derived class of sixDoFRigidBodyMotionSolver and allows to register as an object in object registry

SourceFiles
    sixDoFRigidBodyMotionSolverBis.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotionSolverBis_H
#define sixDoFRigidBodyMotionSolverBis_H

#include "sixDoFRigidBodyMotionSolver.H"
#include "access6DoFRigidBodyMotions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class sixDoFRigidBodyMotionSolverBis Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotionSolverBis
:
    public sixDoFRigidBodyMotionSolver
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDoFRigidBodyMotionSolverBis
        (
            const sixDoFRigidBodyMotionSolverBis&
        );

        //- Disallow default bitwise assignment
        void operator=(const sixDoFRigidBodyMotionSolverBis&);


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyMotionBis");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        sixDoFRigidBodyMotionSolverBis
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~sixDoFRigidBodyMotionSolverBis();

    //- Pure virtual writaData function.
    //  Must be defined in derived types of regIOobject
    //virtual bool writeData(Ostream& os) const
    //{
    //    return os.good();
    //}


    // Member Functions

    const sixDoFRigidBodyMotion& motion() const { return motion_; }	

/*
        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Write state using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;

        //- Read dynamicMeshDict dictionary
        virtual bool read();

*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
