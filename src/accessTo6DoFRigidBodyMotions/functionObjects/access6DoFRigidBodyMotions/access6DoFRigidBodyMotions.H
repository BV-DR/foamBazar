/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::access6DoFRigidBodyMotions

Description
    This class provides access to motion displacement in sixDofRigidMotionSolver 
    through functionObject.

Author
    Charles Monroy, Bureau Veritas.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef access6DoFRigidBodyMotions_H
#define access6DoFRigidBodyMotions_H

#include "functionObjectFile.H"
#include "dictionary.H"
#include "Time.H"
#include "fvMesh.H"
#include "objectRegistry.H"
#include "sixDoFRigidBodyMotionSolver.H"
#include "sixDoFRigidBodyMotion.H"
#include "sixDoFRigidBodyMotionState.H"
#include "sixDoFRigidBodyMotionSolverBis.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sixDoFRigidBodyMotionSolverBis;

//class Time;

/*---------------------------------------------------------------------------*\
                         Class access6DoFRigidBodyMotions Declaration
\*---------------------------------------------------------------------------*/

class access6DoFRigidBodyMotions
:
    public functionObjectFile
{
    // Private data
        //- Name
        const word name_;

        //- Reference to main object registry
        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- for writeFileHeader() only
        label writePrecision_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        access6DoFRigidBodyMotions(const access6DoFRigidBodyMotions&);

        //- Disallow default bitwise assignment
        void operator=(const access6DoFRigidBodyMotions&);


        //- Output file header information
        virtual void writeFileHeader();

public:

    //- Runtime type information
    TypeName("access6DoFRigidBodyMotions");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allows the possibility to load fields from files
        access6DoFRigidBodyMotions
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles=false
        );


    //- Destructor
    //virtual ~access6DoFRigidBodyMotions();


    // Member Functions

        //- Return name of the function object
        virtual const word& name() const
        {
            return name_;
        }

        //- Called at the start of the time-loop
        virtual void start()
        {}

        //- Read and set the function object if its data has changed
        virtual void read(const dictionary& dict)
        {}

        //- execute is called at each ++ or += of the time-loop
        virtual void execute()
        {}

        //- 
        virtual void end()
        {}

        //- 
        virtual void timeSet()
        {}

        //- 
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}
        
        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
