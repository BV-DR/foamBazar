/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    flattenMesh

Description
    Flattens the front and back planes of a 2D cartesian mesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "twoDPointCorrector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("project point(s) on patch(es) on to a plane");
    argList::addOption
    (
        "patches",
        "(patch0 .. patchN)",
        "project point(s) on patch(es) on to a plane"
    );
    argList::addOption
    (
        "planeRef", "(x y z)", "ref. point of the plane"
    );
    argList::addOption
    (
        "planeNormal", "(x y z)", "normal vector of the plane"
    );
   
    #include "setRootCase.H"
    
    if (!args.optionFound("patches")) FatalErrorIn("") << "option -patches not defined" << abort(FatalError);
    if (!args.optionFound("planeRef")) FatalErrorIn("") << "option -planeRef not defined" << abort(FatalError);
    if (!args.optionFound("planeNormal")) FatalErrorIn("") << "option -planeNormal not defined" << abort(FatalError);
    
    #include "createTime.H"
    #include "createPolyMesh.H"   

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    vector planeRef = vector(args.optionLookup("planeRef")());
    vector planeNormal = vector(args.optionLookup("planeNormal")());
    labelHashSet includePatches = bMesh.patchSet(wordReList(args.optionLookup("patches")()));
    
    Info << "Project mesh point(s) onto plane (ref, normal): " << planeRef << " " << planeNormal << endl;

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        const labelList& meshPoints = pp.meshPoints();
        label nPoints(pp.nPoints());
        Info << "patch: " << pp.name() << " (" << returnReduce(nPoints,sumOp<label>()) << ")" << endl;
        forAll(meshPoints, i)
        {
            point& p(points[meshPoints[i]]);
            p -= ((p - planeRef) & planeNormal) * planeNormal;
        }
    }

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));
    Info<< "Writing new points to file: "
    << points.instance()/points.local()/points.filePath().name() << nl << endl;
    points.write();

    Info<< "End\n" << endl;
    
    return 0;
}

// ************************************************************************* //



