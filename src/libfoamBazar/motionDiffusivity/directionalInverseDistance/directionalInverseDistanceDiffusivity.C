/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "directionalInverseDistanceDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "patchWave.H"
#include "HashSet.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directionalInverseDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        directionalInverseDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directionalInverseDistanceDiffusivity::directionalInverseDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    uniformDiffusivity(mesh, mdData),
    dirName_(mdData),
    patchNames_(mdData)
{
    if (dirName_ == "x")
    {
        cmpt_ = 0;
    }
    else if (dirName_ == "y")
    {
        cmpt_ = 1;
    }
    else if (dirName_ == "z")
    {
        cmpt_ = 2;
    }
    else
    {
        FatalErrorInFunction() << "expect x, y or z ... found: " << dirName_ << nl << abort(FatalError);
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directionalInverseDistanceDiffusivity::~directionalInverseDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::directionalInverseDistanceDiffusivity::y() const
{
    labelHashSet patchSet(mesh().boundaryMesh().patchSet(patchNames_));

    if (patchSet.size())
    {
        return tmp<scalarField>
        (
            new scalarField(patchWave(mesh(), patchSet, false).distance())
        );
    }
    else
    {
        return tmp<scalarField>(new scalarField(mesh().nCells(), 1.0));
    }
}


void Foam::directionalInverseDistanceDiffusivity::correct()
{
    volScalarField y_
    (
        IOobject
        (
            "y",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );
    y_.primitiveFieldRef() = y();
    scalar minVal = min(y_.primitiveField());

    scalar bbmin = VGREAT;
    scalar bbmax = -VGREAT;

    labelHashSet patchSet = mesh().boundaryMesh().patchSet(patchNames_);
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        boundBox bb(mesh().boundaryMesh()[iter.key()].localPoints());
        bbmin = min(bb.min()[cmpt_], bbmin);
        bbmax = max(bb.max()[cmpt_], bbmax);
    }

    reduce(bbmin, minOp<scalar>());
    reduce(bbmax, maxOp<scalar>());

    const volVectorField& C = mesh().C();
    forAll(C, i)
    {
        y_[i] = max(minVal, min(fabs(C[i][cmpt_] - bbmin), fabs(C[i][cmpt_] - bbmax)));
    }
    y_.correctBoundaryConditions();

    faceDiffusivity_ = 1.0/fvc::interpolate(y_);
}


// ************************************************************************* //
