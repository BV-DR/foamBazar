/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
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
    accessFluidForces

Author
    Charles Monroy, Bureau Veritas.

\*---------------------------------------------------------------------------*/

#include "access6DoFRigidBodyMotions.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(access6DoFRigidBodyMotions, 0);
};
    


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::access6DoFRigidBodyMotions::access6DoFRigidBodyMotions
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name),
    name_(name),
    obr_(obr),
    active_(true),
    writePrecision_(dict.lookupOrDefault<label>("writePrecision", 6))
{
    if (!active_) return;


    DynamicList<word> names(1);
    const word access6DoFRigidBodyMotionsType(dict.lookup("type"));
    names.append(access6DoFRigidBodyMotionsType);
    functionObjectFile::resetNames(names);

     writeFileHeader();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::access6DoFRigidBodyMotions::write()
{
    if (!active_) return;

    functionObjectFile::write();

    if (Pstream::master())
    {
        const scalar t(obr_.time().value());
         
	const sixDoFRigidBodyMotionSolverBis& aaaa=
	obr_.lookupObject<sixDoFRigidBodyMotionSolverBis>("dynamicMeshDict");


        const sixDoFRigidBodyMotionState myMotionState(aaaa.motion().state());
 
	const point centreOfRotation(myMotionState.centreOfRotation());
	const vector v(myMotionState.v());  
 	const vector a(myMotionState.a());      

        file() << t;
	file() << " " << centreOfRotation.x();
	file() << " " << centreOfRotation.y();
	file() << " " << centreOfRotation.z();
	file() << " " << "0.0";   // to add... rotation...
	file() << " " << "0.0";
	file() << " " << "0.0";

	file() << " " << v.x();
	file() << " " << v.y();
	file() << " " << v.z();
	file() << " " << "0.0";   // to add...
	file() << " " << "0.0";
	file() << " " << "0.0";

	file() << " " << a.x();
	file() << " " << a.y();
	file() << " " << a.z();
	file() << " " << "0.0";   // to add...
	file() << " " << "0.0";
	file() << " " << "0.0";

        file() << endl;
/*
            inline const point& centreOfRotation() const;


        //- Give access to sixDOF object
        inline const fsSixDOFode& sixDOF() const
        {
            return sixDOF_;
        } */
/*
        const floatingBodies& bodies(fltBodies());
        forAll(bodies, i)
        {
            const fsSixDOFode& sixdof(bodies[i].sixDOF());
            const vector Xrel(sixdof.Xrel());
            const vector dotXrel(sixdof.dotXrel());
            const scalarField th(sixdof.theta());
            const vector omg(sixdof.omega());
            const vector rot(sixdof.rot().eulZYX(th) * 180.0/M_PI);
            const vector linAcc(sixdof.linAcc());
            const vector angAcc(sixdof.angAcc());
            
            file(i) << t;
            if (outputDofs_[0]) { file(i) << " " << Xrel.x(); }
            if (outputDofs_[1]) { file(i) << " " << Xrel.y(); }
            if (outputDofs_[2]) { file(i) << " " << Xrel.z(); }
            if (outputDofs_[3]) { file(i) << " " << rot.x(); }
            if (outputDofs_[4]) { file(i) << " " << rot.y(); }
            if (outputDofs_[5]) { file(i) << " " << rot.z(); }
            //
            if (outputDofs_[0]) { file(i) << " " << dotXrel.x(); }
            if (outputDofs_[1]) { file(i) << " " << dotXrel.y(); }
            if (outputDofs_[2]) { file(i) << " " << dotXrel.z(); }
            if (outputDofs_[3]) { file(i) << " " << omg.x(); }
            if (outputDofs_[4]) { file(i) << " " << omg.y(); }
            if (outputDofs_[5]) { file(i) << " " << omg.z(); }
            //
            if (outputDofs_[0]) { file(i) << " " << linAcc.x(); }
            if (outputDofs_[1]) { file(i) << " " << linAcc.y(); }
            if (outputDofs_[2]) { file(i) << " " << linAcc.z(); }
            if (outputDofs_[3]) { file(i) << " " << angAcc.x(); }
            if (outputDofs_[4]) { file(i) << " " << angAcc.y(); }
            if (outputDofs_[5]) { file(i) << " " << angAcc.z(); }
            //
            if (outputDofs_[6] && (sixdof.nModes() !=0 ))
            {
                const scalarField qf(sixdof.modalAmpl());
                const scalarField dqf(sixdof.modalVel());
                const scalarField ddqf(sixdof.modalAcc());
                forAll(qf, k) file(i) << " " << qf[k];
                forAll(dqf, k) file(i) << " " << dqf[k];
                forAll(ddqf, k) file(i) << " " << ddqf[k];
            }
            file(i) << endl;
            
            //file(i) << t
            //<< " " << Xrel.x() << " " << Xrel.y() << " " << Xrel.z()
            //<< " " << rot.x() << " " << rot.y() << " " << rot.z()
            //<< " " << linAcc.x() << " " << linAcc.y() << " " << linAcc.z()
            //<< " " << angAcc.x() << " " << angAcc.y() << " " << angAcc.z() << endl; 
        }
*/
    }
}




//void Foam::access6DoFRigidBodyMotions::writeFileHeader(const label i)
void Foam::access6DoFRigidBodyMotions::writeFileHeader()
{
    // functionObjectFile.createFile() is already on master()
    // so ... parallel sync. is not allowed here


    if (Pstream::master())
    {    
    	//file() << "This is my Header" << endl;

        file().width(0); // << setw(charWidth() - 2);
        //file() << "# motion info (" << bodies[i].name() << ")" << endl;
	file() << "# motion info " << endl;        
	file() << "# time";
        file() << " surge";
        file() << " sway";
        file() << " heave";
        file() << " roll[deg]";
        file() << " pitch[deg]";
        file() << " yaw[deg]";
        //
        file() << " surge.vel";
        file() << " sway.vel";
        file() << " heave.vel";
        file() << " omega.x[rad/s]";
        file() << " omega.y[rad/s]";
        file() << " omega.z[rad/s]";
        //
        file() << " acc(surge)";
        file() << " acc(sway)";
        file() << " acc(heave)";
        file() << " acc(roll,rad/s^2)";
        file() << " acc(pitch,rad/s^2)";
        file() << " acc(yaw,rad/s^2)";
        //

        file() << endl;

      //  file() << "# time surge sway heave roll(deg) pitch(deg) yaw(deg) acc(surge) acc(sway) acc(heave) acc(roll,rad/s^2) acc(pitch,rad/s^2) acc(yaw,rad/s^2)" << endl; 
        file() << setprecision(writePrecision_);

    }
}

// ************************************************************************* //

