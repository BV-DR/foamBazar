Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Modifs Charles: in this way, we don't have to provide both p and p_rgh as initial condition (p is defined from p_rgh)

Info<< "Creating P from P_rgh \n" << endl;
volScalarField p
(
   IOobject
   (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
   ),
   p_rgh //+ rho*gh
);


{  // what is in the scope is only local (except for variables already defined outside the scope -> p)!
    twoPhaseMixture tmp
    (
        mesh,
        IOdictionary
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    volScalarField& alpha1(tmp.alpha1());
    volScalarField& alpha2(tmp.alpha2());

    IOdictionary tmpDict
    (
       IOobject
       	(
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
        )
	);	

	dimensioned<scalar> rho1(tmpDict.subDict(tmp.phase1Name()).lookup("rho"));	
	dimensioned<scalar> rho2(tmpDict.subDict(tmp.phase2Name()).lookup("rho"));	

	#include "readGravitationalAcceleration.H"
	#include "readhRef.H"
	#include "gh.H"

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    );

    p = p_rgh + rho*gh;
}

p.write();
// End modifs



#include "createPhi.H"


Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(mesh);
mixture.thermo1().rho().write();


volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

//volScalarField& p = mixture.p(); // p is already defined correctly
volScalarField& T = mixture.T();
volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& psi1 = mixture.thermo1().psi();
volScalarField& rho2 = mixture.thermo2().rho();
const volScalarField& psi2 = mixture.thermo2().psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);


dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"



// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt
(
    pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
);

// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture);

// Construct compressible turbulence model
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New(rho, U, rhoPhi, mixture)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));
