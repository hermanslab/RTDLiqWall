#include "createRDeltaT.H"

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

    Info<< "Reading field Magnetic Potential UM\n" << endl;
    volScalarField UM
    (
        IOobject
        (
            "UM",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

Info<< "    Reading magnetic properties\n" << nl << endl;
    IOdictionary elPropDict
    (
        IOobject
        (
            "magneticProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
        dimensionedScalar muM1
(
    elPropDict.subDict("phase1").lookup("muM")
);
        dimensionedScalar muM2
(
    elPropDict.subDict("phase2").lookup("muM")
);

        dimensionedScalar Ms
(
    elPropDict.subDict("phase2").lookup("Ms")
);

        dimensionedScalar muM0
(
    elPropDict.subDict("phase2").lookup("muM0")
);

    
    Info<<"Found: muM1=" << muM1 << ",muM2=" << muM2 << endl;
    Info<<"Found: muM0=" << muM0 << "Ms=" << Ms << endl;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

    volScalarField muM
    (
        IOobject
        (
            "muM",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        muM1*muM2/(alpha1*muM2+ (1.0 - alpha1)*muM1)
    );
    muM.oldTime();

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


       IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    

dimensionedScalar Dw
(
    transportProperties.lookup("Dw")
);

dimensionedScalar Do
(
    transportProperties.lookup("Do")
);

dimensionedScalar K
(
    transportProperties.lookup("K")
);



    IOdictionary kineticProperties
(
    IOobject
    (
        "kineticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
      
dimensionedScalar k1
(
    kineticProperties.lookup("k1")
);
    

dimensionedScalar km1
(
    kineticProperties.lookup("km1")
);
  

Info<< "Reading field c1\n" << endl;
volScalarField c1
(
    IOobject
    (
        "c1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field c2\n" << endl;
volScalarField c2
(
    IOobject
    (
        "c2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field c3\n" << endl;
volScalarField c3
(
    IOobject
    (
        "c3",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (1 - K) / (alpha1 + K * (1 - alpha1))
);


surfaceScalarField phir
(
    IOobject
    (
        "phir",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);
    
volScalarField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (Dw * Do) / (alpha1 * Do + (1 - alpha1) * Dw)
);

volVectorField Hmag
(
    IOobject
    (
        "Hmag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(UM)
);

volVectorField M
(
    IOobject
    (
        "M",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (muM/muM0-1)*Hmag
);



volVectorField Fmag
(
    IOobject
    (
        "Fmag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   muM0*( (fvc::grad(Hmag) & M) -scalar(0.5)*fvc::grad(M & Hmag) )
);


typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


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
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"
