if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
dimensionedScalar gamma = 3.0*(muM2/muM0-1.0)/Ms; //m/A
dimensionedScalar SMALLH = 1e-20*Ms;

muM == muM1*(((Ms/(SMALLH+mag(fvc::grad(UM)))*(1.0/tanh(SMALL+gamma*mag(fvc::grad(UM)))-1.0/(SMALL+gamma*mag(fvc::grad(UM)))))+1.0)*muM0)/(alpha1*(((Ms/(SMALLH+mag(fvc::grad(UM)))*(1.0/tanh(SMALL+gamma*mag(fvc::grad(UM)))-1.0/(SMALL+gamma*mag(fvc::grad(UM)))))+1.0)*muM0)+ (1.0 - alpha1)*muM1);



