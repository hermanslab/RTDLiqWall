    MRF.correctBoundaryVelocity(U);

    Hmag = -fvc::grad(UM);

    M =  (muM/muM0-1)*Hmag;
    volTensorField gradH = fvc::grad(Hmag);
    volScalarField MdH = (M & Hmag);
    Fmag =  muM0*( (gradH & M) -scalar(0.5)*fvc::grad(MdH) );
    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U) 
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U) + Fmag
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
