
      {
        
        D=(Dw * Do) / (alpha1 * Do + (1 - alpha1) * Dw);
        H=(1 - K) / (alpha1 + K * (1 - alpha1));
        
        // Standard face-flux compression coefficient
        surfaceScalarField phic = mag(phi/mesh.magSf());
        phic = min(mixture.cAlpha()*phic, max(phic));
        
        // Add the optional isotropic compression contribution
        if (icAlpha > 0)
        {
          phic *= (1.0 - icAlpha);
          phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
        }
        
        // Add the optional shear compression contribution
        if (scAlpha > 0)
        {
          phic +=
          scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
        }
        
        
        surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();
        
        // Do not compress interface at non-coupled boundary faces
        // (inlets, outlets etc.)
        forAll(phic.boundaryField(), patchi)
        {
          fvsPatchScalarField& phicp = phicBf[patchi];
          
          if (!phicp.coupled())
          {
            phicp == 0;
          }
        }
        
        surfaceScalarField phir = phic*mixture.nHatf();
        
        fvScalarMatrix c1Eqn
        (
          fvm::ddt(c1)
          + fvm::div(phi, c1)
          - fvm::laplacian(D, c1)
          + fvc::laplacian(D*H*c1, alpha1)
          + fvc::div(phir, H*c1*alpha1*(1-alpha1))
        );
        
        c1Eqn.relax();
        c1Eqn.solve();
        
        fvScalarMatrix c2Eqn
        (
          fvm::ddt(c2)
          + fvm::div(phi, c2)
          - fvm::laplacian(D, c2)
          + fvc::laplacian(D*H*c2, alpha1)
          + fvc::div(phir, H*c2*alpha1*(1-alpha1))
        );
        
        c2Eqn.relax();
        c2Eqn.solve();
        
        fvScalarMatrix c3Eqn
        (
          fvm::ddt(c3)
          + fvm::div(phi, c3)
          - fvm::laplacian(D, c3)
          + fvc::laplacian(D*H*c3, alpha1)
          + fvc::div(phir, H*c3*alpha1*(1-alpha1))
        );
        
        c3Eqn.relax();
        c3Eqn.solve();
        
        double dc1,dc2,dc3;
        
        forAll(mesh.cells(), i) {
          
          dc1 = (-k1.value()*c1[i]*c2[i]+km1.value()*c3[i])*runTime.deltaTValue();
          dc2 = (-k1.value()*c1[i]*c2[i]+km1.value()*c3[i])*runTime.deltaTValue();
          dc3 = (k1.value()*c1[i]*c2[i]-km1.value()*c3[i])*runTime.deltaTValue();
          
          c1[i] += dc1;
          c2[i] += dc2;
          c3[i] += dc3;
        }
        
      }
      
