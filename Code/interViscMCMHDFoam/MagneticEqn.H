
      	    // solve
            fvScalarMatrix UMEqn
            (
                fvm::laplacian(muM,UM)
            );
            UMEqn.solve();

