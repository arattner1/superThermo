{
    bool closedVolume = p_rgh.needReference();
    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);

    rho = thermo.rho();

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
      + phig
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    tmp<fvScalarMatrix> p_rghDDtEqn
    (
        new fvScalarMatrix(p_rgh, dimMass/dimTime)
    );

    {
        if (compressible)
        {
            p_rghDDtEqn =
            (
                fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
             ==
                fvOptions(psi, p_rgh, rho.name())
            );
        }

        // Thermodynamic density needs to be updated by psi*d(p) after the
        // pressure solution - done in 2 parts. Part 1:
        thermo.rho() -= psi*p_rgh;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix p_rghEqn
            (
                p_rghDDtEqn()
              + fvc::div(phiHbyA)
              - fvm::laplacian(rhorAUf, p_rgh)
            );

            p_rghEqn.solve
            (
                mesh.solver
                (
                    p_rgh.select
                    (
                        (
                           oCorr == nOuterCorr-1
                        && corr == nCorr-1
                        && nonOrth == nNonOrthCorr
                        )
                    )
                )
            );

            if (nonOrth == nNonOrthCorr)
            {
                phi = phiHbyA + p_rghEqn.flux();
                U = HbyA
                  + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
                U.correctBoundaryConditions();
            //    fvOptions.correct(U);

                   rAGave = 0.0;
                   magGbarAve = 0.0;

                   rAGave = gSum( ( 1.0/UEqn.A() ) * rho.field() * mesh.V() ) / gSum( mesh.V() );

	           GbarAve = gSum( rho.field() * U.field() * mesh.V() ) / gSum( mesh.V() );
		   magGbarAve = mag(GbarAve[2]);

 forAll(fluidRegions, i)
 {
		 IOdictionary GAdaptDictionary
		(
		    IOobject
		    (
			"GAdaptDictionary",
			runTime.constant(),
			fluidRegions[i],
			IOobject::MUST_READ_IF_MODIFIED,
			IOobject::NO_WRITE
		    )
		);

		//Target value of G
		const scalar Gbar( readScalar(GAdaptDictionary.lookup("Gbar")) );

		//Over-relaxation parameter
		const scalar GRelax ( readScalar(GAdaptDictionary.lookup("Relaxation")) );

		dGradP = GRelax*(mag(Gbar) - magGbarAve)/rAGave;
}

		   //correct U
                   scalarField rAU2 = ( 1.0/UEqn.A() * mesh.V() ) / ( mesh.V());
		   U.field() += flowDir*rAU2*dGradP;
		   U.correctBoundaryConditions();

		   gradP = gradP0 + dGradP;   
		    
	           dimensionedVector SuUpdate("SuUpdate", dimensionSet( 1, -2, -2, 0, 0, 0, 0 ), flowDir*gradP ); 
		   Su = SuUpdate;  

                K = 0.5*magSqr(U);
            }
        }

        // Second part of thermodynamic density update
        thermo.rho() += psi*p_rgh;
    }

    p = p_rgh + rho*gh;

    // Update pressure time derivative if needed
    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }

    if (compressible)
    {
        // Solve continuity
        #include "rhoEqn.H"
    }

    // Update continuity errors
    #include "compressibleContinuityErrors.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume && compressible)
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /compressibility;
        rho = thermo.rho();
        p_rgh = p - rho*gh;
    }
}
