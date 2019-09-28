     if ( n > 1.0)
     {
       gradP0 += dGradP;
     } 

   // Solve the modified Momentum equation
    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turb.divDevRhoReff(U)
     ==
        fvOptions(rho, U) + Su
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

 //   fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            ),
            mesh.solver(U.select(finalIter))
        );

    }

   rAGave = gSum( ( 1.0/UEqn.A() ) * rho.field() * mesh.V() ) / gSum( mesh.V() );

   vector Ubar = gSum( U.field() * mesh.V() ) / gSum( mesh.V() );
   vector GbarAve = gSum( rho.field() * U.field() * mesh.V() ) / gSum( mesh.V() );
   scalar magGbarAve = mag(GbarAve);

     vector flowDir = Ubar/mag(Ubar);

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
   
   gradP = gradP0 + dGradP;   
 
 forAll(fluidRegions, i)
 {
	dimensionedScalar gradPFieldUpdate("gradPFieldUpdate",dimensionSet(1,-2,-2,0,0,0,0),gradP);
   	gradPField[i] = gradPFieldUpdate; 
 }
   K = 0.5*magSqr(U);

   n = n + 1.0; 

