
     if ( n > 1.0)
     {
       gradP0 += dGradP;
     }
  
   // Solve the modified Momentum equation
    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U) 
     ==
        fvOptions(rho, U) + Su
    );

    UEqn.relax();


//  fvOptions.constrain(UEqn);

  if (pimple.momentumPredictor())
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
            ) 
        );
    } 

   rAGave = gSum( ( 1.0/UEqn.A() ) * rho.field() * mesh.V() ) / gSum( mesh.V() );

   vector Ubar = gSum( U.field() * mesh.V() ) / gSum( mesh.V() );
   vector GbarAve = gSum( rho.field() * U.field() * mesh.V() ) / gSum( mesh.V() );
   scalar magGbarAve = mag(GbarAve);

   vector flowDir = Ubar/mag(Ubar);
   //vector flowDir = vector(1.0, 0.0, 1.0); 

   dGradP = GRelax*(mag(Gbar) - magGbarAve)/rAGave;

   //correct U
   scalarField rAU2 = ( 1.0/UEqn.A() * mesh.V() ) / ( mesh.V());
   U.field() += flowDir*rAU2*dGradP;
   
   gradP = gradP0 + dGradP;   

   gradPField.field() = gradP; 
 
   K = 0.5*magSqr(U);
   n = n + 1.0; 


