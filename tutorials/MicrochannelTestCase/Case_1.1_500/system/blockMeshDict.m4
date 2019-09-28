/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.7.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//For m4 preprocessing
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(roundcalc, [esyscmd(perl -e 'printf ( int( ($1) + 0.5 ) )')])

//Some dimensions here
define(W, 750.02E-6)             //Square channel width
define(H, 737.32E-6)             //Square channel height
define(L1,  2.5E-3)	      //Development length
define(L2,  5E-3)	      //Total Domain length
define(NXY, 48)                  //Number of cells across (radial direction)
define(NZ1, 48)               //Number of cells across unheated section (radial direction)
define(NZ2, 48)               //Number of cells across heated section (radial direction)
define(ExpRat, 20)           //Expansion ratio normal to the walls

convertToMeters 1;

vertices        
(
	//Back face
	(0          0            0)     //00
	(W          0            0)     //01
	(W          H            0)     //02
	(0          H            0)     //03
	//Mid face
	(0          0            L1)     //04
	(W          0            L1)     //05
	(W          H            L1)     //06
	(0          H            L1)     //07
    //Front face
	(0          0            L2)     //08
	(W          0            L2)     //09
	(W          H            L2)     //10
	(0          H            L2)     //11
);

blocks          
(
	hex ( 0  1  2  3  4  5  6  7) (calc(2*NXY)  calc(2*NXY)  NZ1) 
	    simpleGrading ( ( (0.5 0.5 ExpRat) (0.5 0.5 calc(1/ExpRat)) ) 
	                    ( (0.5 0.5 ExpRat) (0.5 0.5 calc(1/ExpRat)) )
	                    1 )
	hex ( 4  5  6  7  8  9 10 11) (calc(2*NXY)  calc(2*NXY)  NZ2) 
	    simpleGrading ( ( (0.5 0.5 ExpRat) (0.5 0.5 calc(1/ExpRat)) ) 
	                    ( (0.5 0.5 ExpRat) (0.5 0.5 calc(1/ExpRat)) )
	                    1 )                    
);

edges           
(
);

boundary         
(
	BottWalls
	{
		type     wall;
		faces
		(
		    ( 0  4  5  1)
		    ( 4  8  9  5)
		);
	}

        SideWalls
	{
		type     wall;
		faces
		(
		    ( 1  5  6  2)
  		    ( 5  9 10  6)
		    ( 0  3  7  4)
  		    ( 4  7 11  8)
		);
	}

	TopWalls
	{
		type     wall;
		faces
		(
		    ( 2  6  7  3)
   		    ( 6 10 11  7)
		);
	}

	Inlet
	{
		type            cyclic;
		neighbourPatch  Outlet;
	    	faces
		(   ( 0  3  2  1) 	);
	}

	Outlet
	{
		type            cyclic;
		neighbourPatch  Inlet;
		faces
		(	( 8 11 10  9)	);
	}

);

mergePatchPairs 
();

// ************************************************************************* //
