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
define(W, 500E-6)             //Square channel width
define(L,  calc(30*W))	      //Domain length
define(NXY, 25)               //Number of cells across (radial direction)
define(NZ, 300)               //Number of cells across (radial direction)
define(ExpRat, 5)           //Expansion ratio in x-direction

convertToMeters 1;

vertices        
(
	//Back face
	(0                   calc(-W/2)           0)      //00
	(calc(W/2)           calc(-W/2)           0)      //01
	(calc(W/2)           0                    0)      //02
	(0                   0                    0)      //03
	(calc(W/2)           calc(W/2)            0)      //04
	(0                   calc(W/2)            0)      //05
	//Front face
	(0                   calc(-W/2)           L)      //06
	(calc(W/2)           calc(-W/2)           L)      //07
	(calc(W/2)           0                    L)      //08
	(0                   0                    L)      //09
	(calc(W/2)           calc(W/2)            L)      //10
	(0                   calc(W/2)            L)      //11
);

blocks          
(
	hex ( 0  1  2  3  6  7  8  9) (NXY  NXY  NZ) simpleGrading (calc(1/ExpRat) ExpRat 1)  //00
	hex ( 3  2  4  5  9  8 10 11) (NXY  NXY  NZ) simpleGrading (calc(1/ExpRat) calc(1/ExpRat) 1)  //01
);

edges           
(
);

boundary         
(
	BotWall
	{
		type     wall;
		faces
		(
			( 0  1  7  6)
		);
	}

	OtherWalls
	{
		type     wall;
		faces
		(
			( 1  2  8  7)
			( 2  4 10  8)
			( 4  5 11 10)
		);
	}

	Inlet
	{
		type             patch;
		faces
		(
			( 0  3  2  1)
			( 2  3  5  4)
		);
	}

	Outlet
	{
		type             patch;
		faces
		(
			( 6  7  8  9)
			( 8 10 11  9)
		);
	}


	CenterPlane
	{
		type             symmetryPlane;
		faces
		(
			( 0  6  9  3)
			( 3  9 11  5)
		);
	}
);

mergePatchPairs 
();

// ************************************************************************* //
