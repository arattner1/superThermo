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
define(W, 750E-6)             //Square channel width
define(L1,  30E-3)	      //Development length
define(L2,  50E-3)	      //Total Domain length
define(NXY, 30)               //Number of cells across (radial direction)
define(NZ, 1000)               //Number of cells across (radial direction)
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
	//Middle face
	(0                   calc(-W/2)          L1)      //06
	(calc(W/2)           calc(-W/2)          L1)      //07
	(calc(W/2)           0                   L1)      //08
	(0                   0                   L1)      //09
	(calc(W/2)           calc(W/2)           L1)      //10
	(0                   calc(W/2)           L1)      //11
	//Front face
	(0                   calc(-W/2)           L2)      //12
	(calc(W/2)           calc(-W/2)           L2)      //13
	(calc(W/2)           0                    L2)      //14
	(0                   0                    L2)      //15
	(calc(W/2)           calc(W/2)            L2)      //16
	(0                   calc(W/2)            L2)      //17
);

blocks          
(
	hex ( 0  1  2  3  6  7  8  9) (NXY  NXY  NZ) simpleGrading (calc(1/ExpRat) ExpRat 1)  //00
	hex ( 3  2  4  5  9  8 10 11) (NXY  NXY  NZ) simpleGrading (calc(1/ExpRat) calc(1/ExpRat) 1)  //01
	hex ( 6  7  8  9 12 13 14 15) (NXY  NXY  NZ) simpleGrading (calc(1/ExpRat) ExpRat 1)  //02
	hex ( 9  8 10 11 15 14 16 17) (NXY  NXY  NZ) simpleGrading (calc(1/ExpRat) calc(1/ExpRat) 1)  //03
);

edges           
(
);

boundary         
(
	BotWall1
	{
		type     wall;
		faces
		(
			( 0  1  7 6)
		);
	}

	BotWall2
	{
		type     wall;
		faces
		(
			( 6 7 13 12)
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
			( 7  8 14 13)
			( 8  10 16 14)
			( 10 11 17 16)
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
			( 12 13 14 15)
			( 14 16 17 15)
		);
	}


	CenterPlane
	{
		type             symmetryPlane;
		faces
		(
			( 0  6  9  3)
			( 3  9 11  5)
			( 6 12 15 9)
			(9 15 17 11)
		);
	}
);

mergePatchPairs 
();

// ************************************************************************* //
