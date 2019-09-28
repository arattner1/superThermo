/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.0                                   |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;
    instance        "system";
    class           dictionary;
    object          perturbUDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Expected turbulent Re(tau) Reynolds number with respect to utau
Retau 572.7;

// Streamwise component of flow. 0=x, 1=y, 2=z
streamwise 2;

// Spanwise component of flow. 0=x, 1=y, 2=z
spanwise 1;

// Halfheight of channel. This is the direction normal to both streamwise and
// spanwise directions.
h 368.66e-6;

// Set (overwrite) velocity to laminar profile
setBulk true;

// Perturb velocity with some cosine like perturbations
perturb true;

// Perturbation properties
//wall normal circulation as a fraction of Ubar/utau
duplus 0.25;
//spanwise perturbation spacing in wall units
betaPlus 100;
//streamwise perturbation spacing in wall units
alphaPlus 500;
//transverse decay
sigma 0.00055;
//linear perturbation amplitude as a fraction of Ubar
epsilon 0.05;


// ************************************************************************* //
