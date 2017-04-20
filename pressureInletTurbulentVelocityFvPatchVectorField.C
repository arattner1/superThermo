/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pressureInletTurbulentVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureInletTurbulentVelocityFvPatchVectorField::
pressureInletTurbulentVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    UName_("U"),
    ranGen_(label(0)),
    fluctuationScale_(p.size()),
    alpha_(0.1),
    curTimeIndex_(-1)
{}


Foam::pressureInletTurbulentVelocityFvPatchVectorField::
pressureInletTurbulentVelocityFvPatchVectorField
(
    const pressureInletTurbulentVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    UName_(ptf.UName_),
    ranGen_(label(0)),
    fluctuationScale_(ptf.fluctuationScale_, mapper),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


Foam::pressureInletTurbulentVelocityFvPatchVectorField::
pressureInletTurbulentVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
	UName_(dict.lookupOrDefault<word>("U", "U")),
    ranGen_(label(0)),
    fluctuationScale_("fluctuationScale", dict, p.size()),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.1)),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::pressureInletTurbulentVelocityFvPatchVectorField::
pressureInletTurbulentVelocityFvPatchVectorField
(
    const pressureInletTurbulentVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    UName_(pivpvf.UName_),
    ranGen_(pivpvf.ranGen_),
    fluctuationScale_(pivpvf.fluctuationScale_),
    alpha_(pivpvf.alpha_),
    curTimeIndex_(pivpvf.curTimeIndex_)
{}


Foam::pressureInletTurbulentVelocityFvPatchVectorField::
pressureInletTurbulentVelocityFvPatchVectorField
(
    const pressureInletTurbulentVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    UName_(pivpvf.UName_),
    ranGen_(pivpvf.ranGen_),
    fluctuationScale_(pivpvf.fluctuationScale_),
    alpha_(pivpvf.alpha_),
    curTimeIndex_(pivpvf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureInletTurbulentVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    tmp<vectorField> n = patch().nf();
    const Field<scalar>& magS = patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(n*phip/magS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(n*phip/(rhop*magS));
    }
    else
    {
        FatalErrorIn("pressureInletTurbulentVelocityFvPatchVectorField::updateCoeffs()")
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

	if (curTimeIndex_ != this->db().time().timeIndex())
    {

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

        vectorField& patchField = *this;

        vectorField randomField(this->size());

        forAll(patchField, facei)
        {
            ranGen_.randomise(randomField[facei]);
        }

        // Correction-factor to compensate for the loss of RMS fluctuation
        // due to the temporal correlation introduced by the alpha parameter.
        scalar rmsCorr = sqrt(12*(2*alpha_ - sqr(alpha_)))/alpha_;

        patchField =
            (1 - alpha_)*patchField
          + alpha_*
            (
                Up
              + rmsCorr*cmptMultiply
                (
                    randomField - 0.5*pTraits<vector>::one,
                    fluctuationScale_
                )*mag(Up)
            );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::pressureInletTurbulentVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    os.writeKeyword("fluctuationScale")
    << fluctuationScale_ << token::END_STATEMENT << nl;
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureInletTurbulentVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureInletTurbulentVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
