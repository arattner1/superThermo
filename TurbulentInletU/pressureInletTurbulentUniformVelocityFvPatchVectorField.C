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

#include "pressureInletTurbulentUniformVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::
pressureInletTurbulentUniformVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletVelocityFvPatchVectorField(p, iF),
    UName_("U"),
    ranGen_(label(0)),
    fluctuationScale_(p.size()),
    alpha_(0.1),
    curTimeIndex_(-1)
{}


Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::
pressureInletTurbulentUniformVelocityFvPatchVectorField
(
    const pressureInletTurbulentUniformVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    pressureInletVelocityFvPatchVectorField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    ranGen_(label(0)),
    fluctuationScale_(ptf.fluctuationScale_, mapper),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::
pressureInletTurbulentUniformVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    pressureInletVelocityFvPatchVectorField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    ranGen_(label(0)),
    fluctuationScale_("fluctuationScale", dict, p.size()),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.1)),
    curTimeIndex_(-1)
{}


Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::
pressureInletTurbulentUniformVelocityFvPatchVectorField
(
    const pressureInletTurbulentUniformVelocityFvPatchVectorField& pivpvf
)
:
    pressureInletVelocityFvPatchVectorField(pivpvf),
    UName_(pivpvf.UName_),
    ranGen_(pivpvf.ranGen_),
    fluctuationScale_(pivpvf.fluctuationScale_),
    alpha_(pivpvf.alpha_),
    curTimeIndex_(pivpvf.curTimeIndex_)
{}


Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::
pressureInletTurbulentUniformVelocityFvPatchVectorField
(
    const pressureInletTurbulentUniformVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletVelocityFvPatchVectorField(pivpvf, iF),
    UName_(pivpvf.UName_),
    ranGen_(pivpvf.ranGen_),
    fluctuationScale_(pivpvf.fluctuationScale_),
    alpha_(pivpvf.alpha_),
    curTimeIndex_(pivpvf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/* void Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    pressureInletVelocityFvPatchVectorField::autoMap(m);
    fluctuationScale_.autoMap(m);
}


void Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    pressureInletVelocityFvPatchVectorField::rmap(ptf, addr);

    const pressureInletTurbulentUniformVelocityFvPatchVectorField& tiptf =
        refCast<const pressureInletTurbulentUniformVelocityFvPatchVectorField >(ptf);

    fluctuationScale_.rmap(tiptf.fluctuationScale_, addr);
} */


void Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	if (curTimeIndex_ != this->db().time().timeIndex())
    {

    pressureInletVelocityFvPatchVectorField::updateCoeffs();

    operator==(patch().nf()*gSum(patch().Sf() & *this)/gSum(patch().magSf()));

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

}


void Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
	this->writeEntry("fluctuationScale", os);
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureInletTurbulentUniformVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    operator==(patch().nf()*gSum(patch().Sf() & pvf)/gSum(patch().magSf()));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureInletTurbulentUniformVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
