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

#include "pressureDirectedTurbulentInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::
pressureDirectedTurbulentInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    inletDir_(p.size()),
    referenceField_(p.size()),
    ranGen_(label(0)),
    fluctuationScale_(p.size()),
    alpha_(0.1),
    curTimeIndex_(-1)
{}


Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::
pressureDirectedTurbulentInletVelocityFvPatchVectorField
(
    const pressureDirectedTurbulentInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    inletDir_(ptf.inletDir_, mapper),
    referenceField_(ptf.referenceField_, mapper),
    ranGen_(label(0)),
    fluctuationScale_(ptf.fluctuationScale_, mapper),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1)
{}


Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::
pressureDirectedTurbulentInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    inletDir_("inletDirection", dict, p.size()),
    referenceField_("referenceField", dict, p.size()),
    ranGen_(label(0)),
    fluctuationScale_("fluctuationScale", dict, p.size()),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.1)),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
    	fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
    	fvPatchVectorField::operator=(vectorField("referenceField", dict, p.size()));
    }

}


Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::
pressureDirectedTurbulentInletVelocityFvPatchVectorField
(
    const pressureDirectedTurbulentInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_),
    referenceField_(pivpvf.referenceField_),
    ranGen_(pivpvf.ranGen_),
    fluctuationScale_(pivpvf.fluctuationScale_),
    alpha_(pivpvf.alpha_),
    curTimeIndex_(pivpvf.curTimeIndex_)
{}


Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::
pressureDirectedTurbulentInletVelocityFvPatchVectorField
(
    const pressureDirectedTurbulentInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    inletDir_(pivpvf.inletDir_),
    referenceField_(pivpvf.referenceField_),
    ranGen_(pivpvf.ranGen_),
    fluctuationScale_(pivpvf.fluctuationScale_),
    alpha_(pivpvf.alpha_),
    curTimeIndex_(pivpvf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    inletDir_.autoMap(m);
    referenceField_.autoMap(m);
}


void Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const pressureDirectedTurbulentInletVelocityFvPatchVectorField& tiptf =
        refCast<const pressureDirectedTurbulentInletVelocityFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
    referenceField_.rmap(tiptf.referenceField_, addr);
}


void Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::updateCoeffs()
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
    tmp<scalarField> ndmagS = (n & inletDir_)*patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(inletDir_*phip/ndmagS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(inletDir_*phip/(rhop*ndmagS));
    }
    else
    {
        FatalErrorIn
        (
            "pressureDirectedTurbulentInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

	if (curTimeIndex_ != this->db().time().timeIndex())
    {
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
                referenceField_
              + rmsCorr*cmptMultiply
                (
                    randomField - 0.5*pTraits<vector>::one,
                    fluctuationScale_
                )*mag(referenceField_)
            );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    inletDir_.writeEntry("inletDirection", os);
    os.writeKeyword("fluctuationScale")
    << fluctuationScale_ << token::END_STATEMENT << nl;
    referenceField_.writeEntry("referenceField", os);
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(inletDir_*(inletDir_ & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureDirectedTurbulentInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
