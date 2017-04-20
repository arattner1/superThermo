/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::pressureDirectedTurbulentInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This velocity inlet boundary condition is applied to patches where the
    pressure is specified.  The inflow velocity is obtained from the flux with
    the specified inlet direction" direction.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        inletDirection | inlet direction per patch face | yes |
        U            | velocity field name     | no          | U
        fluctuationScale | RMS fluctuation scale (fraction of mean) | yes |
        alpha | fraction of new random component added to previous| no| 0.1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            pressureDirectedTurbulentInletVelocity;
        phi             phi;
        U               U;
        rho             rho;
        inletDirection  uniform (1 0 0);
        value           uniform 0;
        fluctuationScale 0.1;
        alpha           0.1;
    }
    \endverbatim

Note
    If reverse flow is possible or expected use the 
    pressureDirectedInletOutletVelocityFvPatchVectorField condition instead.


SeeAlso
    Foam::fixedValueFvPatchField
    Foam::pressureDirectedInletOutletVelocityFvPatchVectorField

SourceFiles
    pressureDirectedTurbulentInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureDirectedTurbulentInletVelocityFvPatchVectorField_H
#define pressureDirectedTurbulentInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class pressureDirectedTurbulentInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressureDirectedTurbulentInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Flux field name
        word phiName_;

        //- Density field name
        word rhoName_;

        //- Inlet direction
        vectorField inletDir_;

        //- Random number generator
        Random ranGen_;

        //- Name of the velocity field
        word UName_;

        //- Fluctuation scake
        vectorField fluctuationScale_;

        //- Fraction of RMS component to apply to last time step values
        scalar alpha_;

        //- Current time index (used for updating)
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("pressureDirectedTurbulentInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureDirectedTurbulentInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureDirectedTurbulentInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureDirectedTurbulentInletVelocityFvPatchVectorField
        //  onto a new patch
        pressureDirectedTurbulentInletVelocityFvPatchVectorField
        (
            const pressureDirectedTurbulentInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureDirectedTurbulentInletVelocityFvPatchVectorField
        (
            const pressureDirectedTurbulentInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureDirectedTurbulentInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureDirectedTurbulentInletVelocityFvPatchVectorField
        (
            const pressureDirectedTurbulentInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureDirectedTurbulentInletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the name of rho
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Return the reference field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the reference field to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the fluctuation scale
            const vectorField& fluctuationScale() const
            {
                return fluctuationScale_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            vectorField& fluctuationScale()
            {
                return fluctuationScale_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
