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
    Foam::pressureInletUniformVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This velocity inlet boundary condition is applied to patches where the
    pressure is specified.  The uniform inflow velocity is obtained by
    averaging the flux over the patch, and then applying it in the direction
    normal to the patch faces.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            pressureInletTurbulentUniformVelocity;
        U               U;
        fluctuationScale 0.1;
        alpha           0.1;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    pressureInletTurbulentUniformVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInletTurbulentUniformVelocityFvPatchVectorField_H
#define pressureInletTurbulentUniformVelocityFvPatchVectorField_H

#include "Random.H"
#include "pressureInletVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class pressureInletTurbulentUniformVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressureInletTurbulentUniformVelocityFvPatchVectorField
:
    public pressureInletVelocityFvPatchVectorField
{

    // Private data
 		
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
    TypeName("pressureInletTurbulentUniformVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureInletTurbulentUniformVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureInletTurbulentUniformVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureInletTurbulentUniformVelocityFvPatchVectorField onto a new patch
        pressureInletTurbulentUniformVelocityFvPatchVectorField
        (
            const pressureInletTurbulentUniformVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureInletTurbulentUniformVelocityFvPatchVectorField
        (
            const pressureInletTurbulentUniformVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureInletTurbulentUniformVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureInletTurbulentUniformVelocityFvPatchVectorField
        (
            const pressureInletTurbulentUniformVelocityFvPatchVectorField&,
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
                new pressureInletTurbulentUniformVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

		// Access

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
         /*  virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            ); */

        // Evaluation functions

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
