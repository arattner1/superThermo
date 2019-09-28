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
    Foam::superEOS

Description
    Generalized supercritical EOS.

SourceFiles
    superEOSI.H
    superEOS.C

\*---------------------------------------------------------------------------*/

#ifndef superEOS_H
#define superEOS_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class superEOS;

template<class Specie>
inline superEOS<Specie> operator+
(
    const superEOS<Specie>&,
    const superEOS<Specie>&
);

template<class Specie>
inline superEOS<Specie> operator-
(
    const superEOS<Specie>&,
    const superEOS<Specie>&
);

template<class Specie>
inline superEOS<Specie> operator*
(
    const scalar,
    const superEOS<Specie>&
);

template<class Specie>
inline superEOS<Specie> operator==
(
    const superEOS<Specie>&,
    const superEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const superEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class superEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class superEOS
:
    public Specie
{

    // Private data

        //- Density
        scalar rhoL_;
        scalar rhoG_;
        scalar psiL_;
        scalar psiG_;
        scalar rhoVariationFactor_;
		scalar psiSpike_;
		scalar SpikeWidth_;
        scalar Tband_;
        scalar Pband_;

public:

    // Constructors

        //- Construct from components
        inline superEOS( const Specie& sp,
                         const scalar rhoL,
                         const scalar rhoG,
                         const scalar psiL,
                         const scalar psiG,
                         const scalar rhoVariationFactor,
                         const scalar psiSpike,
                         const scalar SpikeWidth,
                         const scalar Tband,
                         const scalar Pband );
        
        
        //- Construct from Istream
        superEOS(Istream&);

        //- Construct from dictionary
        superEOS(const dictionary& dict);

        //- Construct as named copy
        inline superEOS(const word& name, const superEOS&);

        //- Construct and return a clone
        inline autoPtr<superEOS> clone() const;

        // Selector from Istream
        inline static autoPtr<superEOS> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<superEOS> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "superEOS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const superEOS&);
        inline void operator-=(const superEOS&);

        inline void operator*=(const scalar);


    // Friend operators

        friend superEOS operator+ <Specie>
        (
            const superEOS&,
            const superEOS&
        );

        friend superEOS operator- <Specie>
        (
            const superEOS&,
            const superEOS&
        );

        friend superEOS operator* <Specie>
        (
            const scalar s,
            const superEOS&
        );

        friend superEOS operator== <Specie>
        (
            const superEOS&,
            const superEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const superEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "superEOSI.H"

#ifdef NoRepository
#   include "superEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
