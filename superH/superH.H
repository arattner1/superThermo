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
    Foam::superH

Description
    Generic enthalpy definition for supercritical fluids

SourceFiles
    superHI.H
    superH.C

\*---------------------------------------------------------------------------*/

#ifndef superH_H
#define superH_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class superH;

template<class EquationOfState>
inline superH<EquationOfState> operator+
(
    const superH<EquationOfState>&,
    const superH<EquationOfState>&
);

template<class EquationOfState>
inline superH<EquationOfState> operator-
(
    const superH<EquationOfState>&,
    const superH<EquationOfState>&
);

template<class EquationOfState>
inline superH<EquationOfState> operator*
(
    const scalar,
    const superH<EquationOfState>&
);

template<class EquationOfState>
inline superH<EquationOfState> operator==
(
    const superH<EquationOfState>&,
    const superH<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const superH<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class superH Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class superH
:
    public EquationOfState
{
    // Private data

    scalar cpL;
    scalar cpG;
    scalar TSpike;
    scalar SpikeWidth1;
    scalar cpSpike1;
    scalar SpikeWidth2;
    scalar cpSpike2;
    scalar SpikeWidth3;
    scalar cpSpike3;
    scalar SpikeWidth4;
    scalar cpSpike4;


    // Private Member Functions

    //- Construct from components
    inline superH
    (
        const EquationOfState& st,
        const scalar cpL,
        const scalar cpG,
        const scalar TSpike,
        const scalar SpikeWidth1,
        const scalar cpSpike1,
        const scalar SpikeWidth2,
        const scalar cpSpike2,
        const scalar SpikeWidth3,
        const scalar cpSpike3,
        const scalar SpikeWidth4,
        const scalar cpSpike4
    );


public:

    // Constructors

        //- Construct from Istream
        superH(Istream&);

        //- Construct from dictionary
        superH(const dictionary& dict);

        //- Construct as named copy
        inline superH(const word&, const superH&);

        //- Construct and return a clone
        inline autoPtr<superH> clone() const;

        //- Selector from Istream
        inline static autoPtr<superH> New(Istream& is);

        //- Selector from dictionary
        inline static autoPtr<superH> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "superH<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const superH&);
        inline void operator-=(const superH&);


    // Friend operators

        friend superH operator+ <EquationOfState>
        (
            const superH&,
            const superH&
        );

        friend superH operator- <EquationOfState>
        (
            const superH&,
            const superH&
        );

        friend superH operator* <EquationOfState>
        (
            const scalar,
            const superH&
        );

        friend superH operator== <EquationOfState>
        (
            const superH&,
            const superH&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const superH&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "superHI.H"

#ifdef NoRepository
#   include "superH.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
