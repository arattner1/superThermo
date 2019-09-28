/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::PengRobinsonGasPConst

Group
    grpSpecieEquationOfState

Description
    PengRobinsonGasPConst gas equation of state.
    Alex Rattner, 2017-06-27: Variant which assumes an effectively constant 
    pressure to stabilize solves for near-critical-point simulations

SourceFiles
    PengRobinsonGasPConstI.H
    PengRobinsonGasPConst.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonGasPConst_H
#define PengRobinsonGasPConst_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class PengRobinsonGasPConst;

template<class Specie>
inline PengRobinsonGasPConst<Specie> operator+
(
    const PengRobinsonGasPConst<Specie>&,
    const PengRobinsonGasPConst<Specie>&
);

template<class Specie>
inline PengRobinsonGasPConst<Specie> operator-
(
    const PengRobinsonGasPConst<Specie>&,
    const PengRobinsonGasPConst<Specie>&
);

template<class Specie>
inline PengRobinsonGasPConst<Specie> operator*
(
    const scalar,
    const PengRobinsonGasPConst<Specie>&
);

template<class Specie>
inline PengRobinsonGasPConst<Specie> operator==
(
    const PengRobinsonGasPConst<Specie>&,
    const PengRobinsonGasPConst<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const PengRobinsonGasPConst<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class PengRobinsonGasPConst Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class PengRobinsonGasPConst
:
    public Specie
{
    // Private data

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;
        
        //- Assumed constant pressure for density calculations
        scalar PConst_;


public:

    // Constructors

        //- Construct from components
        inline PengRobinsonGasPConst
        (
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Zc,
            const scalar& Pc,
            const scalar& omega,
            const scalar& PConst
        );

        //- Construct from Istream
        PengRobinsonGasPConst(Istream&);

        //- Construct from dictionary
        PengRobinsonGasPConst(const dictionary& dict);

        //- Construct as named copy
        inline PengRobinsonGasPConst(const word& name, const PengRobinsonGasPConst&);

        //- Construct and return a clone
        inline autoPtr<PengRobinsonGasPConst> clone() const;

        // Selector from Istream
        inline static autoPtr<PengRobinsonGasPConst> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<PengRobinsonGasPConst> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "PengRobinsonGasPConst<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kmol]
            inline scalar h(const scalar p, const scalar T) const;

            //- Return cp departure [J/(kmol K]
            inline scalar cp(scalar p, scalar T) const;

            //- Return entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor [-]
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const PengRobinsonGasPConst&);
        inline void operator-=(const PengRobinsonGasPConst&);

        inline void operator*=(const scalar);


    // Friend operators

        friend PengRobinsonGasPConst operator+ <Specie>
        (
            const PengRobinsonGasPConst&,
            const PengRobinsonGasPConst&
        );

        friend PengRobinsonGasPConst operator- <Specie>
        (
            const PengRobinsonGasPConst&,
            const PengRobinsonGasPConst&
        );

        friend PengRobinsonGasPConst operator* <Specie>
        (
            const scalar s,
            const PengRobinsonGasPConst&
        );

        friend PengRobinsonGasPConst operator== <Specie>
        (
            const PengRobinsonGasPConst&,
            const PengRobinsonGasPConst&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const PengRobinsonGasPConst&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PengRobinsonGasPConstI.H"

#ifdef NoRepository
    #include "PengRobinsonGasPConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
