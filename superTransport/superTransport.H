/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::superTransport

Description
    superCritical properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    superTransportI.H
    superTransport.C

\*---------------------------------------------------------------------------*/

#ifndef superTransport_H
#define superTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class superTransport;

template<class Thermo>
inline superTransport<Thermo> operator+
(
    const superTransport<Thermo>&,
    const superTransport<Thermo>&
);

template<class Thermo>
inline superTransport<Thermo> operator-
(
    const superTransport<Thermo>&,
    const superTransport<Thermo>&
);

template<class Thermo>
inline superTransport<Thermo> operator*
(
    const scalar,
    const superTransport<Thermo>&
);

template<class Thermo>
inline superTransport<Thermo> operator==
(
    const superTransport<Thermo>&,
    const superTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const superTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class superTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class superTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
		  scalar muLcoeff0_;
		  scalar muLcoeff1_;
		  scalar muLcoeff2_;
		  scalar muLcoeff3_;
		  scalar muLcoeff4_;
		  scalar muLcoeff5_;
		  scalar muVcoeff0_;
		  scalar muVcoeff1_;
		  scalar muVcoeff2_;
		  scalar muVcoeff3_;
		  scalar muVcoeff4_;
		  scalar muVcoeff5_;

        //- Reciprocal Prandtl Number []
		  scalar PrLcoeff0_;
		  scalar PrLcoeff1_;
		  scalar PrLcoeff2_;
		  scalar PrLcoeff3_;
		  scalar PrLcoeff4_;
		  scalar PrLcoeff5_;
		  scalar PrVcoeff0_;
		  scalar PrVcoeff1_;
		  scalar PrVcoeff2_;
		  scalar PrVcoeff3_;
		  scalar PrVcoeff4_;
		  scalar PrVcoeff5_;

		// Other Constants
		  scalar Trcoeff0_;
	      scalar Trcoeff1_;
	      scalar Trcoeff2_;
		  scalar Pcr_;

    // Private Member Functions

        //- Construct from components
        inline superTransport
        (
            const Thermo& t,
			const scalar muLcoeff0,
			const scalar muLcoeff1,
			const scalar muLcoeff2,
			const scalar muLcoeff3,
			const scalar muLcoeff4,
			const scalar muLcoeff5,
			const scalar muVcoeff0,
			const scalar muVcoeff1,
			const scalar muVcoeff2,
			const scalar muVcoeff3,
			const scalar muVcoeff4,
			const scalar muVcoeff5,
			const scalar PrLcoeff0,
			const scalar PrLcoeff1,
			const scalar PrLcoeff2,
			const scalar PrLcoeff3,
			const scalar PrLcoeff4,
			const scalar PrLcoeff5,
			const scalar PrVcoeff0,
			const scalar PrVcoeff1,
			const scalar PrVcoeff2,
			const scalar PrVcoeff3,
			const scalar PrVcoeff4,
			const scalar PrVcoeff5,
            const scalar Trcoeff0,
            const scalar Trcoeff1,
            const scalar Trcoeff2,
			const scalar Pcr 
        );


public:

    // Constructors

        //- Construct as named copy
        inline superTransport(const word&, const superTransport&);

        //- Construct from Istream
        superTransport(Istream&);

        //- Construct from dictionary
        superTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<superTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<superTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<superTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "superTransport<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar Pr(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline superTransport& operator=(const superTransport&);

        inline void operator+=(const superTransport&);

        inline void operator-=(const superTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend superTransport operator+ <Thermo>
        (
            const superTransport&,
            const superTransport&
        );

        friend superTransport operator- <Thermo>
        (
            const superTransport&,
            const superTransport&
        );

        friend superTransport operator* <Thermo>
        (
            const scalar,
            const superTransport&
        );

        friend superTransport operator== <Thermo>
        (
            const superTransport&,
            const superTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const superTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "superTransportI.H"

#ifdef NoRepository
#   include "superTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
