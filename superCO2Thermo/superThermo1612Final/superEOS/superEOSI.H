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

\*---------------------------------------------------------------------------*/

#include "superEOS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::superEOS<Specie>::superEOS
(
	const Specie& sp,
	const scalar rhoL,
	const scalar rhoG,
	const scalar psiL,
	const scalar psiG,
	const scalar rhoVariationFactor,
	const scalar psiSpike,
	const scalar SpikeWidth,
	const scalar Tband,
	const scalar Pband
)
:
    Specie(sp),
    rhoL_(rhoL),
    rhoG_(rhoG),
    psiL_(psiL),
    psiG_(psiG),
    rhoVariationFactor_(rhoVariationFactor),
    psiSpike_(psiSpike),
    SpikeWidth_(SpikeWidth),
    Tband_(Tband),
    Pband_(Pband)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::superEOS<Specie>::superEOS
(
    const word& name,
    const superEOS<Specie>& pg
)
:
    Specie(name, pg),
    rhoL_(pg.rhoL_),
    rhoG_(pg.rhoG_),
    psiL_(pg.psiL_),
    psiG_(pg.psiG_),
    rhoVariationFactor_(pg.rhoVariationFactor_),
    psiSpike_(pg.psiSpike_),
    SpikeWidth_(pg.SpikeWidth_),
    Tband_(pg.Tband_),
    Pband_(pg.Pband_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::superEOS<Specie> >
Foam::superEOS<Specie>::clone() const
{
    return autoPtr<superEOS<Specie> >(new superEOS<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::superEOS<Specie> >
Foam::superEOS<Specie>::New(Istream& is)
{
    return autoPtr<superEOS<Specie> >(new superEOS<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::superEOS<Specie> > Foam::superEOS<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<superEOS<Specie> >(new superEOS<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::superEOS<Specie>::rho(scalar p, scalar T) const //rho as a function of p and T
{
	const scalar Tclamp = min(max(T,-Tband_),Tband_);
	const scalar Pclamp = min(max(p,-Pband_),Pband_);
	
	const scalar PComp = Pclamp*( 
	                              ( psiL_*Foam::exp(-Tclamp*rhoVariationFactor_/Tband_) + psiG_*Foam::exp(Tclamp*rhoVariationFactor_/Tband_) )
	                              / ( Foam::exp(-Tclamp*rhoVariationFactor_/Tband_) + Foam::exp(Tclamp*rhoVariationFactor_/Tband_) )
	                              + psiSpike_/( 1 + pow((Tclamp*SpikeWidth_/Tband_),2) ) );
	const scalar TComp = ( rhoL_*Foam::exp(-Tclamp*rhoVariationFactor_/Tband_) + rhoG_*Foam::exp(Tclamp*rhoVariationFactor_/Tband_) )
	                     / ( Foam::exp(-Tclamp*rhoVariationFactor_/Tband_) + Foam::exp(Tclamp*rhoVariationFactor_/Tband_) );
	                              
	
	
	return TComp + PComp;
}


template<class Specie>
inline Foam::scalar Foam::superEOS<Specie>::psi(scalar, scalar T) const //d(rho) / d(p)
{
	const scalar Tclamp = min(max(T,-Tband_),Tband_);
	//const scalar Pclamp = min(max(p,-Pband_),Pband_);

	const scalar psi    = ( psiL_*Foam::exp(-Tclamp*rhoVariationFactor_/Tband_) + psiG_*Foam::exp(Tclamp*rhoVariationFactor_/Tband_) )
	                        / ( Foam::exp(-Tclamp*rhoVariationFactor_/Tband_) + Foam::exp(Tclamp*rhoVariationFactor_/Tband_) )
	                        + psiSpike_/( 1 + pow((Tclamp*SpikeWidth_/Tband_),2) );


    return psi;
}


template<class Specie>
inline Foam::scalar Foam::superEOS<Specie>::Z(scalar, scalar) const //I think this is: p/(rho*R*T)
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::superEOS<Specie>::cpMcv(scalar, scalar) const //This is c_p - c_v
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::superEOS<Specie>::operator+=(const superEOS<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::superEOS<Specie>::operator-=(const superEOS<Specie>& pg)
{
    Specie::operator-=(pg);
}


template<class Specie>
inline void Foam::superEOS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::superEOS<Specie> Foam::operator+
(
    const superEOS<Specie>& pg1,
    const superEOS<Specie>& pg2
)
{
    return superEOS<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::superEOS<Specie> Foam::operator-
(
    const superEOS<Specie>& pg1,
    const superEOS<Specie>& pg2
)
{
    return superEOS<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::superEOS<Specie> Foam::operator*
(
    const scalar s,
    const superEOS<Specie>& pg
)
{
    return superEOS<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::superEOS<Specie> Foam::operator==
(
    const superEOS<Specie>& pg1,
    const superEOS<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
