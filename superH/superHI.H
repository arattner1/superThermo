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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
inline Foam::superH<equationOfState>::superH
(
    const equationOfState& st,
    const scalar cpL,
    const scalar cpG,
    const scalar cpVF,
    const scalar TSpike,
    const scalar SpikeWidth,
    const scalar Tband
)
:
    equationOfState(st),
    cpL_(cpL),
    cpG_(cpG),
    cpVF_(cpVF),
    TSpike_(TSpike),
    SpikeWidth_(SpikeWidth),
    Tband_(Tband)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::superH<equationOfState>::superH
(
    const word& name,
    const superH& ct
)
:
    equationOfState(name, ct),
    cpL_(ct.cpL),
    cpG_(ct.cpG),
    cpVF_(ct.cpVF),
    TSpike_(ct.TSpike),
    SpikeWidth_(ct.SpikeWidth),
    Tband_(ct.Tband)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::superH<equationOfState> >
Foam::superH<equationOfState>::clone() const
{
    return autoPtr<superH<equationOfState> >
    (
        new superH<equationOfState>(*this)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::superH<equationOfState> >
Foam::superH<equationOfState>::New(Istream& is)
{
    return autoPtr<superH<equationOfState> >
    (
        new superH<equationOfState>(is)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::superH<equationOfState> >
Foam::superH<equationOfState>::New(const dictionary& dict)
{
    return autoPtr<superH<equationOfState> >
    (
        new superH<equationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::superH<EquationOfState>::limit
(
    const scalar T
) const
{
    	return T;
	//return max(min(Tband_,T),-Tband_);
}


template<class equationOfState>
inline Foam::scalar Foam::superH<equationOfState>::cp
(
    const scalar p,
    const scalar T
) const
{
	const scalar cp = ( cpL_*Foam::exp(-cpVF_*T/Tband_) + cpG_*Foam::exp(cpVF_*T/Tband_) )
	                  / ( Foam::exp(-cpVF_*T/Tband_) + Foam::exp(cpVF_*T/Tband_) ) + TSpike_/(1 + pow((SpikeWidth_*T/Tband_),2) );

    return cp;
}


template<class equationOfState>
inline Foam::scalar Foam::superH<equationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    const scalar c = cpVF_/Tband_;
    const scalar B = SpikeWidth_/Tband_;
    const scalar H_main = cpL_*T + (cpL_ - cpG_)*Foam::log( Foam::exp(2*c*T) + 1 )/(2*c);
    const scalar H_spike = TSpike_*Foam::atan(B*T)/B;
    return (H_main + H_spike);
}


template<class equationOfState>
inline Foam::scalar Foam::superH<equationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    return ha(p, T);
}


template<class equationOfState>
inline Foam::scalar Foam::superH<equationOfState>::hc() const
{
    return 0.0;
}


template<class equationOfState>
inline Foam::scalar Foam::superH<equationOfState>::s
(
    const scalar p, const scalar T
) const
{
    notImplemented
    (
        "scalar superH<equationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    return T;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//Never tested...
template<class equationOfState>
inline void Foam::superH<equationOfState>::operator+=
(
    const superH<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    cpL_    = molr1*cpL_ + molr2*ct.cpL_;
    cpG_    = molr1*cpG_ + molr2*ct.cpG_;
    cpVF_   = molr1*cpVF_ + molr2*ct.cpVF_;
    TSpike_ = molr1*TSpike_ + molr2*ct.TSpike_;
    SpikeWidth_ = molr1*SpikeWidth_ + molr2*ct.SpikeWidth_;
    Tband_ = molr1*Tband_ + molr2*ct.Tband_;
}


template<class equationOfState>
inline void Foam::superH<equationOfState>::operator-=
(
    const superH<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    cpL_    = molr1*cpL_ - molr2*ct.cpL_;
    cpG_    = molr1*cpG_ - molr2*ct.cpG_;
    cpVF_   = molr1*cpVF_ - molr2*ct.cpVF_;
    TSpike_ = molr1*TSpike_ - molr2*ct.TSpike_;
    SpikeWidth_ = molr1*SpikeWidth_ - molr2*ct.SpikeWidth_;
    Tband_ = molr1*Tband_ - molr2*ct.Tband_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::superH<equationOfState> Foam::operator+
(
    const superH<equationOfState>& ct1,
    const superH<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    return superH<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.CpL_
      + ct2.nMoles()/eofs.nMoles()*ct2.CpL_,
        ct1.nMoles()/eofs.nMoles()*ct1.CpG_
      + ct2.nMoles()/eofs.nMoles()*ct2.CpG_,
        ct1.nMoles()/eofs.nMoles()*ct1.CpVF_
      + ct2.nMoles()/eofs.nMoles()*ct2.CpVF_,
        ct1.nMoles()/eofs.nMoles()*ct1.TSpike_
      + ct2.nMoles()/eofs.nMoles()*ct2.TSpike_,
        ct1.nMoles()/eofs.nMoles()*ct1.SpikeWidth_
      + ct2.nMoles()/eofs.nMoles()*ct2.SpikeWidth_,
        ct1.nMoles()/eofs.nMoles()*ct1.Tband_
      + ct2.nMoles()/eofs.nMoles()*ct2.Tband_
    );
}


template<class equationOfState>
inline Foam::superH<equationOfState> Foam::operator-
(
    const superH<equationOfState>& ct1,
    const superH<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      - static_cast<const equationOfState&>(ct2)
    );

    return superH<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.CpL_
      - ct2.nMoles()/eofs.nMoles()*ct2.CpL_,
        ct1.nMoles()/eofs.nMoles()*ct1.CpG_
      - ct2.nMoles()/eofs.nMoles()*ct2.CpG_,
        ct1.nMoles()/eofs.nMoles()*ct1.CpVF_
      - ct2.nMoles()/eofs.nMoles()*ct2.CpVF_,
        ct1.nMoles()/eofs.nMoles()*ct1.TSpike_
      - ct2.nMoles()/eofs.nMoles()*ct2.TSpike_,
        ct1.nMoles()/eofs.nMoles()*ct1.SpikeWidth_
      - ct2.nMoles()/eofs.nMoles()*ct2.SpikeWidth_,
        ct1.nMoles()/eofs.nMoles()*ct1.Tband_
      - ct2.nMoles()/eofs.nMoles()*ct2.Tband_
    );
}


template<class equationOfState>
inline Foam::superH<equationOfState> Foam::operator*
(
    const scalar s,
    const superH<equationOfState>& ct
)
{
    return superH<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.cpL_,
        ct.cpG_,
        ct.cpVF_,
        ct.TSpike_,
        ct.SpikeWidth_,
        ct.Tband_
    );
}


template<class equationOfState>
inline Foam::superH<equationOfState> Foam::operator==
(
    const superH<equationOfState>& ct1,
    const superH<equationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
