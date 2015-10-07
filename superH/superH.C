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

#include "superH.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::superH<equationOfState>::superH(Istream& is)
:
    equationOfState(is),
    cpL_(readScalar(is)),
    cpG_(readScalar(is)),
    cpVF_(readScalar(is)),
    TSpike_(readScalar(is)),
    SpikeWidth_(readScalar(is)),
    Tband_(readScalar(is))
{
    is.check("superH::superH(Istream& is)");

    cpL_ *= this->W();
    cpG_ *= this->W();
    cpVF_ *= this->W();
    TSpike_ *= this->W();
    SpikeWidth_ *= this->W();
    Tband_ *= this->W();
}


template<class equationOfState>
Foam::superH<equationOfState>::superH(const dictionary& dict)
:
    equationOfState(dict),
    cpL_(readScalar(dict.subDict("thermodynamics").lookup("cpL"))),
    cpG_(readScalar(dict.subDict("thermodynamics").lookup("cpG"))),
    cpVF_(readScalar(dict.subDict("thermodynamics").lookup("cpVariationFactor"))),
    TSpike_(readScalar(dict.subDict("thermodynamics").lookup("cpSpike"))),
    SpikeWidth_(readScalar(dict.subDict("thermodynamics").lookup("SpikeWidth"))),
    Tband_(readScalar(dict.subDict("thermodynamics").lookup("Tband")))
{
    cpL_ *= this->W();
    cpG_ *= this->W();
    cpVF_ *= this->W();
    TSpike_ *= this->W();
    SpikeWidth_ *= this->W();
    Tband_ *= this->W();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
void Foam::superH<equationOfState>::write(Ostream& os) const
{
    equationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("cpL", cpL_/this->W());
    dict.add("cpG", cpG_/this->W());
    dict.add("cpVariationFactor", cpVF_/this->W());
    dict.add("cpSpike", TSpike_/this->W());
    dict.add("SpikeWidth", SpikeWidth_/this->W());
    dict.add("Tband", Tband_/this->W());
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const superH<equationOfState>& ct
)
{
    os  << static_cast<const equationOfState&>(ct)
        << tab << ct.cpL_/ct.W()
        << tab << ct.cpG_/ct.W()
        << tab << ct.cpVF_/ct.W()
        << tab << ct.TSpike_/ct.W()
        << tab << ct.SpikeWidth_/ct.W()
        << tab << ct.Tband_/ct.W();

    os.check("Ostream& operator<<(Ostream& os, const superH& ct)");
    return os;
}


// ************************************************************************* //
