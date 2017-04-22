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
    cpL(readScalar(is)),
    cpG(readScalar(is)),
    TSpike(readScalar(is)),
    SpikeWidth(readScalar(is)),
    cpSpike(readScalar(is))
{
    is.check("superH::superH(Istream& is)");

    cpL  *= this->W();
    cpG  *= this->W();
    cpSpike *= this->W();
}


template<class equationOfState>
Foam::superH<equationOfState>::superH(const dictionary& dict)
:
    equationOfState(dict),
    cpL(readScalar(dict.subDict("thermodynamics").lookup("cpL"))),
    cpG(readScalar(dict.subDict("thermodynamics").lookup("cpG"))),
    TSpike(readScalar(dict.subDict("thermodynamics").lookup("TSpike"))),
    SpikeWidth(readScalar(dict.subDict("thermodynamics").lookup("SpikeWidth"))),
    cpSpike(readScalar(dict.subDict("thermodynamics").lookup("cpSpike")))
{
    cpL  *= this->W();
    cpG  *= this->W();
    cpSpike *= this->W();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
void Foam::superH<equationOfState>::write(Ostream& os) const
{
    equationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("cpL", cpL/this->W());
    dict.add("cpG", cpG/this->W());
    dict.add("TSpike", TSpike);
    dict.add("SpikeWidth", SpikeWidth);
    dict.add("cpSpike", cpSpike/this->W());
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
        << tab << ct.cpL/ct.W()
        << tab << ct.cpG/ct.W()
        << tab << ct.TSpike
        << tab << ct.SpikeWidth
        << tab << ct.cpSpike/ct.W();

    os.check("Ostream& operator<<(Ostream& os, const superH& ct)");
    return os;
}


// ************************************************************************* //
