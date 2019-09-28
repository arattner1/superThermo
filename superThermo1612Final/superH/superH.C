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
    SpikeWidth1(readScalar(is)),
    cpSpike1(readScalar(is)),
    SpikeWidth2(readScalar(is)),
    cpSpike2(readScalar(is)),
    SpikeWidth3(readScalar(is)),
    cpSpike3(readScalar(is)),
    SpikeWidth4(readScalar(is)),
    cpSpike4(readScalar(is))
{
    is.check("superH::superH(Istream& is)");

    cpL  *= this->W();
    cpG  *= this->W();
    cpSpike1 *= this->W();
    cpSpike2 *= this->W();
    cpSpike3 *= this->W();
    cpSpike4 *= this->W();
}


template<class equationOfState>
Foam::superH<equationOfState>::superH(const dictionary& dict)
:
    equationOfState(dict),
    cpL(readScalar(dict.subDict("thermodynamics").lookup("cpL"))),
    cpG(readScalar(dict.subDict("thermodynamics").lookup("cpG"))),
    TSpike(readScalar(dict.subDict("thermodynamics").lookup("TSpike"))),
    SpikeWidth1(readScalar(dict.subDict("thermodynamics").lookup("SpikeWidth1"))),
    cpSpike1(readScalar(dict.subDict("thermodynamics").lookup("cpSpike1"))),
    SpikeWidth2(readScalar(dict.subDict("thermodynamics").lookup("SpikeWidth2"))),
    cpSpike2(readScalar(dict.subDict("thermodynamics").lookup("cpSpike2"))),
    SpikeWidth3(readScalar(dict.subDict("thermodynamics").lookup("SpikeWidth3"))),
    cpSpike3(readScalar(dict.subDict("thermodynamics").lookup("cpSpike3"))),
    SpikeWidth4(readScalar(dict.subDict("thermodynamics").lookup("SpikeWidth4"))),
    cpSpike4(readScalar(dict.subDict("thermodynamics").lookup("cpSpike4")))
{
    cpL  *= this->W();
    cpG  *= this->W();
    cpSpike1 *= this->W();
    cpSpike2 *= this->W();
    cpSpike3 *= this->W();
    cpSpike4 *= this->W();
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
    dict.add("SpikeWidth1", SpikeWidth1);
    dict.add("cpSpike1", cpSpike1/this->W());
    dict.add("SpikeWidth2", SpikeWidth2);
    dict.add("cpSpike2", cpSpike2/this->W());
    dict.add("SpikeWidth3", SpikeWidth3);
    dict.add("cpSpike3", cpSpike3/this->W());
    dict.add("SpikeWidth4", SpikeWidth4);
    dict.add("cpSpike4", cpSpike4/this->W());
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
        << tab << ct.SpikeWidth1
        << tab << ct.cpSpike1/ct.W()
        << tab << ct.SpikeWidth2
        << tab << ct.cpSpike2/ct.W()
        << tab << ct.SpikeWidth3
        << tab << ct.cpSpike3/ct.W()
        << tab << ct.SpikeWidth4
        << tab << ct.cpSpike4/ct.W();

    os.check("Ostream& operator<<(Ostream& os, const superH& ct)");
    return os;
}


// ************************************************************************* //
