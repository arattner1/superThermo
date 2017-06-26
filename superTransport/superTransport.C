/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "superTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::superTransport<Thermo>::superTransport(Istream& is)
:
    Thermo(is),
    muLcoeff0_(readScalar(is)),
    muLcoeff1_(readScalar(is)),
    muLcoeff2_(readScalar(is)),
    muLcoeff3_(readScalar(is)),
    muLcoeff4_(readScalar(is)),
    muLcoeff5_(readScalar(is)),
    muVcoeff0_(readScalar(is)),
    muVcoeff1_(readScalar(is)),
    muVcoeff2_(readScalar(is)),
    muVcoeff3_(readScalar(is)),
    muVcoeff4_(readScalar(is)),
    muVcoeff5_(readScalar(is)),
    PrLcoeff0_(readScalar(is)),
    PrLcoeff1_(readScalar(is)),
    PrLcoeff2_(readScalar(is)),
    PrLcoeff3_(readScalar(is)),
    PrLcoeff4_(readScalar(is)),
    PrLcoeff5_(readScalar(is)),
    PrVcoeff0_(readScalar(is)),
    PrVcoeff1_(readScalar(is)),
    PrVcoeff2_(readScalar(is)),
    PrVcoeff3_(readScalar(is)),
    PrVcoeff4_(readScalar(is)),
    PrVcoeff5_(readScalar(is)),
	Trcoeff0_(readScalar(is)),
	Trcoeff1_(readScalar(is)),
	Trcoeff2_(readScalar(is)),
	Pcr_(readScalar(is))

{
    is.check("superTransport::superTransport(Istream& is)");
}


template<class Thermo>
Foam::superTransport<Thermo>::superTransport(const dictionary& dict)
:
    Thermo(dict),
    muLcoeff0_(readScalar(dict.subDict("transport").lookup("muLcoeff0"))),
    muLcoeff1_(readScalar(dict.subDict("transport").lookup("muLcoeff1"))),
    muLcoeff2_(readScalar(dict.subDict("transport").lookup("muLcoeff2"))),
    muLcoeff3_(readScalar(dict.subDict("transport").lookup("muLcoeff3"))),
    muLcoeff4_(readScalar(dict.subDict("transport").lookup("muLcoeff4"))),
    muLcoeff5_(readScalar(dict.subDict("transport").lookup("muLcoeff5"))),
    muVcoeff0_(readScalar(dict.subDict("transport").lookup("muVcoeff0"))),
    muVcoeff1_(readScalar(dict.subDict("transport").lookup("muVcoeff1"))),
    muVcoeff2_(readScalar(dict.subDict("transport").lookup("muVcoeff2"))),
    muVcoeff3_(readScalar(dict.subDict("transport").lookup("muVcoeff3"))),
    muVcoeff4_(readScalar(dict.subDict("transport").lookup("muVcoeff4"))),
    muVcoeff5_(readScalar(dict.subDict("transport").lookup("muVcoeff5"))),
    PrLcoeff0_(readScalar(dict.subDict("transport").lookup("PrLcoeff0"))),
    PrLcoeff1_(readScalar(dict.subDict("transport").lookup("PrLcoeff1"))),
    PrLcoeff2_(readScalar(dict.subDict("transport").lookup("PrLcoeff2"))),
    PrLcoeff3_(readScalar(dict.subDict("transport").lookup("PrLcoeff3"))),
    PrLcoeff4_(readScalar(dict.subDict("transport").lookup("PrLcoeff4"))),
    PrLcoeff5_(readScalar(dict.subDict("transport").lookup("PrLcoeff5"))),
    PrVcoeff0_(readScalar(dict.subDict("transport").lookup("PrVcoeff0"))),
    PrVcoeff1_(readScalar(dict.subDict("transport").lookup("PrVcoeff1"))),
    PrVcoeff2_(readScalar(dict.subDict("transport").lookup("PrVcoeff2"))),
    PrVcoeff3_(readScalar(dict.subDict("transport").lookup("PrVcoeff3"))),
    PrVcoeff4_(readScalar(dict.subDict("transport").lookup("PrVcoeff4"))),
    PrVcoeff5_(readScalar(dict.subDict("transport").lookup("PrVcoeff5"))),
	Trcoeff0_(readScalar(dict.subDict("transport").lookup("Trcoeff0"))),
	Trcoeff1_(readScalar(dict.subDict("transport").lookup("Trcoeff1"))),
	Trcoeff2_(readScalar(dict.subDict("transport").lookup("Trcoeff2"))),
    Pcr_(readScalar(dict.subDict("transport").lookup("Pcr")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::superTransport<Thermo>::superTransport::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("muLcoeff0", muLcoeff0_);
    dict.add("muLcoeff1", muLcoeff1_);
    dict.add("muLcoeff2", muLcoeff2_);
    dict.add("muLcoeff3", muLcoeff3_);
    dict.add("muLcoeff4", muLcoeff4_);
    dict.add("muLcoeff5", muLcoeff5_);
    dict.add("muVcoeff0", muVcoeff0_);
    dict.add("muVcoeff1", muVcoeff1_);
    dict.add("muVcoeff2", muVcoeff2_);
    dict.add("muVcoeff3", muVcoeff3_);
    dict.add("muVcoeff4", muVcoeff4_);
    dict.add("muVcoeff5", muVcoeff5_);
    dict.add("PrLcoeff0", PrLcoeff0_);
    dict.add("PrLcoeff1", PrLcoeff1_);
    dict.add("PrLcoeff2", PrLcoeff2_);
    dict.add("PrLcoeff3", PrLcoeff3_);
    dict.add("PrLcoeff4", PrLcoeff4_);
    dict.add("PrLcoeff5", PrLcoeff5_);
    dict.add("PrVcoeff0", PrVcoeff0_);
    dict.add("PrVcoeff1", PrVcoeff1_);
    dict.add("PrVcoeff2", PrVcoeff2_);
    dict.add("PrVcoeff3", PrVcoeff3_);
    dict.add("PrVcoeff4", PrVcoeff4_);
    dict.add("PrVcoeff5", PrVcoeff5_);
	dict.add("Trcoeff0", Trcoeff0_);
	dict.add("Trcoeff1", Trcoeff1_);
	dict.add("Trcoeff2", Trcoeff2_);
    dict.add("Pcr", Pcr_);

    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const superTransport<Thermo>& ct)
{
    operator<<(os, static_cast<const Thermo&>(ct));
    os << tab << ct.muLcoeff0_;
    os << tab << ct.muLcoeff1_;
    os << tab << ct.muLcoeff2_;
    os << tab << ct.muLcoeff3_;
    os << tab << ct.muLcoeff4_;
    os << tab << ct.muLcoeff5_;
    os << tab << ct.muVcoeff0_;
    os << tab << ct.muVcoeff1_;
    os << tab << ct.muVcoeff2_;
    os << tab << ct.muVcoeff3_;
    os << tab << ct.muVcoeff4_;
    os << tab << ct.muVcoeff5_;
    os << tab << ct.PrLcoeff0_;
    os << tab << ct.PrLcoeff1_;
    os << tab << ct.PrLcoeff2_;
    os << tab << ct.PrLcoeff3_;
    os << tab << ct.PrLcoeff4_;
    os << tab << ct.PrLcoeff5_;
    os << tab << ct.PrVcoeff0_;
    os << tab << ct.PrVcoeff1_;
    os << tab << ct.PrVcoeff2_;
    os << tab << ct.PrVcoeff3_;
    os << tab << ct.PrVcoeff4_;
    os << tab << ct.PrVcoeff5_;
	os << tab << ct.Trcoeff0_;
	os << tab << ct.Trcoeff1_;
    os << tab << ct.Trcoeff2_;
	os << tab << ct.Pcr_;


    os.check("Ostream& operator<<(Ostream&, const superTransport&)");

    return os;
}


// ************************************************************************* //
