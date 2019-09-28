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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::superEOS<Specie>::superEOS(Istream& is)
:
    Specie(is),
	rhoL_(readScalar(is)),
	rhoG_(readScalar(is)),
	psiL_(readScalar(is)),
	psiG_(readScalar(is)),
	rhoVariationFactor_(readScalar(is)),
	psiSpike_(readScalar(is)),
	SpikeWidth_(readScalar(is)),
	Tband_(readScalar(is)),
	Pband_(readScalar(is))
{
    is.check("superEOS<Specie>::superEOS(Istream& is)");
}


template<class Specie>
Foam::superEOS<Specie>::superEOS(const dictionary& dict)
:
    Specie(dict),
    rhoL_(readScalar(dict.subDict("equationOfState").lookup("rhoL"))),
    rhoG_(readScalar(dict.subDict("equationOfState").lookup("rhoG"))),
    psiL_(readScalar(dict.subDict("equationOfState").lookup("psiL"))),
    psiG_(readScalar(dict.subDict("equationOfState").lookup("psiG"))),
    rhoVariationFactor_(readScalar(dict.subDict("equationOfState").lookup("rhoVariationFactor"))),
    psiSpike_(readScalar(dict.subDict("equationOfState").lookup("psiSpike"))),
    SpikeWidth_(readScalar(dict.subDict("equationOfState").lookup("SpikeWidth"))),
    Tband_(readScalar(dict.subDict("equationOfState").lookup("Tband"))),
    Pband_(readScalar(dict.subDict("equationOfState").lookup("Pband")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::superEOS<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    
    dictionary dict("equationOfState");
    dict.add("rhoL", rhoL_);
    dict.add("rhoG", rhoG_);
    dict.add("psiL", psiL_);
    dict.add("psiG", psiG_);
    dict.add("rhoVariationFactor", rhoVariationFactor_);
    dict.add("psiSpike", psiSpike_);
    dict.add("SpikeWidth", SpikeWidth_);
    dict.add("Tband", Tband_);
    dict.add("Pband", Pband_);

    
    os  << indent << dict.dictName() << dict;
    
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const superEOS<Specie>& pg)
{
    os  << static_cast<const Specie&>(pg)
        << token::SPACE << pg.rhoL_
        << token::SPACE << pg.rhoG_
        << token::SPACE << pg.psiL_
        << token::SPACE << pg.psiG_
		<< token::SPACE << pg.rhoVariationFactor_
		<< token::SPACE << pg.psiSpike_
		<< token::SPACE << pg.SpikeWidth_
		<< token::SPACE << pg.Tband_
        << token::SPACE << pg.Pband_;
        
    os.check("Ostream& operator<<(Ostream& os, const superEOS<Specie>& st)");
    return os;
}


// ************************************************************************* //
