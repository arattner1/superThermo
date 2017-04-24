/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::superTransport<Thermo>::superTransport
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
)
:
    Thermo(t),
    muLcoeff0_(muLcoeff0),
    muLcoeff1_(muLcoeff1),
    muLcoeff2_(muLcoeff2),
    muLcoeff3_(muLcoeff3),
    muLcoeff4_(muLcoeff4),
    muLcoeff5_(muLcoeff5),
    muVcoeff0_(muVcoeff0),
    muVcoeff1_(muVcoeff1),
    muVcoeff2_(muVcoeff2),
    muVcoeff3_(muVcoeff3),
    muVcoeff4_(muVcoeff4),
    muVcoeff5_(muVcoeff5),
    PrLcoeff0_(PrLcoeff0),
    PrLcoeff1_(PrLcoeff1),
    PrLcoeff2_(PrLcoeff2),
    PrLcoeff3_(PrLcoeff3),
    PrLcoeff4_(PrLcoeff4),
    PrLcoeff5_(PrLcoeff5),
    PrVcoeff0_(PrVcoeff0),
    PrVcoeff1_(PrVcoeff1),
    PrVcoeff2_(PrVcoeff2),
    PrVcoeff3_(PrVcoeff3),
    PrVcoeff4_(PrVcoeff4),
    PrVcoeff5_(PrVcoeff5),
	Trcoeff0_(Trcoeff0),
	Trcoeff1_(Trcoeff1),
	Trcoeff2_(Trcoeff2),
	Pcr_(Pcr)
{}


template<class Thermo>
inline Foam::superTransport<Thermo>::superTransport
(
    const word& name,
    const superTransport& ct
)
:
    Thermo(name, ct),
    muLcoeff0_(ct.muLcoeff0_),
    muLcoeff1_(ct.muLcoeff1_),
    muLcoeff2_(ct.muLcoeff2_),
    muLcoeff3_(ct.muLcoeff3_),
    muLcoeff4_(ct.muLcoeff4_),
    muLcoeff5_(ct.muLcoeff5_),
    muVcoeff0_(ct.muVcoeff0_),
    muVcoeff1_(ct.muVcoeff1_),
    muVcoeff2_(ct.muVcoeff2_),
    muVcoeff3_(ct.muVcoeff3_),
    muVcoeff4_(ct.muVcoeff4_),
    muVcoeff5_(ct.muVcoeff5_),
    PrLcoeff0_(ct.PrLcoeff0_),
    PrLcoeff1_(ct.PrLcoeff1_),
    PrLcoeff2_(ct.PrLcoeff2_),
    PrLcoeff3_(ct.PrLcoeff3_),
    PrLcoeff4_(ct.PrLcoeff4_),
    PrLcoeff5_(ct.PrLcoeff5_),
    PrVcoeff0_(ct.PrVcoeff0_),
    PrVcoeff1_(ct.PrVcoeff1_),
    PrVcoeff2_(ct.PrVcoeff2_),
    PrVcoeff3_(ct.PrVcoeff3_),
    PrVcoeff4_(ct.PrVcoeff4_),
    PrVcoeff5_(ct.PrVcoeff5_),
	Trcoeff0_(ct.Trcoeff0_),
	Trcoeff1_(ct.Trcoeff1_),
	Trcoeff2_(ct.Trcoeff2_),
	Pcr_(ct.Pcr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::superTransport<Thermo> >
Foam::superTransport<Thermo>::clone() const
{
    return autoPtr<superTransport<Thermo> >
    (
        new superTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::superTransport<Thermo> >
Foam::superTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<superTransport<Thermo> >
    (
        new superTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::superTransport<Thermo> >
Foam::superTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<superTransport<Thermo> >
    (
        new superTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::superTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
 	const scalar Tr = Trcoeff0_+Trcoeff1_*T+Trcoeff2_*p;
	const scalar pr = p/Pcr_;

	if ( (Tr <= 1) && (Tr >= -1))
		{
			const scalar f=(Tr+1.0)/2.0;
			const scalar muL = muLcoeff0_ + muLcoeff1_*Tr + muLcoeff2_*pow(Tr,2) + muLcoeff3_*pow(Tr,3) + muLcoeff4_/(pr+muLcoeff5_);
			const scalar muV = muVcoeff0_ + muVcoeff1_*Tr + muVcoeff2_*pow(Tr,2) + muVcoeff3_*pow(Tr,3) + muVcoeff4_/(pr+muVcoeff5_);
		    const scalar mu = f*muV + (1.0-f)*muL;
			Info<< " first condition for mu" << endl;
			Info<< "mu = " << mu << endl;
    		return mu; 
		}

	else if (Tr > 1)
		{
			const scalar mu = muVcoeff0_ + muVcoeff1_*Tr + muVcoeff2_*pow(Tr,2) + muVcoeff3_*pow(Tr,3) + muVcoeff4_/(pr+muVcoeff5_);
		//	Info<< " second condition for mu" << endl;
		//	Info<< "mu = " << mu << endl;
    		return mu;  
		}

	else
		{
			const scalar mu = muLcoeff0_ + muLcoeff1_*Tr + muLcoeff2_*pow(Tr,2) + muLcoeff3_*pow(Tr,3) + muLcoeff4_/(pr+muLcoeff5_);
		//	Info<< " third condition for mu" << endl;
		//	Info<< "mu = " << mu << endl;
    		return mu;  
		} 
}

template<class Thermo>
inline Foam::scalar Foam::superTransport<Thermo>::Pr
(
    const scalar p,
    const scalar T
) const
{
 	const scalar Tr = Trcoeff0_+Trcoeff1_*T+Trcoeff2_*p;
	const scalar pr = p/Pcr_;

	if ( (Tr <= 1) && (Tr >= -1))
		{
			const scalar f=(Tr+1.0)/2.0;
			const scalar PrL = PrLcoeff0_*Foam::exp(PrLcoeff1_*Tr) + PrLcoeff2_*Foam::exp(PrLcoeff3_*Tr) + PrLcoeff4_/(pr+PrLcoeff5_);
			const scalar PrV = PrVcoeff0_*Foam::exp(PrVcoeff1_*Tr) + PrVcoeff2_*Foam::exp(PrVcoeff3_*Tr) + PrVcoeff4_/(pr+PrVcoeff5_);
		    const scalar Pr = f*PrV + (1.0-f)*PrL;
			Info<< " first condition for Pr" << endl;
			Info<< "Pr = " << Pr << endl;
    		return Pr; 
		}

	else if (Tr > 1)
		{
			const scalar Pr = PrVcoeff0_*Foam::exp(PrVcoeff1_*Tr) + PrVcoeff2_*Foam::exp(PrVcoeff3_*Tr) + PrVcoeff4_/(pr+PrVcoeff5_);
		//	Info<< " second condition for Pr" << endl;
		//	Info<< "Pr = " << Pr << endl;
    		return Pr;  
		}

	else
		{
			const scalar Pr = PrLcoeff0_*Foam::exp(PrLcoeff1_*Tr) + PrLcoeff2_*Foam::exp(PrLcoeff3_*Tr) + PrLcoeff4_/(pr+PrLcoeff5_);
		//	Info<< " third condition for Pr" << endl;
		//	Info<< "Pr = " << Pr << endl;
    		return Pr;  
		} 
}
//Thermal Conductivity
template<class Thermo>
inline Foam::scalar Foam::superTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
	//Printing Cp values for checking
    Info<< "from superTransport library Cp =  " << (this->Cp)(p, T) << endl;
    return this->Cp(p, T)*mu(p, T)*(1.0/Pr(p, T));
}


template<class Thermo>
inline Foam::scalar Foam::superTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*(1.0/Pr(p, T));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::superTransport<Thermo>& Foam::superTransport<Thermo>::operator=
(
    const superTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

	  muLcoeff0_ = ct.muLcoeff0_;
	  muLcoeff1_ = ct.muLcoeff1_;
	  muLcoeff2_ = ct.muLcoeff2_;
	  muLcoeff3_ = ct.muLcoeff3_;
	  muLcoeff4_ = ct.muLcoeff4_;
	  muLcoeff5_ = ct.muLcoeff5_;
	  muVcoeff0_ = ct.muVcoeff0_;
	  muVcoeff1_ = ct.muVcoeff1_;
	  muVcoeff2_ = ct.muVcoeff2_;
	  muVcoeff3_ = ct.muVcoeff3_;
	  muVcoeff4_ = ct.muVcoeff4_;
	  muVcoeff5_ = ct.muVcoeff5_;
	  PrLcoeff0_ = ct.PrLcoeff0_;
	  PrLcoeff1_ = ct.PrLcoeff1_;
	  PrLcoeff2_ = ct.PrLcoeff2_;
	  PrLcoeff3_ = ct.PrLcoeff3_;
	  PrLcoeff4_ = ct.PrLcoeff4_;
	  PrLcoeff5_ = ct.PrLcoeff5_;
	  PrVcoeff0_ = ct.PrVcoeff0_;
	  PrVcoeff1_ = ct.PrVcoeff1_;
	  PrVcoeff2_ = ct.PrVcoeff2_;
	  PrVcoeff3_ = ct.PrVcoeff3_;
	  PrVcoeff4_ = ct.PrVcoeff4_;
	  PrVcoeff5_ = ct.PrVcoeff5_;
	  Trcoeff0_ = ct.Trcoeff0_;
	  Trcoeff1_ = ct.Trcoeff1_;
	  Trcoeff2_ = ct.Trcoeff2_;
	  Pcr_ = ct.Pcr_;

    return *this;
}


template<class Thermo>
inline void Foam::superTransport<Thermo>::operator+=
(
    const superTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

		muLcoeff0_    = molr1*muLcoeff0_ + molr2*st.muLcoeff0_;
		muLcoeff1_    = molr1*muLcoeff1_ + molr2*st.muLcoeff1_;
		muLcoeff2_    = molr1*muLcoeff2_ + molr2*st.muLcoeff2_;
		muLcoeff3_    = molr1*muLcoeff3_ + molr2*st.muLcoeff3_;
		muLcoeff4_    = molr1*muLcoeff4_ + molr2*st.muLcoeff4_;
		muLcoeff5_    = molr1*muLcoeff5_ + molr2*st.muLcoeff5_;
		muVcoeff0_    = molr1*muVcoeff0_ + molr2*st.muVcoeff0_;
		muVcoeff1_    = molr1*muVcoeff1_ + molr2*st.muVcoeff1_;
		muVcoeff2_    = molr1*muVcoeff2_ + molr2*st.muVcoeff2_;
		muVcoeff3_    = molr1*muVcoeff3_ + molr2*st.muVcoeff3_;
		muVcoeff4_    = molr1*muVcoeff4_ + molr2*st.muVcoeff4_;
		muVcoeff5_    = molr1*muVcoeff5_ + molr2*st.muVcoeff5_;
		PrLcoeff0_    = molr1*PrLcoeff0_ + molr2*st.PrLcoeff0_;
		PrLcoeff1_    = molr1*PrLcoeff1_ + molr2*st.PrLcoeff1_;
		PrLcoeff2_    = molr1*PrLcoeff2_ + molr2*st.PrLcoeff2_;
		PrLcoeff3_    = molr1*PrLcoeff3_ + molr2*st.PrLcoeff3_;
		PrLcoeff4_    = molr1*PrLcoeff4_ + molr2*st.PrLcoeff4_;
		PrLcoeff5_    = molr1*PrLcoeff5_ + molr2*st.PrLcoeff5_;
		PrVcoeff0_    = molr1*PrVcoeff0_ + molr2*st.PrVcoeff0_;
		PrVcoeff1_    = molr1*PrVcoeff1_ + molr2*st.PrVcoeff1_;
		PrVcoeff2_    = molr1*PrVcoeff2_ + molr2*st.PrVcoeff2_;
		PrVcoeff3_    = molr1*PrVcoeff3_ + molr2*st.PrVcoeff3_;
		PrVcoeff4_    = molr1*PrVcoeff4_ + molr2*st.PrVcoeff4_;
		PrVcoeff5_    = molr1*PrVcoeff5_ + molr2*st.PrVcoeff5_;
		Trcoeff0_    = molr1*Trcoeff0_ + molr2*st.Trcoeff0_;
		Trcoeff1_    = molr1*Trcoeff1_ + molr2*st.Trcoeff1_;
		Trcoeff2_    = molr1*Trcoeff2_ + molr2*st.Trcoeff2_;
		Pcr_    = molr1*Pcr_ + molr2*st.Pcr_;
    }
}


template<class Thermo>
inline void Foam::superTransport<Thermo>::operator-=
(
    const superTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

		muLcoeff0_    = molr1*muLcoeff0_ - molr2*st.muLcoeff0_;
		muLcoeff1_    = molr1*muLcoeff1_ - molr2*st.muLcoeff1_;
		muLcoeff2_    = molr1*muLcoeff2_ - molr2*st.muLcoeff2_;
		muLcoeff3_    = molr1*muLcoeff3_ - molr2*st.muLcoeff3_;
		muLcoeff4_    = molr1*muLcoeff4_ - molr2*st.muLcoeff4_;
		muLcoeff5_    = molr1*muLcoeff5_ - molr2*st.muLcoeff5_;
		muVcoeff0_    = molr1*muVcoeff0_ - molr2*st.muVcoeff0_;
		muVcoeff1_    = molr1*muVcoeff1_ - molr2*st.muVcoeff1_;
		muVcoeff2_    = molr1*muVcoeff2_ - molr2*st.muVcoeff2_;
		muVcoeff3_    = molr1*muVcoeff3_ - molr2*st.muVcoeff3_;
		muVcoeff4_    = molr1*muVcoeff4_ - molr2*st.muVcoeff4_;
		muVcoeff5_    = molr1*muVcoeff5_ - molr2*st.muVcoeff5_;
		PrLcoeff0_    = molr1*PrLcoeff0_ - molr2*st.PrLcoeff0_;
		PrLcoeff1_    = molr1*PrLcoeff1_ - molr2*st.PrLcoeff1_;
		PrLcoeff2_    = molr1*PrLcoeff2_ - molr2*st.PrLcoeff2_;
		PrLcoeff3_    = molr1*PrLcoeff3_ - molr2*st.PrLcoeff3_;
		PrLcoeff4_    = molr1*PrLcoeff4_ - molr2*st.PrLcoeff4_;
		PrLcoeff5_    = molr1*PrLcoeff5_ - molr2*st.PrLcoeff5_;
		PrVcoeff0_    = molr1*PrVcoeff0_ - molr2*st.PrVcoeff0_;
		PrVcoeff1_    = molr1*PrVcoeff1_ - molr2*st.PrVcoeff1_;
		PrVcoeff2_    = molr1*PrVcoeff2_ - molr2*st.PrVcoeff2_;
		PrVcoeff3_    = molr1*PrVcoeff3_ - molr2*st.PrVcoeff3_;
		PrVcoeff4_    = molr1*PrVcoeff4_ - molr2*st.PrVcoeff4_;
		PrVcoeff5_    = molr1*PrVcoeff5_ - molr2*st.PrVcoeff5_;
		Trcoeff0_    = molr1*Trcoeff0_ - molr2*st.Trcoeff0_;
		Trcoeff1_    = molr1*Trcoeff1_ - molr2*st.Trcoeff1_;
		Trcoeff2_    = molr1*Trcoeff2_ - molr2*st.Trcoeff2_;
		Pcr_    = molr1*Pcr_ - molr2*st.Pcr_;
    }
}


template<class Thermo>
inline void Foam::superTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::superTransport<Thermo> Foam::operator+
(
    const superTransport<Thermo>& ct1,
    const superTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(ct1.nMoles()) + mag(ct2.nMoles()) < SMALL)
    {
        return superTransport<Thermo>
        (
            t,
            0,
			ct1.muLcoeff0_,
			ct1.muLcoeff1_,
			ct1.muLcoeff2_,
			ct1.muLcoeff3_,
			ct1.muLcoeff4_,
			ct1.muLcoeff5_,
			ct1.muVcoeff0_,
			ct1.muVcoeff1_,
			ct1.muVcoeff2_,
			ct1.muVcoeff3_,
			ct1.muVcoeff4_,
			ct1.muVcoeff5_,
			ct1.PrLcoeff0_,
			ct1.PrLcoeff1_,
			ct1.PrLcoeff2_,
			ct1.PrLcoeff3_,
			ct1.PrLcoeff4_,
			ct1.PrLcoeff5_,
			ct1.PrVcoeff0_,
			ct1.PrVcoeff1_,
			ct1.PrVcoeff2_,
			ct1.PrVcoeff3_,
			ct1.PrVcoeff4_,
			ct1.PrVcoeff5_,
			ct1.Trcoeff0_,
			ct1.Trcoeff1_,
			ct1.Trcoeff2_,
			ct1.Pcr_
        );
    }
    else
    {
        scalar molr1 = ct1.nMoles()/t.nMoles();
        scalar molr2 = ct2.nMoles()/t.nMoles();

        return superTransport<Thermo>
        (
            t,
            molr1*ct1.muLcoeff0_ + molr2*ct2.muLcoeff0_,
            molr1*ct1.muLcoeff1_ + molr2*ct2.muLcoeff1_,
            molr1*ct1.muLcoeff2_ + molr2*ct2.muLcoeff2_,
            molr1*ct1.muLcoeff3_ + molr2*ct2.muLcoeff3_,
            molr1*ct1.muLcoeff4_ + molr2*ct2.muLcoeff4_,
            molr1*ct1.muLcoeff5_ + molr2*ct2.muLcoeff5_,
            molr1*ct1.muVcoeff0_ + molr2*ct2.muVcoeff0_,
            molr1*ct1.muVcoeff1_ + molr2*ct2.muVcoeff1_,
            molr1*ct1.muVcoeff2_ + molr2*ct2.muVcoeff2_,
            molr1*ct1.muVcoeff3_ + molr2*ct2.muVcoeff3_,
            molr1*ct1.muVcoeff4_ + molr2*ct2.muVcoeff4_,
            molr1*ct1.muVcoeff5_ + molr2*ct2.muVcoeff5_,
            molr1*ct1.PrLcoeff0_ + molr2*ct2.PrLcoeff0_,
            molr1*ct1.PrLcoeff1_ + molr2*ct2.PrLcoeff1_,
            molr1*ct1.PrLcoeff2_ + molr2*ct2.PrLcoeff2_,
            molr1*ct1.PrLcoeff3_ + molr2*ct2.PrLcoeff3_,
            molr1*ct1.PrLcoeff4_ + molr2*ct2.PrLcoeff4_,
            molr1*ct1.PrLcoeff5_ + molr2*ct2.PrLcoeff5_,
            molr1*ct1.PrVcoeff0_ + molr2*ct2.PrVcoeff0_,
            molr1*ct1.PrVcoeff1_ + molr2*ct2.PrVcoeff1_,
            molr1*ct1.PrVcoeff2_ + molr2*ct2.PrVcoeff2_,
            molr1*ct1.PrVcoeff3_ + molr2*ct2.PrVcoeff3_,
            molr1*ct1.PrVcoeff4_ + molr2*ct2.PrVcoeff4_,
            molr1*ct1.PrVcoeff5_ + molr2*ct2.PrVcoeff5_,
            molr1*ct1.Trcoeff0_ + molr2*ct2.Trcoeff0_,
            molr1*ct1.Trcoeff1_ + molr2*ct2.Trcoeff1_,
            molr1*ct1.Trcoeff2_ + molr2*ct2.Trcoeff2_,
            molr1*ct1.Pcr_ + molr2*ct2.Pcr_
        );
    }
}


template<class Thermo>
inline Foam::superTransport<Thermo> Foam::operator-
(
    const superTransport<Thermo>& ct1,
    const superTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    if (mag(ct1.nMoles()) + mag(ct2.nMoles()) < SMALL)
    {
        return superTransport<Thermo>
        (
            t,
            0,
			ct1.muLcoeff0_,
			ct1.muLcoeff1_,
			ct1.muLcoeff2_,
			ct1.muLcoeff3_,
			ct1.muLcoeff4_,
			ct1.muLcoeff5_,
			ct1.muVcoeff0_,
			ct1.muVcoeff1_,
			ct1.muVcoeff2_,
			ct1.muVcoeff3_,
			ct1.muVcoeff4_,
			ct1.muVcoeff5_,
			ct1.PrLcoeff0_,
			ct1.PrLcoeff1_,
			ct1.PrLcoeff2_,
			ct1.PrLcoeff3_,
			ct1.PrLcoeff4_,
			ct1.PrLcoeff5_,
			ct1.PrVcoeff0_,
			ct1.PrVcoeff1_,
			ct1.PrVcoeff2_,
			ct1.PrVcoeff3_,
			ct1.PrVcoeff4_,
			ct1.PrVcoeff5_,
			ct1.Trcoeff0_,
			ct1.Trcoeff1_,
			ct1.Trcoeff2_,
			ct1.Pcr_
        );
    }
    else
    {
        scalar molr1 = ct1.nMoles()/t.nMoles();
        scalar molr2 = ct2.nMoles()/t.nMoles();

        return superTransport<Thermo>
        (
            t,
            molr1*ct1.muLcoeff0_ - molr2*ct2.muLcoeff0_,
            molr1*ct1.muLcoeff1_ - molr2*ct2.muLcoeff1_,
            molr1*ct1.muLcoeff2_ - molr2*ct2.muLcoeff2_,
            molr1*ct1.muLcoeff3_ - molr2*ct2.muLcoeff3_,
            molr1*ct1.muLcoeff4_ - molr2*ct2.muLcoeff4_,
            molr1*ct1.muLcoeff5_ - molr2*ct2.muLcoeff5_,
            molr1*ct1.muVcoeff0_ - molr2*ct2.muVcoeff0_,
            molr1*ct1.muVcoeff1_ - molr2*ct2.muVcoeff1_,
            molr1*ct1.muVcoeff2_ - molr2*ct2.muVcoeff2_,
            molr1*ct1.muVcoeff3_ - molr2*ct2.muVcoeff3_,
            molr1*ct1.muVcoeff4_ - molr2*ct2.muVcoeff4_,
            molr1*ct1.muVcoeff5_ - molr2*ct2.muVcoeff5_,
            molr1*ct1.PrLcoeff0_ - molr2*ct2.PrLcoeff0_,
            molr1*ct1.PrLcoeff1_ - molr2*ct2.PrLcoeff1_,
            molr1*ct1.PrLcoeff2_ - molr2*ct2.PrLcoeff2_,
            molr1*ct1.PrLcoeff3_ - molr2*ct2.PrLcoeff3_,
            molr1*ct1.PrLcoeff4_ - molr2*ct2.PrLcoeff4_,
            molr1*ct1.PrLcoeff5_ - molr2*ct2.PrLcoeff5_,
            molr1*ct1.PrVcoeff0_ - molr2*ct2.PrVcoeff0_,
            molr1*ct1.PrVcoeff1_ - molr2*ct2.PrVcoeff1_,
            molr1*ct1.PrVcoeff2_ - molr2*ct2.PrVcoeff2_,
            molr1*ct1.PrVcoeff3_ - molr2*ct2.PrVcoeff3_,
            molr1*ct1.PrVcoeff4_ - molr2*ct2.PrVcoeff4_,
            molr1*ct1.PrVcoeff5_ - molr2*ct2.PrVcoeff5_,
            molr1*ct1.Trcoeff0_ - molr2*ct2.Trcoeff0_,
            molr1*ct1.Trcoeff1_ - molr2*ct2.Trcoeff1_,
            molr1*ct1.Trcoeff2_ - molr2*ct2.Trcoeff2_,
            molr1*ct1.Pcr_ - molr2*ct2.Pcr_
        );
    }
}


template<class Thermo>
inline Foam::superTransport<Thermo> Foam::operator*
(
    const scalar s,
    const superTransport<Thermo>& ct
)
{
    return superTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
			ct.muLcoeff0_,
			ct.muLcoeff1_,
			ct.muLcoeff2_,
			ct.muLcoeff3_,
			ct.muLcoeff4_,
			ct.muLcoeff5_,
			ct.muVcoeff0_,
			ct.muVcoeff1_,
			ct.muVcoeff2_,
			ct.muVcoeff3_,
			ct.muVcoeff4_,
			ct.muVcoeff5_,
			ct.PrLcoeff0_,
			ct.PrLcoeff1_,
			ct.PrLcoeff2_,
			ct.PrLcoeff3_,
			ct.PrLcoeff4_,
			ct.PrLcoeff5_,
			ct.PrVcoeff0_,
			ct.PrVcoeff1_,
			ct.PrVcoeff2_,
			ct.PrVcoeff3_,
			ct.PrVcoeff4_,
			ct.PrVcoeff5_,
			ct.Trcoeff0_,
			ct.Trcoeff1_,
			ct.Trcoeff2_,
			ct.Pcr_
    );
}


template<class Thermo>
inline Foam::superTransport<Thermo> Foam::operator==
(
    const superTransport<Thermo>& ct1,
    const superTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
