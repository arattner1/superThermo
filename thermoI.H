/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "thermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type>::thermo
(
    const Thermo& sp
)
:
    Thermo(sp)
{}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::T
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (thermo<Thermo, Type>::*F)(const scalar, const scalar) const,
    scalar (thermo<Thermo, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (thermo<Thermo, Type>::*limit)(const scalar) const
) const
{
//Info<< endl <<"f (h): " << f << endl;
//Info<<"T0: " << T0 << endl;

    // Brent method, based on implementation here:
    // https://gist.github.com/thomasantony/6471125
    scalar a = 1;    //Lower bound for temperature
    scalar b = 1000; //Upper bound for temperature
    const scalar tol = T0*tol_;

    const scalar f_a = (this->*F)(p, a) - f;
    const scalar f_b = (this->*F)(p, b) - f;

    scalar s, c, f_s, f_c;

    //Counter added by MN
	int	count = 0;

	do{

		c = (a+b)/2.0;
		f_c = (this->*F)(p, c) - f;
//Info<< "c: " << c << endl;
//Info<< "f_a, f_c, f_b: " << f_a << ", " << f_c << ", " << f_b << endl;

		if(mag(f_a-f_c) > tol && mag(f_b-f_c) > tol) 	// f(a)!=f(c) and f(b)!=f(c)
		{
			// Inverse quadratic interpolation
			s = a*f_b*f_c/((f_a-f_b)*(f_a-f_c)) + b*f_a*f_c/((f_b-f_a)*(f_b-f_c)) + c*f_a*f_b/((f_c-f_a)*(f_c-f_b));
//Info<< "Quadratic, s: " << s << endl;
		}else{
			// Secant rule
			s = b - f_b*(b-a)/(f_b-f_a);
//Info<< "Secant, s: " << s << endl;
		}
		f_s = (this->*F)(p, s) - f;
//Info<< "f_s: " << f_s << endl;
		if(c>s)
		{
            scalar tmp = s;
            s = c;
            c = tmp;
            //Also need to flip f_c and f_s here?
            tmp = f_s;
            f_s = f_c;
            f_c = tmp;
		}
		if(f_c*f_s < 0)
		{
			a = c;
			b = s;
		}else if(f_s*f_b < 0){
			a = s;
		}else{
			b = c;
		}


		count = count + 1;
		if(count > 500)
		{ 
			Info<< "count: " << count << endl;
			Info<< "Iterations on temperature not successful" << endl;
		}

//Info<< "a, b: " << a << ", " << b << endl;
	}while(mag(f_s) > tol && mag(b-a)>tol);	// Convergence conditions
//Info<< "T: " << s << endl;
	return s;

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type>::thermo
(
    const word& name,
    const thermo& st
)
:
    Thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::word
Foam::species::thermo<Thermo, Type>::heName()
{
    return Type<thermo<Thermo, Type>>::name();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::he(const scalar p, const scalar T) const
{
    return Type<thermo<Thermo, Type>>::he(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::cv(const scalar p, const scalar T) const
{
    return this->cp(p, T) - this->cpMcv(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::cpv(const scalar p, const scalar T) const
{
    return Type<thermo<Thermo, Type>>::cpv(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::gamma(const scalar p, const scalar T) const
{
    scalar cp = this->cp(p, T);
    return cp/(cp - this->cpMcv(p, T));
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::cpBycpv
(
    const scalar p,
    const scalar T
) const
{
    return Type<thermo<Thermo, Type>>::cpBycpv(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::es(const scalar p, const scalar T) const
{
    return this->hs(p, T) - p*this->W()/this->rho(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::ea(const scalar p, const scalar T) const
{
    return this->ha(p, T) - p*this->W()/this->rho(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::g(const scalar p, const scalar T) const
{
    return this->ha(p, T) - T*this->s(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::a(const scalar p, const scalar T) const
{
    return this->ea(p, T) - T*this->s(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Cpv(const scalar p, const scalar T) const
{
    return this->cpv(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Cp(const scalar p, const scalar T) const
{
    return this->cp(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Cv(const scalar p, const scalar T) const
{
    return this->cv(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::HE(const scalar p, const scalar T) const
{
    return Type<thermo<Thermo, Type>>::HE(*this, p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Hs(const scalar p, const scalar T) const
{
    return this->hs(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Hc() const
{
    return this->hc()/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Ha(const scalar p, const scalar T) const
{
    return this->ha(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::S(const scalar p, const scalar T) const
{
    return this->s(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::E(const scalar p, const scalar T) const
{
    return this->e(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Es(const scalar p, const scalar T) const
{
    return this->es(p, T)/this->W();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Ea(const scalar p, const scalar T) const
{
    return this->ea(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::G(const scalar p, const scalar T) const
{
    return this->g(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::A(const scalar p, const scalar T) const
{
    return this->a(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -this->nMoles()*this->g(Pstd, T)/(RR*T);

    if (arg < 600.0)
    {
        return exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo<Thermo, Type>::Kc(const scalar p, const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/(RR*T), this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::Kx
(
    const scalar p,
    const scalar T
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/p, this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::Kn
(
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(n*Pstd/p, this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::THE
(
    const scalar he,
    const scalar p,
    const scalar T0
) const
{
    return Type<thermo<Thermo, Type>>::THE(*this, he, p, T0);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::THs
(
    const scalar hs,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hs,
        p,
        T0,
        &thermo<Thermo, Type>::Hs,
        &thermo<Thermo, Type>::Cp,
        &thermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::THa
(
    const scalar ha,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ha,
        p,
        T0,
        &thermo<Thermo, Type>::Ha,
        &thermo<Thermo, Type>::Cp,
        &thermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::TEs
(
    const scalar es,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        es,
        p,
        T0,
        &thermo<Thermo, Type>::Es,
        &thermo<Thermo, Type>::Cv,
        &thermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo<Thermo, Type>::TEa
(
    const scalar ea,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ea,
        p,
        T0,
        &thermo<Thermo, Type>::Ea,
        &thermo<Thermo, Type>::Cv,
        &thermo<Thermo, Type>::limit
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline void Foam::species::thermo<Thermo, Type>::operator+=
(
    const thermo<Thermo, Type>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::species::thermo<Thermo, Type>::operator-=
(
    const thermo<Thermo, Type>& st
)
{
    Thermo::operator-=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::species::thermo<Thermo, Type>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type> Foam::species::operator+
(
    const thermo<Thermo, Type>& st1,
    const thermo<Thermo, Type>& st2
)
{
    return thermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type> Foam::species::operator-
(
    const thermo<Thermo, Type>& st1,
    const thermo<Thermo, Type>& st2
)
{
    return thermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type> Foam::species::operator*
(
    const scalar s,
    const thermo<Thermo, Type>& st
)
{
    return thermo<Thermo, Type>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo<Thermo, Type> Foam::species::operator==
(
    const thermo<Thermo, Type>& st1,
    const thermo<Thermo, Type>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
