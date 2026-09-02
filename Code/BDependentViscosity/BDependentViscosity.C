/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd
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

#include "BDependentViscosity.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(BDependentViscosity, 0);
    addToRunTimeSelectionTable(viscosityModel, BDependentViscosity, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::BDependentViscosity::calcNu() const
{
    
    const volScalarField& UM= U_.mesh().lookupObject<volScalarField>("UM");
    
    dimensionedScalar gamma = 3.0*(muM2_/muM0_-1.0)/Ms_; //m/A
    dimensionedScalar SMALLH = 1e-20*Ms_;
    
    volVectorField B = (((Ms_/(SMALLH+mag(fvc::grad(UM)))*(1.0/tanh(SMALL+gamma*mag(fvc::grad(UM)))-1.0/(SMALL+gamma*mag(fvc::grad(UM)))))+1.0)*muM0_)*fvc::grad(UM);
    
    return a_*mag(B)/(mag(B)+b_)+c_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::BDependentViscosity::BDependentViscosity
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    polynomialCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    
    a_(polynomialCoeffs_.lookup("a")),
    b_(polynomialCoeffs_.lookup("b")),
    c_(polynomialCoeffs_.lookup("c")),
    muM2_(polynomialCoeffs_.lookup("muM2")),
    muM0_(polynomialCoeffs_.lookup("muM0")),
    Ms_(polynomialCoeffs_.lookup("Ms")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::BDependentViscosity::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    polynomialCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    polynomialCoeffs_.lookup("a") >> a_;
    polynomialCoeffs_.lookup("b") >> b_;
    polynomialCoeffs_.lookup("c") >> c_;
    polynomialCoeffs_.lookup("muM2") >> muM2_;
    polynomialCoeffs_.lookup("muM0") >> muM0_;
    polynomialCoeffs_.lookup("Ms") >> Ms_;

    return true;
}


// ************************************************************************* //
