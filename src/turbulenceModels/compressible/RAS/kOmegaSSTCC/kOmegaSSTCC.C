/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013  Dmitry Bogdanov
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "kOmegaSSTCC.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaSSTCC, 0);
addToRunTimeSelectionTable(RASModel, kOmegaSSTCC, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSSTCC::kOmegaSSTCC
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, rho, U, phi, thermophysicalModel, turbulenceModelName),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr1",
            coeffDict_,
            1.0
        )
    ),
    cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr2",
            coeffDict_,
            2.0
        )
    ),
    cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cr3",
            coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),

    y_(mesh_),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateMut("mut", mesh_)
    ),
    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateAlphat("alphat", mesh_)
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    mut_ =
    (
        a1_*rho_*k_
      / max
        (
            a1_*omega_,
            F2()*sqrt(2.0*magSqr(symm(fvc::grad(U_))))
        )
    );
    mut_.correctBoundaryConditions();

    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void kOmegaSSTCC::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ =
            a1_*rho_*k_
           /max(a1_*omega_, F2()*sqrt(2.0*magSqr(symm(fvc::grad(U_)))));
        mut_.correctBoundaryConditions();

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    volScalarField divU(fvc::div(phi_/fvc::interpolate(rho_)));

    if (mesh_.changing())
    {
        y_.correct();
    }

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    tmp<volTensorField> tSkew = skew(tgradU());
    tmp<volSymmTensorField> tSymm = symm(tgradU());
    volScalarField symInnerProduct = 2. * tSymm() && tSymm();
    volScalarField asymInnerProduct = max(2. * tSkew() && tSkew(), dimensionedScalar("1e-16", dimensionSet(0, 0, -2, 0, 0), 1e-10) );
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField GbyMu((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField rStar = sqrt(symInnerProduct/asymInnerProduct);
    volScalarField G("RASModel::G", mut_*GbyMu);
    tgradU.clear();
    volScalarField D = sqrt(max(asymInnerProduct, 0.09*omega_*omega_));
    omega_.boundaryField().updateCoeffs();
    tmp<volSymmTensorField> divS = fvc::ddt(tSymm()) + fvc::div(surfaceScalarField("phiU",phi_/fvc::interpolate(rho_)), tSymm());
    volScalarField rT = tSkew().component(0)*tSymm().component(0)*divS().component(0) +
                            tSkew().component(0)*tSymm().component(1)*divS().component(1) +
                            tSkew().component(0)*tSymm().component(2)*divS().component(2) +
                            tSkew().component(3)*tSymm().component(0)*divS().component(3) +
                            tSkew().component(3)*tSymm().component(1)*divS().component(4) +
                            tSkew().component(3)*tSymm().component(2)*divS().component(5) +
                            tSkew().component(6)*tSymm().component(0)*divS().component(6) +
                            tSkew().component(6)*tSymm().component(1)*divS().component(7) +
                            tSkew().component(6)*tSymm().component(2)*divS().component(8) +
                            tSkew().component(1)*tSymm().component(3)*divS().component(0) +
                            tSkew().component(1)*tSymm().component(4)*divS().component(1) +
                            tSkew().component(1)*tSymm().component(5)*divS().component(2) +
                            tSkew().component(4)*tSymm().component(3)*divS().component(3) +
                            tSkew().component(4)*tSymm().component(4)*divS().component(4) +
                            tSkew().component(4)*tSymm().component(5)*divS().component(5) +
                            tSkew().component(7)*tSymm().component(3)*divS().component(6) +
                            tSkew().component(7)*tSymm().component(4)*divS().component(7) +
                            tSkew().component(7)*tSymm().component(5)*divS().component(8) +
                            tSkew().component(2)*tSymm().component(6)*divS().component(0) +
                            tSkew().component(2)*tSymm().component(7)*divS().component(1) +
                            tSkew().component(2)*tSymm().component(8)*divS().component(2) +
                            tSkew().component(7)*tSymm().component(6)*divS().component(3) +
                            tSkew().component(7)*tSymm().component(7)*divS().component(4) +
                            tSkew().component(7)*tSymm().component(8)*divS().component(5) +
                            tSkew().component(8)*tSymm().component(6)*divS().component(6) +
                            tSkew().component(8)*tSymm().component(7)*divS().component(7) +
                            tSkew().component(8)*tSymm().component(8)*divS().component(8);
    divS.clear();
    tSkew.clear();
    tSymm.clear();
    volScalarField rTilda = 2. * rT/sqrt(asymInnerProduct)/D/D/D;
    volScalarField Frot(
        IOobject(
            "Frot",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        max(min((scalar(1) + cr1_)*2*rStar/(scalar(1)+rStar)*(scalar(1)-cr3_*atan(cr2_*rTilda))-cr1_, scalar(1.25)), scalar(0))
    );
    if(runTime_.outputTime())
    {
        Frot.write();
    }
    rStar.clear();
    rTilda.clear();
    rT.clear();
    D.clear();
    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField rhoGammaF1(rho_*gamma(F1));

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(rho_, omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        rhoGammaF1*GbyMu*Frot
      - fvm::SuSp((2.0/3.0)*rhoGammaF1*divU, omega_)
      - fvm::Sp(rho_*beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            rho_*(F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G*Frot, (c1_*betaStar_)*rho_*k_*omega_)
      - fvm::SuSp(2.0/3.0*rho_*divU, k_)
      - fvm::Sp(rho_*betaStar_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    mut_ = a1_*rho_*k_/max(a1_*omega_, F2()*sqrt(S2));
    mut_.correctBoundaryConditions();

    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
