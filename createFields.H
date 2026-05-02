#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();

// The benchmark paper does not solve the gas energy equation.
// Thermo is still created because rho, psi, mu, T and p are needed by the
// compressible gas solver, turbulence model, Brownian diffusion coefficient,
// and Lagrangian spray cloud.
//
// This validation line may be required by psiReactionThermo even though EEqn.H
// is not solved.
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

volScalarField& p = thermo.p();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

mesh.setFluxRequired(p.name());

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Keep only if K is used by your pEqn.H or other solver files.
#include "createK.H"

#include "createClouds.H"
#include "createFvOptions.H"

// ---- Eulerian aerosol fields ---------------------------------------------
// Nine size bins AP2-AP10 from Liang et al. benchmark.
// These fields are aerosol mass fractions Yi [-], not number concentrations.
// The transported quantity is rho*Yi [kg/m3].
//
// Mapping:
//   aerosolBin0 = AP2  = 0.198 um
//   aerosolBin1 = AP3  = 0.305 um
//   aerosolBin2 = AP4  = 0.407 um
//   aerosolBin3 = AP5  = 0.505 um
//   aerosolBin4 = AP6  = 0.583 um
//   aerosolBin5 = AP7  = 0.724 um
//   aerosolBin6 = AP8  = 0.778 um
//   aerosolBin7 = AP9  = 0.899 um
//   aerosolBin8 = AP10 = 1.038 um

const label nBins = 9;

const scalarList aerosolDp
{
    0.198e-6,
    0.305e-6,
    0.407e-6,
    0.505e-6,
    0.583e-6,
    0.724e-6,
    0.778e-6,
    0.899e-6,
    1.038e-6
};

PtrList<volScalarField> Y_aerosol(nBins);

for (label bi = 0; bi < nBins; bi++)
{
    Y_aerosol.set
    (
        bi,
        new volScalarField
        (
            IOobject
            (
                "aerosolBin" + Foam::name(bi),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

// First-order removal coefficient K_i [1/s]. Rebuilt every time step from
// parcel positions, diameters and velocities in computeAerosolSink.H.
PtrList<volScalarField> aerosolKeff(nBins);

for (label bi = 0; bi < nBins; bi++)
{
    aerosolKeff.set
    (
        bi,
        new volScalarField
        (
            IOobject
            (
                "aerosolKeff" + Foam::name(bi),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless/dimTime, Zero)
        )
    );
}

// Local liquid volume fraction reconstructed from Lagrangian droplets.
// Used only in the Park et al. interception efficiency correlation.
volScalarField alphaL
(
    IOobject
    (
        "alphaL",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, Zero)
);
