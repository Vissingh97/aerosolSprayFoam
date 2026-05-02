// -------------------------------------------------------------------------
// Eulerian aerosol-bin passive scalar transport, paper-consistent form.
//
// Liang et al. solve aerosol particle bins as gas species/passive scalars:
//
//   d(rho*Yi)/dt + div(rho*U*Yi) - div(muEff grad Yi) = rhoDot_s,i
//
// Here the droplet scavenging sink is written as a first-order loss:
//
//   rhoDot_s,i = -rho*K_i*Yi
//
// where:
//   Yi          aerosol mass fraction for bin i [-]
//   rho*Yi      aerosol mass concentration [kg/m3]
//   K_i         droplet scavenging coefficient [1/s]
//
// This replaces the previous number-concentration formulation and is the
// strict paper-consistent equation for benchmark reproduction.
// -------------------------------------------------------------------------

const volScalarField aerosolMuEff
(
    IOobject
    (
        "aerosolMuEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*turbulence->nuEff()
);

for (label bi = 0; bi < nBins; bi++)
{
    fvScalarMatrix YAEqn
    (
        fvm::ddt(rho, Y_aerosol[bi])
      + fvm::div(phi, Y_aerosol[bi])
      - fvm::laplacian(aerosolMuEff, Y_aerosol[bi])
     ==
      - fvm::Sp(rho*aerosolKeff[bi], Y_aerosol[bi])
      + fvOptions(rho, Y_aerosol[bi])
    );

    YAEqn.relax();
    fvOptions.constrain(YAEqn);
    YAEqn.solve();
    fvOptions.correct(Y_aerosol[bi]);

    Y_aerosol[bi].max(0.0);
}

if (runTime.writeTime())
{
    scalar Mtotal = 0.0;

    Info<< "Aerosol bin mass from rho*Yi" << nl;

    for (label bi = 0; bi < nBins; bi++)
    {
        const scalar Mbin =
            gSum
            (
                rho.primitiveField()
               *Y_aerosol[bi].primitiveField()
               *mesh.V()
            );

        Mtotal += Mbin;

        Info<< "  aerosolBin" << bi
            << " dp=" << aerosolDp[bi]*1e6 << " um: "
            << Mbin*1e6 << " mg" << nl;
    }

    Info<< "  Total aerosol mass = " << Mtotal*1e6 << " mg" << nl << endl;
}
