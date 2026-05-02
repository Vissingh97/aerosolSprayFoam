// Compute aerosol removal rate coefficients from Lagrangian spray parcels.
//
// Paper equations implemented here:
//
// Inertial impaction:
//   St = rhoP*dP^2*|UG - UD|/(9*muG*dD)
//
//   eta_imp = 0                                               St <= 0.0833
//           = 8.57*(St/(St+0.5))^2*(St - 0.08336)             0.0833 < St < 0.2
//           = (St/(St+0.5))^2                                 St >= 0.2
//
// Interception:
//   eta_int = (1-alphaL)/(J + sigma*K) *
//             [ R/(1+R) + 0.5*(R/(1+R))^2*(3*sigma + 4) ]
//
//   J     = 1 - (6/5)*alphaL^(1/3) + (1/5)*alphaL^2
//   K     = 1 - (9/5)*alphaL^(1/3) + (1/5)*alphaL^2
//   sigma = muD/muG
//   R     = dP/dD
//
// Brownian diffusion:
//   eta_diff = (2*Pe*dD)^(-1/2)
//
//   Pe    = dD*|UG - UD|/Ddiff
//   Ddiff = kB*T*C/(3*pi*nuG*rhoG*dP)
//         = kB*T*C/(3*pi*muG*dP)
//   C     = 2.609*sqrt(2*l/dP), for 0.05 um < dP < 1.0 um
//
// Combined collection efficiency:
//   eta_total = 1 - (1 - eta_imp)*(1 - eta_int)*(1 - eta_diff)
//
// Paper Eq. (15) uses droplet velocity magnitude UD for swept volume:
//   dV/dt = pi/4*dD^2*|UD|
//
// The local first-order scavenging coefficient is:
//   K_i = sum_parcels(eta_total * pi/4*dD^2 * |UD| * nParticle / Vcell)
//
// The coefficient K_i [1/s] is later used in aerosolTransport.H as:
//   rhoDot_s,i = -rho*K_i*Yi
//
// Note:
//   The paper defines rhoDot_s,i as the aerosol sink term but does not
//   explicitly write rhoDot_s,i = -rho*K_i*Yi. That is the local first-order
//   implementation closure derived from the droplet capture rate.
// -------------------------------------------------------------------------

const scalar piVal = constant::mathematical::pi;

for (label bi = 0; bi < nBins; bi++)
{
    aerosolKeff[bi] = dimensionedScalar("zero", dimless/dimTime, Zero);
}

alphaL = dimensionedScalar("zero", dimless, Zero);

// Physical constants
const scalar rhoP      = 5890.0;        // ZrO2 density [kg/m3]
const scalar muD       = 1.0e-3;        // water dynamic viscosity [Pa s]
const scalar kB        = 1.380649e-23;  // Boltzmann constant [J/K]
const scalar lambdaMFP = 68.0e-9;       // air mean free path near room T, 1 atm [m]

const volScalarField muVol = thermo.mu();

const scalarField& muField  = muVol.primitiveField();
const scalarField& TField   = T.primitiveField();
const vectorField& UField   = U.primitiveField();
const scalarField& VField   = mesh.V();

// First pass: local liquid volume fraction from parcels.
forAllConstIter(basicSprayCloud, parcels, pIter)
{
    const auto& parcel = pIter();
    const label ci = parcel.cell();

    if (ci < 0) continue;

    const scalar dD = max(parcel.d(), SMALL);
    const scalar nP = max(parcel.nParticle(), 0.0);

    alphaL.primitiveFieldRef()[ci] +=
        nP*(piVal/6.0)*pow3(dD)/max(VField[ci], VSMALL);
}

alphaL.correctBoundaryConditions();

// Second pass: capture coefficient for each aerosol bin.
forAllConstIter(basicSprayCloud, parcels, pIter)
{
    const auto& parcel = pIter();
    const label ci = parcel.cell();

    if (ci < 0)
    {
        continue;
    }

    const scalar dD = max(parcel.d(), SMALL);
    const scalar nP = max(parcel.nParticle(), scalar(0));
    const scalar Vc = max(VField[ci], VSMALL);

    const scalar vRel   = max(mag(UField[ci] - parcel.U()), SMALL);
    // Strict paper Eq. (15) swept volume uses droplet speed magnitude.
    // For a more physical general form, vRel could be used here instead.
    const scalar vSweep = max(mag(parcel.U()), SMALL);  // paper Eq. (15)

    const scalar muG  = max(muField[ci], VSMALL);
    const scalar Tc   = max(TField[ci], SMALL);

    const scalar aL = min(max(alphaL.primitiveField()[ci], scalar(0)), scalar(0.999));
    const scalar aL13 = scalar(std::cbrt(static_cast<double>(aL)));

    const scalar J = 1.0 - (6.0/5.0)*aL13 + (1.0/5.0)*sqr(aL);
    const scalar K = 1.0 - (9.0/5.0)*aL13 + (1.0/5.0)*sqr(aL);
    const scalar sigma = muD/muG;

    for (label bi = 0; bi < nBins; bi++)
    {
        const scalar dP = aerosolDp[bi];

        // Eq. (9): modified Stokes number
        const scalar St = rhoP*sqr(dP)*vRel/(9.0*muG*dD);

        // Eq. (8): inertial impaction
        scalar etaImp = 0.0;
        if (St <= 0.0833)
        {
            etaImp = 0.0;
        }
        else if (St < 0.2)
        {
            etaImp = 8.57*sqr(St/(St + 0.5))*(St - 0.08336);
        }
        else
        {
            etaImp = sqr(St/(St + 0.5));
        }
        etaImp = min(max(etaImp, scalar(0)), scalar(1));

        // Eq. (10): interception
        const scalar R = dP/dD;
        const scalar RR = R/(1.0 + R);
        const scalar denom = max(J + sigma*K, VSMALL);

        scalar etaInt = ((1.0 - aL)/denom)
            *(RR + 0.5*sqr(RR)*(3.0*sigma + 4.0));
        etaInt = min(max(etaInt, scalar(0)), scalar(1));

        // Eq. (12): Cunningham correction as printed in the paper.
        const scalar cArg = max(2.0*lambdaMFP/dP, VSMALL);

        const scalar Cc = 2.609*scalar(std::sqrt(static_cast<double>(cArg)));

        // Eq. (11): Brownian diffusion
        const scalar Ddiff = kB*Tc*Cc/(3.0*piVal*muG*dP);
        const scalar Pe = dD*vRel/max(Ddiff, VSMALL);

        const scalar diffArg = max(2.0*Pe*dD, VSMALL);

        scalar etaDiff = 1.0/scalar(std::sqrt(static_cast<double>(diffArg)));

        etaDiff = min(max(etaDiff, scalar(0)), scalar(1));

        // Eq. (13): combined independent collection efficiency
        const scalar etaTotal = min
        (
            max(1.0 - (1.0 - etaImp)*(1.0 - etaInt)*(1.0 - etaDiff), scalar(0)),
            scalar(1)
        );

        // K_i [1/s] = eta_total * swept-area-rate * droplets-per-parcel / cell-volume
        aerosolKeff[bi].primitiveFieldRef()[ci] +=
            etaTotal*0.25*piVal*sqr(dD)*vSweep*nP/Vc;
    }
}

for (label bi = 0; bi < nBins; bi++)
{
    aerosolKeff[bi].correctBoundaryConditions();
}
// Diagnostic: domain-integrated aerosol mass removal rate
{
    scalar totalMassRemovalRate = 0.0;
    for (label bi = 0; bi < nBins; bi++)
    {
        totalMassRemovalRate += gSum
        (
            aerosolKeff[bi].primitiveField()
           *rho.primitiveField()
           *Y_aerosol[bi].primitiveField()
           *mesh.V()
        );
    }

    Info<< "Aerosol removal rate, all bins = "
        << totalMassRemovalRate*1e6 << " mg/s" << endl;
}

