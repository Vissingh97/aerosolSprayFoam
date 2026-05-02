/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    aerosolSprayFoam

Description
    Minimal Euler-Lagrange aerosol washout solver for reproducing Liang et al.
    (2020) aerosol removal by water spray droplets.

    Paper-consistent benchmark simplifications:
      - no energy equation
      - no combustion
      - no radiation
      - no thermophoresis
      - no diffusiophoresis
      - dynamic mesh retained for future moving/deforming-mesh cases
      - no MRF

    Gas is solved with compressible PIMPLE equations. Spray droplets are tracked
    with basicSprayCloud. Aerosol bins are Eulerian passive scalars with a
    first-order capture sink from the Lagrangian droplets.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "basicSprayCloud.H"
#include "psiReactionThermo.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Aerosol washout by Lagrangian spray droplets: no energy, "
        "no combustion, no radiation."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readDyMControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Optional dynamic-mesh update retained.
        // For the Liang et al. benchmark use staticFvMesh in dynamicMeshDict.
        mesh.update();

        if (mesh.changing())
        {
         Info<< "Mesh changed. Flux correction is currently disabled in this solver."
            << nl
            << "For the Liang benchmark, use staticFvMesh so this branch is not active."
             << endl;
        }

        parcels.evolve();

        // Compute Eulerian aerosol sink coefficients from current droplets.
        #include "computeAerosolSink.H"

        if (pimple.solveFlow())
        {
            #include "rhoEqn.H"

            while (pimple.loop())
            {
                #include "UEqn.H"

                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();

            // Aerosol passive-scalar transport with droplet-capture sink.
            #include "aerosolTransport.H"
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
