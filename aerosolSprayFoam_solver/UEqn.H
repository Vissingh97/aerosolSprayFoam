// Momentum equation.  Dynamic-mesh-capable solver, no MRF.
// Mesh-motion flux effects are handled through phi after mesh.update().

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
 ==
    rho*g
  + parcels.SU(U)
  + fvOptions(rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
