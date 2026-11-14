/*---------------------------------------------------------------------------*\
  Optional gas mass accounting. This does not change the flow solution.

  The reported residual includes the custom gas-to-parcel carbon sink.
  Other mass sources (native cloud, surface film, fvOptions) are not counted.
  It is not a complete multiphase conservation test when these are active.
\*---------------------------------------------------------------------------*/

if
(
    carbonMassAuditEnabled
 && (runTime.timeIndex() % carbonMassAuditEveryN == 0 || runTime.outputTime())
)
{
    const scalar Mnew = fvc::domainIntegrate(rho).value();
    const scalar Mold = fvc::domainIntegrate(rho.oldTime()).value();
    const scalar dt = runTime.deltaTValue();
    scalar boundaryPhi = 0;
    scalar inletPhi = 0;
    scalar outletPhi = 0;
    forAll(mesh.boundary(), patchi)
    {
        // Global coupled processor boundaries are not external surfaces.
        if (!mesh.boundary()[patchi].coupled())
        {
            const scalar flux = sum(phi.boundaryField()[patchi]);
            boundaryPhi += flux;
            if (mesh.boundary()[patchi].name() == carbonMassAuditInlet)
            {
                inletPhi += flux;
            }
            if (mesh.boundary()[patchi].name() == carbonMassAuditOutlet)
            {
                outletPhi += flux;
            }
        }
    }
    reduce(boundaryPhi, sumOp<scalar>());
    reduce(inletPhi, sumOp<scalar>());
    reduce(outletPhi, sumOp<scalar>());

    const scalar carbonSinkRate =
        gSum(carbonTransferMassRate.primitiveField()*mesh.V());
    const scalar dMdt = (Mnew - Mold)/dt;
    const scalar residual = dMdt + boundaryPhi + carbonSinkRate;
    const label eIdx = carbonUsePhysicalTime ? staggeredCarbonEventIndex
    (
        runTime.value(), dt, carbonFirstTransferTime,
        carbonTransferInterval, carbonStaggerGroups
    ) : -1;

    Info<< "MassAudit time=" << runTime.value()
        << " Mnew=" << Mnew
        << " Mold=" << Mold
        << " dMdt=" << dMdt
        << " boundaryPhi=" << boundaryPhi
        << " inletPhi=" << inletPhi
        << " outletPhi=" << outletPhi
        << " carbonSinkRate=" << carbonSinkRate
        << " eventGroup=" << (eIdx >= 0 ? eIdx % carbonStaggerGroups : -1)
        << " residual=" << residual
        << " kg/s" << endl;
}
