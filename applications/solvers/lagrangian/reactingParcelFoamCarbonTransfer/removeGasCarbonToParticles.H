// Remove the mass represented by newly injected carbon parcels.
// The carrier-phase equations apply these sources during this time step.

carbonTransferMassRate.primitiveFieldRef() = scalar(0);
carbonTransferMomentumRate.primitiveFieldRef() = vector::zero;
carbonTransferEnergyRate.primitiveFieldRef() = scalar(0);

if (carbonToParticleCorrection)
{
    const label timeIdx = runTime.timeIndex();
    bool scheduled = false;
    label staggerEventIndex = -1;

    if (!carbonUsePhysicalTime)
    {
        scheduled =
            timeIdx >= carbonRemovalOffset
         && ((timeIdx - carbonRemovalOffset) % carbonRemovalInterval == 0);
    }
    else
    {
        staggerEventIndex = staggeredCarbonEventIndex
        (
            runTime.value(), runTime.deltaTValue(),
            carbonFirstTransferTime, carbonTransferInterval,
            carbonStaggerGroups
        );
        if (staggerEventIndex == -2)
        {
            FatalErrorInFunction
                << "Time step is larger than staggered sub-event spacing."
                << exit(FatalError);
        }
        scheduled = staggerEventIndex >= 0;
    }

    scalar localInjected = 0;
    forAll(carbonInjectedMass, celli)
    {
        localInjected += carbonInjectedMass[celli];
    }
    const scalar globalInjected = returnReduce(localInjected, sumOp<scalar>());

    if (!scheduled && globalInjected > ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Parcels injected at t=" << runTime.value()
            << " but carbon removal is not scheduled at this step."
            << exit(FatalError);
    }

    if (scheduled && globalInjected > ROOTVSMALL)
    {
        const volScalarField& Yc = Y[carbonSpecieIndex];
        const scalar invDt = 1.0/runTime.deltaTValue();
        scalar localRemoved = 0;
        scalar localUnmatched = 0;
        scalar localMaxRemovedFraction = 0;
        label localNCells = 0;

        forAll(carbonInjectedMass, celli)
        {
            const scalar injected = max(carbonInjectedMass[celli], scalar(0));
            if (injected <= ROOTVSMALL)
            {
                continue;
            }

            const scalar rho0 = rho[celli];
            const scalar cellV = mesh.V()[celli];
            const scalar yC0 = max(min(Yc[celli], scalar(1)), scalar(0));
            const scalar available = rho0*cellV*yC0;
            const scalar removed = min(injected, available);
            localUnmatched += injected - removed;

            if (removed <= ROOTVSMALL)
            {
                continue;
            }

            const scalar massFraction = removed/max(rho0*cellV, VSMALL);
            if (scalar(1) - massFraction <= minCarrierMassFraction)
            {
                FatalErrorInFunction
                    << "Insufficient carrier gas in cell " << celli
                    << exit(FatalError);
            }

            const scalar parcelEp =
                carbonInjectedInternalEnergy[celli]/(injected + VSMALL);
            const vector parcelUp =
                carbonInjectedMomentum[celli]/(injected + VSMALL);
            const scalar rate = removed/cellV*invDt;

            carbonTransferMassRate[celli] = rate;
            carbonTransferMomentumRate[celli] = rate*parcelUp;
            carbonTransferEnergyRate[celli] =
                rate*(parcelEp + 0.5*magSqr(parcelUp));

            localRemoved += removed;
            localMaxRemovedFraction = max(localMaxRemovedFraction, massFraction);
            ++localNCells;
        }

        const scalar globallyRemoved =
            returnReduce(localRemoved, sumOp<scalar>());
        const scalar globallyUnmatched =
            returnReduce(localUnmatched, sumOp<scalar>());

        if
        (
            globallyUnmatched >
            max(1e-9*globalInjected, scalar(1e-18))
        )
        {
            FatalErrorInFunction
                << "Carbon transfer unmatched mass: "
                << globallyUnmatched << " kg"
                << exit(FatalError);
        }

        Info<< "StaggerTransfer time=" << runTime.value()
            << " subEvent=" << staggerEventIndex
            << " group="
            << (staggerEventIndex >= 0
               ? staggerEventIndex % carbonStaggerGroups : 0)
            << " groups=" << carbonStaggerGroups
            << " cells=" << returnReduce(localNCells, sumOp<label>())
            << " injected=" << globalInjected
            << " removed=" << globallyRemoved
            << " unmatched=" << globallyUnmatched
            << " rate=" << globallyRemoved*invDt << " kg/s"
            << " maxCellRemovedFraction="
            << returnReduce(localMaxRemovedFraction, maxOp<scalar>())
            << nl;

    }
}

// Clear the parcel record after the source terms have been set.
carbonInjectedMass.primitiveFieldRef() = scalar(0);
carbonInjectedInternalEnergy.primitiveFieldRef() = scalar(0);
carbonInjectedMomentum.primitiveFieldRef() = vector::zero;
