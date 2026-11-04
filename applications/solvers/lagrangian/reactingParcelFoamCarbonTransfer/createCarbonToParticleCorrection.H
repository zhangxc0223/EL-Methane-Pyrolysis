/*---------------------------------------------------------------------------*\
  Set up carbon transfer between the gas and the particle cloud.
\*---------------------------------------------------------------------------*/

IOdictionary carbonToParticleProperties
(
    IOobject
    (
        "carbonToParticleProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const Switch carbonToParticleCorrection
(
    carbonToParticleProperties.lookupOrDefault<Switch>("active", false)
);

const word carbonSpecieName
(
    carbonToParticleProperties.lookupOrDefault<word>("carbonSpecie", "C")
);

const label carbonRemovalInterval
(
    carbonToParticleProperties.lookupOrDefault<label>("removeEveryNTimeSteps", 0)
);

const label carbonRemovalOffset
(
    carbonToParticleProperties.lookupOrDefault<label>("triggerOffset", 0)
);

// Use physical time to keep transfer events aligned after a restart.
const Switch carbonUsePhysicalTime
(
    carbonToParticleProperties.lookupOrDefault<Switch>
    (
        "usePhysicalTransferSchedule", false
    )
);
const scalar carbonFirstTransferTime
(
    carbonToParticleProperties.lookupOrDefault<scalar>
    (
        "firstTransferTime", 0.1
    )
);
const scalar carbonTransferInterval
(
    carbonToParticleProperties.lookupOrDefault<scalar>
    (
        "transferInterval", 0.1
    )
);

const scalar minCarrierMassFraction
(
    carbonToParticleProperties.lookupOrDefault<scalar>("minCarrierMassFraction", VSMALL)
);

// Use this reference temperature for the particle internal energy.
const scalar carbonParticleCp
(
    carbonToParticleProperties.lookupOrDefault<scalar>("particleCp", 700.0)
);
const scalar carbonParticleEnergyReferenceT
(
    carbonToParticleProperties.lookupOrDefault<scalar>
    (
        "particleEnergyReferenceTemperature", 298.15
    )
);
// Carbon mass added to each cell as parcels during this time step (kg).
volScalarField carbonInjectedMass
(
    IOobject
    (
        "carbonInjectedMass",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroCarbonInjectedMass", dimMass, 0.0)
);

// Mass-weighted particle energy and momentum from the injection model.
volScalarField carbonInjectedInternalEnergy
(
    IOobject
    (
        "carbonInjectedInternalEnergy", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroCarbonParticleE", dimEnergy, 0)
);

volVectorField carbonInjectedMomentum
(
    IOobject
    (
        "carbonInjectedMomentum", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zeroCarbonParticleMomentum", dimMass*dimVelocity, Zero)
);

// Carbon transfer sources for the gas equations.
volScalarField carbonTransferMassRate
(
    IOobject
    (
        "carbonTransferMassRate", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroCarbonMassRate", dimDensity/dimTime, 0)
);

volVectorField carbonTransferMomentumRate
(
    IOobject
    (
        "carbonTransferMomentumRate", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "zeroCarbonMomentumRate", dimDensity*dimVelocity/dimTime, Zero
    )
);

volScalarField carbonTransferEnergyRate
(
    IOobject
    (
        "carbonTransferEnergyRate", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zeroCarbonEnergyRate", dimEnergy/dimVolume/dimTime, 0
    )
);

label carbonSpecieIndex = -1;

if (carbonToParticleCorrection)
{
    if (!composition.species().found(carbonSpecieName))
    {
        FatalIOErrorIn(args.executable().c_str(), carbonToParticleProperties)
            << "carbonSpecie " << carbonSpecieName
            << " not found in available species "
            << composition.species()
            << exit(FatalIOError);
    }

    carbonSpecieIndex = composition.species()[carbonSpecieName];

    if (carbonSpecieIndex == inertIndex)
    {
        FatalIOErrorIn(args.executable().c_str(), carbonToParticleProperties)
            << "carbonSpecie " << carbonSpecieName
            << " is the inertSpecie. The carbon species must be solved as an "
            << "active species. Select the inertSpecie according to the gas mixture."
            << exit(FatalIOError);
    }

    if (!composition.active(carbonSpecieIndex))
    {
        FatalIOErrorIn(args.executable().c_str(), carbonToParticleProperties)
            << "carbonSpecie " << carbonSpecieName
            << " exists, but is not active in the reacting mixture."
            << exit(FatalIOError);
    }

    if (carbonParticleCp <= 0)
    {
        FatalErrorInFunction
            << "particleCp must be positive."
            << exit(FatalError);
    }

    if (carbonUsePhysicalTime && carbonTransferInterval <= 0)
    {
        FatalErrorInFunction
            << "transferInterval must be positive."
            << exit(FatalError);
    }

    if (!carbonUsePhysicalTime && carbonRemovalInterval <= 0)
    {
        FatalIOErrorIn(args.executable().c_str(), carbonToParticleProperties)
            << "removeEveryNTimeSteps must be a positive integer when "
            << "carbonToParticleCorrection is active."
            << exit(FatalIOError);
    }

    Info<< "Gas-carbon-to-particle correction is active:" << nl
        << "    carbonSpecie          = " << carbonSpecieName << nl
        << "    removeEveryNTimeSteps = " << carbonRemovalInterval << nl
        << "    triggerOffset         = " << carbonRemovalOffset << nl
        << "    particleCp            = " << carbonParticleCp << " J/(kg K)" << nl
        << "    solid energy Tref     = " << carbonParticleEnergyReferenceT << " K" << nl
        << "    physical-time events  = " << carbonUsePhysicalTime << nl
        << endl;
}
