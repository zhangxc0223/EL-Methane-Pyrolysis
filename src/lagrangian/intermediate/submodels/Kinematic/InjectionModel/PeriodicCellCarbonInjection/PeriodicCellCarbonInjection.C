/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This custom model is distributed under the GNU General Public License v3
    or later, consistently with OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "PeriodicCellCarbonInjection.H"
#include "StaggeredCarbonSchedule.H"
#include "volFields.H"
#include "polyMeshTetDecomposition.H"
#include "globalIndex.H"
#include "mathematicalConstants.H"
#include "PstreamReduceOps.H"
#include "IOdictionary.H"
#include <cmath>

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
inline Foam::label Foam::PeriodicCellCarbonInjection<CloudType>::localActiveIndex
(
    const label parcelI
) const
{
    const label globalCellI = parcelI/parcelsPerCell_;
    const label localI = globalCellI - globalCellStart_;
    if (localI < 0 || localI >= activeCells_.size())
    {
        return -1;
    }
    return localI;
}


template<class CloudType>
inline Foam::label Foam::PeriodicCellCarbonInjection<CloudType>::cellGroup
(
    const label celli
) const
{
    return staggeredCarbonCellGroup
    (
        this->owner().mesh().C()[celli],
        staggerGroups_,
        staggerHashResolution_
    );
}


template<class CloudType>
void Foam::PeriodicCellCarbonInjection<CloudType>::updateInjectionData()
{
    const fvMesh& mesh = this->owner().mesh();

    if (useLocalU_ && UFieldPtr_)
    {
        uInterp_ = interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(), *UFieldPtr_
        );
    }
    if (useLocalT_ && TFieldPtr_)
    {
        tInterp_ = interpolation<scalar>::New
        (
            this->owner().solution().interpolationSchemes(), *TFieldPtr_
        );
    }

    // Work with cells on this processor.
    activeCells_.setSize(mesh.nCells());
    activeMassPerParcel_.setSize(mesh.nCells());
    activeDiameters_.setSize(mesh.nCells());

    label nActive = 0;
    scalar localCarbonVolume = 0;
    scalar localFilteredCarbonMass = 0;
    label nFiltered = 0;
    label localGroupCells = 0;

    forAll(referenceField_, celli)
    {
        // Only this group's cells may transfer gas carbon in this step.
        if (cellGroup(celli) != currentStaggerGroup_)
        {
            continue;
        }
        ++localGroupCells;

        const scalar Yc = max(referenceField_[celli], scalar(0));
        const scalar carbonMass =
            factor_*Yc*max(rhoField_[celli], scalar(0))*mesh.V()[celli];

        if
        (
            Yc <= minCarbonMassFraction_
         || carbonMass <= minCarbonMass_
        )
        {
            localFilteredCarbonMass += carbonMass;
            ++nFiltered;
            continue;
        }

        const scalar d = sizeDistribution_->sample();
        const scalar Vp = pi/6.0*pow3(d);
        if
        (
            d <= ROOTVSMALL
         || Vp <= ROOTVSMALL
         || carbonMass/(scalar(parcelsPerCell_)*rhoParticle_*Vp)
            < this->minParticlesPerParcel_
        )
        {
            localFilteredCarbonMass += carbonMass;
            ++nFiltered;
            continue;
        }

        activeCells_[nActive] = celli;
        activeMassPerParcel_[nActive] = carbonMass/scalar(parcelsPerCell_);
        activeDiameters_[nActive] = d;
        localCarbonVolume += carbonMass/rhoParticle_;
        ++nActive;
    }

    activeCells_.setSize(nActive);
    activeMassPerParcel_.setSize(nActive);
    activeDiameters_.setSize(nActive);

    // Number active cells across all processors.
    const globalIndex globalCells(activeCells_.size());
    globalCellStart_ = globalCells.localStart();
    globalActiveCells_ = globalCells.totalSize();
    globalGroupCells_ = returnReduce(localGroupCells, sumOp<label>());

    volumeForEvent_ = returnReduce(localCarbonVolume, sumOp<scalar>());
    this->volumeTotal_ = max(volumeForEvent_, ROOTVSMALL);

    Info<< "StaggerInject group=" << currentStaggerGroup_
        << " groups=" << staggerGroups_
        << " groupCells=" << globalGroupCells_
        << " activeCells=" << globalActiveCells_
        << " parcelsExpected=" << globalActiveCells_*parcelsPerCell_
        << " mass=" << volumeForEvent_*rhoParticle_
        << " kg" << nl
        << "    periodicCellCarbonInjection: active cells = "
        << globalActiveCells_ << ", skipped cells = "
        << returnReduce(nFiltered, sumOp<label>())
        << ", skipped C mass = "
        << returnReduce(localFilteredCarbonMass, sumOp<scalar>())
        << " kg, event mass = " << volumeForEvent_*rhoParticle_
        << " kg" << nl;
}


template<class CloudType>
template<class ParcelType>
auto Foam::PeriodicCellCarbonInjection<CloudType>::setParcelTemperature
(
    ParcelType& parcel,
    const scalar T,
    const int
) -> decltype(parcel.T() = T, void())
{
    parcel.T() = T;
}


template<class CloudType>
template<class ParcelType>
void Foam::PeriodicCellCarbonInjection<CloudType>::setParcelTemperature
(
    ParcelType&,
    const scalar,
    const long
)
{}


template<class CloudType>
template<class ParcelType>
auto Foam::PeriodicCellCarbonInjection<CloudType>::parcelTemperature
(
    const ParcelType& parcel, const int
) -> decltype(parcel.T(), scalar())
{
    return parcel.T();
}


template<class CloudType>
template<class ParcelType>
Foam::scalar Foam::PeriodicCellCarbonInjection<CloudType>::parcelTemperature
(
    const ParcelType&, const long
)
{
    return scalar(0);
}


// * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PeriodicCellCarbonInjection<CloudType>::PeriodicCellCarbonInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    referenceField_
    (
        owner.db().objectRegistry::template lookupObject<volScalarField>
        (
            this->coeffDict().getWord("referenceField")
        )
    ),
    rhoField_
    (
        owner.db().objectRegistry::template lookupObject<volScalarField>
        (
            this->coeffDict().template lookupOrDefault<word>
            (
                "rhoField", "rho"
            )
        )
    ),
    factor_
    (
        this->coeffDict().template lookupOrDefault<scalar>("factor", 1)
    ),
    rhoParticle_
    (
        this->coeffDict().template get<scalar>("rhoParticle")
    ),
    minCarbonMassFraction_
    (
        this->coeffDict().template lookupOrDefault<scalar>
        (
            "minCarbonMassFraction", 0
        )
    ),
    minCarbonMass_
    (
        this->coeffDict().template lookupOrDefault<scalar>
        (
            "minCarbonMass", 0
        )
    ),
    parcelsPerCell_
    (
        this->coeffDict().template lookupOrDefault<label>
        (
            "parcelsPerCell", 1
        )
    ),
    firstInjectionTime_
    (
        this->coeffDict().template get<scalar>("firstInjectionTime")
    ),
    injectionInterval_
    (
        this->coeffDict().template get<scalar>("injectionInterval")
    ),
    staggerGroups_
    (
        this->coeffDict().template lookupOrDefault<label>
        (
            "staggerGroups", 1
        )
    ),
    staggerHashResolution_
    (
        this->coeffDict().template lookupOrDefault<scalar>
        (
            "staggerHashResolution", 1e-8
        )
    ),
    lastEventIndex_(-1),
    currentStaggerGroup_(0),
    activeCells_(),
    activeMassPerParcel_(),
    activeDiameters_(),
    globalCellStart_(0),
    globalActiveCells_(0),
    globalGroupCells_(0),
    volumeForEvent_(0),
    currentActiveCell_(-1),
    U0_(this->coeffDict().template get<vector>("U0")),
    useLocalU_
    (
        this->coeffDict().template lookupOrDefault<Switch>
        (
            "useLocalU", false
        )
    ),
    UFieldPtr_
    (
        useLocalU_
      ? &owner.db().objectRegistry::template lookupObject<volVectorField>
        (
            this->coeffDict().template lookupOrDefault<word>
            (
                "UField", "U"
            )
        )
      : nullptr
    ),
    useLocalT_
    (
        this->coeffDict().template lookupOrDefault<Switch>
        (
            "useLocalT", false
        )
    ),
    TFieldPtr_
    (
        useLocalT_
      ? &owner.db().objectRegistry::template lookupObject<volScalarField>
        (
            this->coeffDict().template lookupOrDefault<word>
            (
                "TField", "T"
            )
        )
      : nullptr
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    uInterp_(),
    tInterp_(),
    currentParcelU_(Zero),
    currentParcelT_(0)
{
    if
    (
        parcelsPerCell_ <= 0
     || injectionInterval_ <= 0
     || staggerGroups_ < 1
     || staggerHashResolution_ <= 0
     || rhoParticle_ <= 0
     || factor_ < 0 || factor_ > 1
    )
    {
        FatalErrorInFunction
            << "Require parcelsPerCell >= 1, injectionInterval > 0, "
            << "rhoParticle > 0, staggerGroups >= 1, "
            << "staggerHashResolution > 0 and 0 <= factor <= 1."
            << exit(FatalError);
    }

    Info<< "    periodic cell-based carbon injection" << nl
        << "    firstInjectionTime = " << firstInjectionTime_ << nl
        << "    injectionInterval = " << injectionInterval_ << nl
        << "    staggerGroups = " << staggerGroups_ << nl
        << "    group event spacing = "
        << injectionInterval_/scalar(staggerGroups_) << " s" << nl
        << "    staggerHashResolution = " << staggerHashResolution_
        << " m" << nl
        << "    minCarbonMassFraction = " << minCarbonMassFraction_ << nl
        << "    minCarbonMass = " << minCarbonMass_ << " kg/cell" << nl
        << "    parcelsPerCell = " << parcelsPerCell_ << nl
        << "    random position: volume-weighted cell tetrahedra" << endl;
}


template<class CloudType>
Foam::PeriodicCellCarbonInjection<CloudType>::PeriodicCellCarbonInjection
(
    const PeriodicCellCarbonInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    referenceField_(im.referenceField_),
    rhoField_(im.rhoField_),
    factor_(im.factor_),
    rhoParticle_(im.rhoParticle_),
    minCarbonMassFraction_(im.minCarbonMassFraction_),
    minCarbonMass_(im.minCarbonMass_),
    parcelsPerCell_(im.parcelsPerCell_),
    firstInjectionTime_(im.firstInjectionTime_),
    injectionInterval_(im.injectionInterval_),
    staggerGroups_(im.staggerGroups_),
    staggerHashResolution_(im.staggerHashResolution_),
    lastEventIndex_(im.lastEventIndex_),
    currentStaggerGroup_(im.currentStaggerGroup_),
    activeCells_(im.activeCells_),
    activeMassPerParcel_(im.activeMassPerParcel_),
    activeDiameters_(im.activeDiameters_),
    globalCellStart_(im.globalCellStart_),
    globalActiveCells_(im.globalActiveCells_),
    globalGroupCells_(im.globalGroupCells_),
    volumeForEvent_(im.volumeForEvent_),
    currentActiveCell_(im.currentActiveCell_),
    U0_(im.U0_),
    useLocalU_(im.useLocalU_),
    UFieldPtr_(im.UFieldPtr_),
    useLocalT_(im.useLocalT_),
    TFieldPtr_(im.TFieldPtr_),
    sizeDistribution_(im.sizeDistribution_.clone()),
    uInterp_(),
    tInterp_(),
    currentParcelU_(im.currentParcelU_),
    currentParcelT_(im.currentParcelT_)
{}


template<class CloudType>
Foam::PeriodicCellCarbonInjection<CloudType>::~PeriodicCellCarbonInjection()
{}


// * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PeriodicCellCarbonInjection<CloudType>::updateMesh()
{
    // Cell IDs are regenerated at every event. Clear stale lists when the
    // mesh changes. The carbon accounting field is mapped by the fvMesh.
    activeCells_.clear();
    activeMassPerParcel_.clear();
    activeDiameters_.clear();
    uInterp_.clear();
    tInterp_.clear();
    globalActiveCells_ = 0;
    globalGroupCells_ = 0;
    currentActiveCell_ = -1;
}


template<class CloudType>
Foam::scalar Foam::PeriodicCellCarbonInjection<CloudType>::timeEnd() const
{
    return GREAT;
}


template<class CloudType>
Foam::label Foam::PeriodicCellCarbonInjection<CloudType>::parcelsToInject
(
    const scalar,
    const scalar
)
{
    volumeForEvent_ = 0;

    const scalar time = this->owner().db().time().value();
    const scalar dt = this->owner().db().time().deltaTValue();
    const label eventIndex = staggeredCarbonEventIndex
    (
        time, dt, firstInjectionTime_, injectionInterval_, staggerGroups_
    );
    if (eventIndex == -2)
    {
        FatalErrorInFunction
            << "Time step " << dt << " s exceeds stagger group spacing "
            << injectionInterval_/scalar(staggerGroups_) << " s. "
            << "Reduce deltaT or staggerGroups."
            << exit(FatalError);
    }
    if (eventIndex < 0 || eventIndex == lastEventIndex_)
    {
        return 0;
    }

    lastEventIndex_ = eventIndex;
    currentStaggerGroup_ = eventIndex % staggerGroups_;
    updateInjectionData();

    if (globalActiveCells_ == 0 || volumeForEvent_ <= ROOTVSMALL)
    {
        return 0;
    }

    Info<< "StaggerEvent time=" << time
        << " subEvent=" << eventIndex
        << " group=" << currentStaggerGroup_
        << " cycle=" << eventIndex/staggerGroups_
        << " groups=" << staggerGroups_
        << " parcelsExpected=" << globalActiveCells_*parcelsPerCell_
        << " mass=" << volumeForEvent_*rhoParticle_ << " kg"
        << endl;

    return globalActiveCells_*parcelsPerCell_;
}


template<class CloudType>
Foam::scalar Foam::PeriodicCellCarbonInjection<CloudType>::volumeToInject
(
    const scalar,
    const scalar
)
{
    return volumeForEvent_;
}


template<class CloudType>
void Foam::PeriodicCellCarbonInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    const label activeI = localActiveIndex(parcelI);

    if (activeI < 0)
    {
        position = Zero;
        cellOwner = -1;
        tetFacei = -1;
        tetPti = -1;
        return;
    }

    const fvMesh& mesh = this->owner().mesh();
    const label celli = activeCells_[activeI];
    const List<tetIndices> cellTets =
        polyMeshTetDecomposition::cellTetIndices(mesh, celli);

    Random& rnd = this->owner().rndGen();
    scalar totalTetVolume = 0;
    forAll(cellTets, tetI)
    {
        totalTetVolume += cellTets[tetI].tet(mesh).mag();
    }

    if (cellTets.empty() || totalTetVolume <= VSMALL)
    {
        FatalErrorInFunction
            << "Cell " << celli << " has no valid tetrahedra"
            << exit(FatalError);
    }

    const scalar pick = rnd.sample01<scalar>()*totalTetVolume;
    scalar accumulated = 0;
    label selectedTet = cellTets.size() - 1;
    forAll(cellTets, tetI)
    {
        accumulated += cellTets[tetI].tet(mesh).mag();
        if (pick <= accumulated)
        {
            selectedTet = tetI;
            break;
        }
    }

    position = cellTets[selectedTet].tet(mesh).randomPoint(rnd);
    cellOwner = celli;
    tetFacei = cellTets[selectedTet].face();
    tetPti = cellTets[selectedTet].tetPt();
}


template<class CloudType>
void Foam::PeriodicCellCarbonInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    const label activeI = localActiveIndex(parcelI);
    currentActiveCell_ = activeI;

    if (activeI < 0)
    {
        return;
    }

    // Read gas velocity and temperature at the sampled parcel position.
    const auto& coordinates = parcel.coordinates();
    const auto& tetIs = parcel.currentTetIndices();

    parcel.U() = useLocalU_ && uInterp_.valid()
      ? uInterp_().interpolate(coordinates, tetIs) : U0_;
    if (useLocalT_ && tInterp_.valid())
    {
        setParcelTemperature
        (
            parcel,
            tInterp_().interpolate(coordinates, tetIs),
            0
        );
    }

    currentParcelU_ = parcel.U();
    currentParcelT_ = parcelTemperature(parcel, 0);
    parcel.d() = activeDiameters_[activeI];
}


template<class CloudType>
Foam::scalar Foam::PeriodicCellCarbonInjection<CloudType>::setNumberOfParticles
(
    const label,
    const scalar,
    const scalar diameter,
    const scalar rho
)
{
    const label i = currentActiveCell_;
    if (i < 0 || i >= activeCells_.size() || rho <= VSMALL)
    {
        return 0;
    }

    const scalar parcelMass = rho*pi/6.0*pow3(diameter);
    const scalar nP = activeMassPerParcel_[i]/(parcelMass + VSMALL);

    if (nP >= this->minParticlesPerParcel_)
    {
        volScalarField& massTransferred = const_cast<volScalarField&>
        (
            this->owner().mesh().template lookupObject<volScalarField>
            (
                "carbonInjectedMass"
            )
        );
        const scalar injectedMass = nP*parcelMass;
        const label celli = activeCells_[i];
        massTransferred[celli] += injectedMass;

        volScalarField& eTransferred = const_cast<volScalarField&>
        (
            this->owner().mesh().template lookupObject<volScalarField>
            (
                "carbonInjectedInternalEnergy"
            )
        );
        volVectorField& momentumTransferred = const_cast<volVectorField&>
        (
            this->owner().mesh().template lookupObject<volVectorField>
            (
                "carbonInjectedMomentum"
            )
        );
        const IOdictionary& transferDict =
            this->owner().mesh().template lookupObject<IOdictionary>
            (
                "carbonToParticleProperties"
            );
        const scalar CpParticle =
            transferDict.lookupOrDefault<scalar>("particleCp", 700.0);
        const scalar Tref =
            transferDict.lookupOrDefault<scalar>
            (
                "particleEnergyReferenceTemperature", 298.15
            );
        eTransferred[celli] +=
            injectedMass*CpParticle*(currentParcelT_ - Tref);
        momentumTransferred[celli] += injectedMass*currentParcelU_;
    }
    return nP;
}


template<class CloudType>
bool Foam::PeriodicCellCarbonInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::PeriodicCellCarbonInjection<CloudType>::validInjection
(
    const label parcelI
)
{
    return localActiveIndex(parcelI) >= 0;
}

// ************************************************************************* //
