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

Class
    Foam::PeriodicCellCarbonInjection

Description
    Periodically inject parcels from every local mesh cell that passes a
    carbon mass-fraction and/or per-cell carbon-mass threshold.

    The injection position is independently randomised for every parcel and
    event, using volume-weighted tetrahedral sampling in the actual cell.
    No positionsFile or precomputed cell centres are required.

    Each selected cell transfers factor*Y_C*rho*V of carbon.
    The mass is split between the specified number of parcels.

    The model records injected carbon mass for removal from the gas phase
    during the same time step.

SourceFiles
    PeriodicCellCarbonInjection.C

\*---------------------------------------------------------------------------*/

#ifndef PeriodicCellCarbonInjection_H
#define PeriodicCellCarbonInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "volFieldsFwd.H"
#include "Switch.H"
#include "interpolation.H"

namespace Foam
{

template<class CloudType>
class PeriodicCellCarbonInjection
:
    public InjectionModel<CloudType>
{
    // Model coefficients
    const volScalarField& referenceField_;
    const volScalarField& rhoField_;
    const scalar factor_;
    const scalar rhoParticle_;
    const scalar minCarbonMassFraction_;
    const scalar minCarbonMass_;
    const label parcelsPerCell_;

    // Injection times, in seconds.
    const scalar firstInjectionTime_;
    const scalar injectionInterval_;
    label lastEventIndex_;

    // Cells with enough carbon to inject.
    labelList activeCells_;
    scalarList activeMassPerParcel_;
    scalarList activeDiameters_;
    label globalCellStart_;
    label globalActiveCells_;
    scalar volumeForEvent_;
    mutable label currentActiveCell_;

    const vector U0_;
    const Switch useLocalU_;
    const volVectorField* UFieldPtr_;
    const Switch useLocalT_;
    const volScalarField* TFieldPtr_;
    const autoPtr<distributionModel> sizeDistribution_;

    // Interpolate gas values at the new particle positions.
    autoPtr<interpolation<vector>> uInterp_;
    autoPtr<interpolation<scalar>> tInterp_;
    vector currentParcelU_;
    scalar currentParcelT_;

    // Find the local cell for a parcel index.
    inline label localActiveIndex(const label parcelI) const;

    // Update the list of cells and injection masses.
    void updateInjectionData();

    template<class ParcelType>
    static auto setParcelTemperature
    (
        ParcelType& parcel,
        const scalar T,
        const int
    ) -> decltype(parcel.T() = T, void());

    template<class ParcelType>
    static void setParcelTemperature
    (
        ParcelType& parcel,
        const scalar T,
        const long
    );

    // Kinematic parcels do not have a temperature field.
    template<class ParcelType>
    static auto parcelTemperature
    (
        const ParcelType& parcel, const int
    ) -> decltype(parcel.T(), scalar());

    template<class ParcelType>
    static scalar parcelTemperature
    (
        const ParcelType& parcel, const long
    );

protected:
    virtual scalar setNumberOfParticles
    (
        const label parcels,
        const scalar volumeFraction,
        const scalar diameter,
        const scalar rho
    );

public:
    TypeName("periodicCellCarbonInjection");

    PeriodicCellCarbonInjection
    (
        const dictionary& dict,
        CloudType& owner,
        const word& modelName
    );

    PeriodicCellCarbonInjection(const PeriodicCellCarbonInjection<CloudType>& im);

    virtual autoPtr<InjectionModel<CloudType>> clone() const
    {
        return autoPtr<InjectionModel<CloudType>>
        (
            new PeriodicCellCarbonInjection<CloudType>(*this)
        );
    }

    virtual ~PeriodicCellCarbonInjection();

    virtual void updateMesh();
    virtual scalar timeEnd() const;
    virtual label parcelsToInject(const scalar time0, const scalar time1);
    virtual scalar volumeToInject(const scalar time0, const scalar time1);

    virtual void setPositionAndCell
    (
        const label parcelI,
        const label nParcels,
        const scalar time,
        vector& position,
        label& cellOwner,
        label& tetFacei,
        label& tetPti
    );

    virtual void setProperties
    (
        const label parcelI,
        const label nParcels,
        const scalar time,
        typename CloudType::parcelType& parcel
    );

    virtual bool fullyDescribed() const;
    virtual bool validInjection(const label parcelI);
};

} // End namespace Foam

#ifdef NoRepository
    #include "PeriodicCellCarbonInjection.C"
#endif

#endif
// ************************************************************************* //
