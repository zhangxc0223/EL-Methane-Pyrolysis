/*---------------------------------------------------------------------------*\
  Shared injection/transfer schedule and decomposition-independent cell groups.

  Each cell is assigned a permanent phase group using its geometric centre.
  Group g injects at firstTime + (cycle*groups + g)*period/groups.
  Every cell therefore keeps the original period while simultaneous global
  mass sinks are split over the period.
\*---------------------------------------------------------------------------*/

#ifndef StaggeredCarbonSchedule_H
#define StaggeredCarbonSchedule_H

#include "label.H"
#include "scalar.H"
#include "vector.H"
#include <cmath>
#include <cstdint>

namespace Foam
{

// Return the group sub-event occurring in the current physical time step.
// -1 means none. The strict lower bound avoids replaying an event after restart.
// Requires dt <= period/groups, so at most one group fires per time step.
inline label staggeredCarbonEventIndex
(
    const scalar now,
    const scalar dt,
    const scalar firstTime,
    const scalar period,
    const label groups
)
{
    if (groups < 1 || period <= 0 || dt <= 0)
    {
        return -1;
    }

    const scalar subPeriod = period/scalar(groups);
    if (dt > subPeriod*(scalar(1) + scalar(1e-7)))
    {
        return -2;  // Unsupported: multiple events in a single time step.
    }

    const scalar tolerance = scalar(1e-8)*dt;
    if (now < firstTime - tolerance)
    {
        return -1;
    }

    const label index = label
    (
        std::floor((now - firstTime + tolerance)/subPeriod)
    );
    if (index < 0)
    {
        return -1;
    }

    const scalar scheduledTime = firstTime + scalar(index)*subPeriod;
    if
    (
        scheduledTime <= now + tolerance
     && scheduledTime > now - dt + tolerance
    )
    {
        return index;
    }
    return -1;
}


// Deterministic spatial hash. No processor-local cell ID is used.
// The resolution is in metres; using the same mesh yields stable groups
// across decompositions. This is not intended for topology-changing meshes.
inline label staggeredCarbonCellGroup
(
    const vector& centre,
    const label groups,
    const scalar resolution
)
{
    if (groups <= 1)
    {
        return 0;
    }

    const std::int64_t x = static_cast<std::int64_t>
    (
        std::llround(centre.x()/resolution)
    );
    const std::int64_t y = static_cast<std::int64_t>
    (
        std::llround(centre.y()/resolution)
    );
    const std::int64_t z = static_cast<std::int64_t>
    (
        std::llround(centre.z()/resolution)
    );

    // FNV-1a and an avalanche finaliser; stable across rank partitions.
    std::uint64_t h = UINT64_C(14695981039346656037);
    h = (h ^ static_cast<std::uint64_t>(x))*UINT64_C(1099511628211);
    h = (h ^ static_cast<std::uint64_t>(y))*UINT64_C(1099511628211);
    h = (h ^ static_cast<std::uint64_t>(z))*UINT64_C(1099511628211);
    h ^= h >> 33;
    h *= UINT64_C(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h *= UINT64_C(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    return label(h % static_cast<std::uint64_t>(groups));
}

} // namespace Foam

#endif
