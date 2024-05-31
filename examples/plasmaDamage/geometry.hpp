#pragma once

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>

#include <vcUtil.hpp>

template <class NumericType, int D>
auto makeFin(viennacore::util::Parameters &params) {
  static_assert(D == 3, "This function only works in 3D");

  using namespace viennals;

  std::vector<SmartPointer<Domain<NumericType, D>>> levelSets;

  const NumericType gridDelta = params.get("gridDelta");

  double bounds[2 * D];
  bounds[0] = -params.get("xExtent") / 2.;
  bounds[1] = params.get("xExtent") / 2.;
  bounds[2] = -params.get("yExtent") / 2.;
  bounds[3] = params.get("yExtent") / 2.;
  bounds[4] = 0.;
  bounds[5] = 1.;

  typename Domain<NumericType, D>::BoundaryType boundaryCons[D];

  for (int i = 0; i < D - 1; i++) {
    boundaryCons[i] = Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = SmartPointer<Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);
  NumericType normal[D] = {0.};
  NumericType origin[D] = {0.};
  normal[D - 1] = 1.;
  MakeGeometry<NumericType, D>(
      substrate, SmartPointer<Plane<NumericType, D>>::New(origin, normal))
      .apply();

  {
    auto fin = SmartPointer<Domain<NumericType, D>>::New(substrate->getGrid());
    NumericType minPoint[D] = {-params.get("width") / 2.0,
                               -params.get("length") / 2.0, 0.};
    NumericType maxPoint[D] = {params.get("width") / 2.0,
                               params.get("length") / 2.0,
                               params.get("height")};
    MakeGeometry<NumericType, D>(
        fin, SmartPointer<Box<NumericType, D>>::New(minPoint, maxPoint))
        .apply();

    BooleanOperation<NumericType, D>(substrate, fin,
                                     BooleanOperationEnum::UNION)
        .apply();
  }

  levelSets.push_back(substrate);
  return levelSets;
}
