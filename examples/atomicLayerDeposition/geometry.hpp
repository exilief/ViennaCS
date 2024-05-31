#pragma once

#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>

#include <vcUtil.hpp>

template <class NumericType, int D>
auto makeLShape(viennacore::util::Parameters &params) {
  static_assert(D == 2, "This function only works in 2D");

  using namespace viennals;

  std::vector<SmartPointer<Domain<NumericType, D>>> levelSets;

  const auto gridDelta = params.get("gridDelta");

  double bounds[2 * D];
  bounds[0] = -params.get("verticalWidth") / 2. - params.get("xPad");
  bounds[1] = -params.get("verticalWidth") / 2. + params.get("xPad") +
              params.get("horizontalWidth");
  bounds[2] = -gridDelta;
  bounds[3] = params.get("verticalDepth") + gridDelta;

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
  origin[D - 1] = params.get("verticalDepth");
  MakeGeometry<NumericType, D>(
      substrate, SmartPointer<Plane<NumericType, D>>::New(origin, normal))
      .apply();

  {
    auto vertBox =
        SmartPointer<Domain<NumericType, D>>::New(substrate->getGrid());
    NumericType minPoint[D] = {-params.get("verticalWidth") / 2.0, 0.};
    NumericType maxPoint[D] = {params.get("verticalWidth") / 2.0,
                               params.get("verticalDepth")};
    MakeGeometry<NumericType, D>(
        vertBox, SmartPointer<Box<NumericType, D>>::New(minPoint, maxPoint))
        .apply();

    BooleanOperation<NumericType, D>(substrate, vertBox,
                                     BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    auto horiBox =
        SmartPointer<Domain<NumericType, D>>::New(substrate->getGrid());
    NumericType minPoint[D] = {-params.get("verticalWidth") / 2.0, 0.};
    NumericType maxPoint[D] = {-params.get("verticalWidth") / 2.0 +
                                   params.get("horizontalWidth"),
                               params.get("horizontalHeight")};

    MakeGeometry<NumericType, D>(
        horiBox, SmartPointer<Box<NumericType, D>>::New(minPoint, maxPoint))
        .apply();

    BooleanOperation<NumericType, D>(substrate, horiBox,
                                     BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  levelSets.push_back(substrate);
  return levelSets;
}
