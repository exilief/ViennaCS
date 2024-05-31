#include <vcUtil.hpp>

#include <csAtomicLayerProcess.hpp>
#include <csMeanFreePath.hpp>
#include <csSegmentCells.hpp>

#include "geometry.hpp"

namespace cs = viennacs;

int main(int argc, char *argv[]) {
  constexpr int D = 2;
  using NumericType = double;

  // Parse the parameters
  cs::util::Parameters params;
  if (argc > 1) {
    params.readConfigFile(argv[1]);
  } else {
    std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
    return 1;
  }
  omp_set_num_threads(params.get<int>("numThreads"));

  // Create a cellSet
  auto cellSet = cs::SmartPointer<cs::DenseCellSet<NumericType, D>>::New();
  cellSet->setCellSetPosition(true);

  auto levelSets = makeLShape<NumericType, D>(params);
  // Generate the cell set from the domain
  cellSet->fromLevelSets(levelSets, nullptr,
                         params.get("verticalDepth") + params.get("topSpace"));
  cellSet->writeVTU("initial.vtu");

  // Segment the cells into surface, material, and gas cells
  cs::SegmentCells<NumericType, D> segmentation(cellSet);
  segmentation.setBulkMaterial(1);
  segmentation.apply();

  cs::Timer timer;
  timer.start();

  // Calculate the mean free path for the gas cells
  cs::MeanFreePath<NumericType, D> mfpCalc(cellSet);
  mfpCalc.setNumRaysPerCell(params.get("raysPerCell"));
  mfpCalc.setReflectionLimit(params.get<int>("reflectionLimit"));
  mfpCalc.setRngSeed(params.get<int>("seed"));
  mfpCalc.setMaterial(1);
  mfpCalc.setBulkLambda(params.get("bulkLambda"));
  mfpCalc.apply();

  timer.finish();
  std::cout << "Mean free path calculation took " << timer.totalDuration * 1e-9
            << " seconds." << std::endl;

  cellSet->writeVTU("initial.vtu");

  cs::AtomicLayerProcess<NumericType, D> model(cellSet);
  model.setMaxLambda(params.get("bulkLambda"));
  model.setPrintInterval(params.get("printInterval"));
  model.setStabilityFactor(params.get("stabilityFactor"));
  model.setFirstPrecursor("H2O", params.get("H2O_meanThermalVelocity"),
                          params.get("H2O_adsorptionRate"),
                          params.get("H2O_desorptionRate"),
                          params.get("p1_time"), params.get("inFlux"));
  model.setSecondPrecursor("TMA", params.get("TMA_meanThermalVelocity"),
                           params.get("TMA_adsorptionRate"),
                           params.get("TMA_desorptionRate"),
                           params.get("p2_time"), params.get("inFlux"));
  model.setPurgeParameters(params.get("purge_meanThermalVelocity"),
                           params.get("purge_time"));
  // The deposition probability is (H2O_cov * TMA_cov)^order
  model.setReactionOrder(params.get("reactionOrder"));

  model.apply();

  cellSet->writeVTU("final.vtu");
}
