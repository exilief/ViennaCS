#include <csDiffusion.hpp>

#include "geometry.hpp"
#include "parameters.hpp"
#include "solver.hpp"

namespace ls = viennals;

using T = double;
constexpr int D = 2;

const int substrateMaterial = 0;
const int maskMaterial = 1;
const int coverMaterial = 2;

template <typename Material> bool isMaterial(Material x, int material) {
  return static_cast<int>(x) == material;
}

template <class BoundaryFunc>
void addBoundaryValue(cs::DenseCellSet<T, D> &cellSet, BoundaryFunc isBoundary,
                      T boundaryValue, std::vector<T> &concentration) {
  auto &materials = *cellSet.getScalarData("Material");

#pragma omp parallel for
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isBoundary(i, cellSet, materials))
      concentration[i] = boundaryValue;
  }
}

int main(int argc, char **argv) {
  cs::Logger::setLogLevel(cs::LogLevel::INTERMEDIATE);

  Parameters<T> params;
  if (argc > 1) {
    auto config = cs::util::readFile(argv[1]);
    if (config.empty()) {
      std::cerr << "Empty config provided" << std::endl;
      return -1;
    }
    params.fromMap(config);
  }
  omp_set_num_threads(params.numThreads);

  auto matMap = cs::SmartPointer<ls::MaterialMap>::New();
  auto levelSets = geometry::makeStructure<T, D>(
      params, matMap, substrateMaterial, maskMaterial);

  cs::DenseCellSet<T, D> cellSet;
  T depth = params.substrateHeight + params.coverHeight + 10.;
  cellSet.setCellSetPosition(true); // isAboveSurface
  cellSet.setCoverMaterial(coverMaterial);
  cellSet.fromLevelSets(levelSets, matMap, depth);

  // We need neighborhood information for solving the diffusion equation
  cellSet.buildNeighborhood();

  // Boundary condition: constant concentration at the top
  auto isDirichletBoundary = [](unsigned cellIdx,
                                const cs::DenseCellSet<T, D> &cellSet,
                                const std::vector<T> &materials) {
    int topNeighbor = cellSet.getNeighbors(cellIdx)[2 * D - 1];
    return isMaterial(materials[cellIdx], substrateMaterial) &&
           topNeighbor != -1 &&
           isMaterial(materials[topNeighbor], coverMaterial);
  };

  auto concentration = cellSet.addScalarData("dopant", 0.);
  addBoundaryValue(cellSet, isDirichletBoundary, params.boundaryValue,
                   *concentration);

  cellSet.writeVTU("initial.vtu");

  if (params.velocity != 0.) {
    const T stability = 2 * params.diffusionCoefficient / params.velocity;
    std::cout << "Stability: " << stability << std::endl;
    if (0.5 * stability <= params.gridDelta)
      std::cout << "Unstable parameters. Reduce grid spacing!" << std::endl;
  }

  std::unordered_map<int, T> coeffs = {
      {substrateMaterial, params.diffusionCoefficient},
      {coverMaterial, 0},
      {maskMaterial, 0}};
  cs::ImplicitSolver solver(
      cellSet, cs::diffusion::MaterialCoeff(coeffs),
      cs::diffusion::makeBoundary<T, D>(isDirichletBoundary), "dopant",
      params.timeStabilityFactor);
  cs::Diffusion diffusion(cellSet, std::move(solver), params.duration);

  diffusion.apply();

  // saveVolumeMesh("final", levelSets, matMap);
  cellSet.writeVTU("final.vtu");
}
