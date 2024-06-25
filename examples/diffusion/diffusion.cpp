#include <csDenseCellSet.hpp>

#include "geometry.hpp"
#include "parameters.hpp"

namespace cs = viennacs;
namespace ls = viennals;

using T = double;
constexpr int D = 2;

const int substrateMaterial = 0;
const int maskMaterial = 1;
const int coverMaterial = 2;

template <typename Material> bool isMaterial(Material x, int material) {
  return static_cast<int>(x) == material;
}

void addConcentration(cs::DenseCellSet<T, D> &cellSet,
                      const Parameters<T> &params) {
  // Add quantity to be diffused (on top of the cell material)
  auto concentration = cellSet.addScalarData("dopant", 0.);
  auto materials = cellSet.getScalarData("Material");

  // Boundary condition: constant concentration at the top (outside)
#pragma omp parallel for
  for (int i = 0; i < cellSet.getNumberOfCells(); ++i) {
    if (isMaterial((*materials)[i], coverMaterial) &&
        cellSet.getCellCenter(i)[D - 1] <
            params.substrateHeight + params.gridDelta) {
      (*concentration)[i] = params.boundaryValue;
    }
  }
}

// Explicit diffusion time-step (forward Euler)
void solveDiffusionStep(cs::DenseCellSet<T, D> &cellSet,
                        const Parameters<T> &params, T dt) {
  auto data = cellSet.getScalarData("dopant");
  auto materials = cellSet.getScalarData("Material");
  std::vector<T> solution(data->size(), 0.);
  const T dx = params.gridDelta;
  const T C = dt * params.diffusionCoefficient / (dx * dx);

#pragma omp parallel for
  for (int e = 0; e < data->size(); e++) {
    if (!isMaterial((*materials)[e], substrateMaterial)) {
      solution[e] = (*data)[e];
      continue;
    }

    auto coord = cellSet.getCellCenter(e);
    int numNeighbors = 0;

    const auto &cellNeighbors = cellSet.getNeighbors(e);
    for (auto n : cellNeighbors) {
      if (n == -1 || isMaterial((*materials)[n], maskMaterial))
        continue;

      solution[e] += (*data)[n];
      numNeighbors++;
    }

    // Diffusion
    solution[e] = (*data)[e] + C * (solution[e] - numNeighbors * (*data)[e]);
  }
  *data = std::move(solution);
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

  addConcentration(cellSet, params);
  cellSet.writeVTU("initial.vtu");

  if (params.velocity != 0.) {
    const T stability = 2 * params.diffusionCoefficient / params.velocity;
    std::cout << "Stability: " << stability << std::endl;
    if (0.5 * stability <= params.gridDelta)
      std::cout << "Unstable parameters. Reduce grid spacing!" << std::endl;
  }

  T duration = params.duration;
  T dx = params.gridDelta;
  T dt = std::min(dx * dx / (params.diffusionCoefficient * 2 * D) *
                      params.timeStabilityFactor,
                  duration);
  T time = 0.;
  while (time < duration) {

    if (time + dt > duration)
      dt = duration - time;

    solveDiffusionStep(cellSet, params, dt);

    time += dt;
  }

  // saveVolumeMesh("final", levelSets, matMap);
  cellSet.writeVTU("final.vtu");
}
