#pragma once

#include <csDenseCellSet.hpp>

namespace viennacs {

using namespace viennacore;

template <class T, int D, class Solver> class Diffusion {
  DenseCellSet<T, D> *cellSet;

  Solver solver;
  T duration = 0;

public:
  Diffusion(DenseCellSet<T, D> &cellSet, Solver solver, T duration)
      : cellSet(&cellSet), solver(std::move(solver)), duration(duration) {}

  void apply() {
    T dt = solver.stableTimeStep();
    T time = 0;

    while (time < duration) {
      if (time + dt > duration)
        dt = duration - time;

      solver.step(dt);

      time += dt;
    }
  }
};

namespace diffusion {

template <class T> T harmonicMean(T c1, T c2) {
  assert(c1 && c2);

  return 2. / (1. / c1 + 1. / c2);
}

template <class T, int D, class Coeff, class BC> class ExplicitSolver {
  DenseCellSet<T, D> *cellSet = nullptr;

  Coeff coeff;
  BC bc;
  std::vector<T> *data;
  T timeStabilityFactor;

public:
  ExplicitSolver(DenseCellSet<T, D> &cellSet, Coeff coeff, BC bc,
                 const std::string &dataId, T timeStabilityFactor = 0.95)
      : cellSet(&cellSet), coeff(std::move(coeff)), bc(std::move(bc)),
        data(cellSet.getScalarData(dataId)),
        timeStabilityFactor(timeStabilityFactor) {}

  // Explicit diffusion time-step (forward Euler)
  void step(T dt) {
    auto &materials = *cellSet->getScalarData("Material");
    std::vector<T> solution(data->size(), 0.);
    const T dx = cellSet->getGridDelta();
    const T dt_dx2 = dt / (dx * dx);

#pragma omp parallel for
    for (int i = 0; i < data->size(); i++) {
      int mat = static_cast<int>(materials[i]);
      bool onBoundary = bc.isDirichletBoundary(i, *cellSet, materials);

      if (onBoundary || coeff(mat) == 0) {
        solution[i] = (*data)[i];
        continue;
      }

      T D_sum = 0;

      for (auto n : cellSet->getNeighbors(i)) {
        if (n == -1)
          continue;

        int matN = static_cast<int>(materials[n]);
        if (coeff(matN) == 0)
          continue;

        T c = coeff.meanValue(mat, matN);

        solution[i] += c * (*data)[n];
        D_sum += c;
      }

      solution[i] = (*data)[i] * (1. - dt_dx2 * D_sum) + dt_dx2 * solution[i];
    }
    *data = std::move(solution);
  }

  void setTimeStep(T dt) {}

  T stableTimeStep(T dtRequested = 0) {
    T dx = cellSet->getGridDelta();
    T dtLimit = dx * dx / (coeff.max() * 2 * D) * timeStabilityFactor;
    if (dtRequested)
      return std::min(dtRequested, dtLimit);
    return dtLimit;
  }
};

template <class T> class ConstCoeff {
  T C;

public:
  ConstCoeff(T C) : C(C) {}

  T operator()(int) const { return C; }

  T operator()(int, int) const { return C; }

  T meanValue(int, int) const { return C; }

  T max() const { return C; }
};

template <class T> class MaterialCoeff {
  std::unordered_map<int, T> C;

public:
  MaterialCoeff(std::unordered_map<int, T> C) : C(std::move(C)) {}

  T operator()(int material) { return C[material]; }

  T operator()(int mat1, int mat2) {
    return C[mat1] == 0 || C[mat2] == 0 ? 0 : meanValue(mat1, mat2);
  }

  // Assume non-zero
  T meanValue(int mat1, int mat2) const {
    return mat1 == mat2 ? C.at(mat1) : harmonicMean(C.at(mat1), C.at(mat2));
  }

  T max() const {
    T val = 0;
    for (const auto &pair : C)
      val = std::max(val, pair.second);
    return val;
  }
};

template <class T, int D,
          class DirichletBC = bool (*)(unsigned, const DenseCellSet<T, D> &,
                                       const std::vector<T> &)>
class Boundary {
  DirichletBC dirichletBC = +[](unsigned, const DenseCellSet<T, D> &,
                                const std::vector<T> &) { return false; };

public:
  Boundary() {}

  Boundary(DirichletBC dirichletBC) : dirichletBC(std::move(dirichletBC)) {}

  bool isDirichletBoundary(unsigned cellIdx, const DenseCellSet<T, D> &cellSet,
                           const std::vector<T> &materials) {
    return dirichletBC(cellIdx, cellSet, materials);
  }
};

template <class T, int D, class DirichletBC>
auto makeBoundary(DirichletBC dirichletBC) {
  return Boundary<T, D, DirichletBC>(dirichletBC);
}

} // namespace diffusion

} // namespace viennacs
