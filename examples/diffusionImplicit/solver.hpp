#pragma once

#include <csDenseCellSet.hpp>

#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseLU>

namespace cs {

using namespace viennacs;

template <class T, int D, class Coeff, class BC> class ImplicitSolver {
  DenseCellSet<T, D> *cellSet = nullptr;

  Coeff coeff;
  BC bc;
  std::vector<T> *data;
  T timeStabilityFactor;
  T dtCached = 0;

  struct System {
    // unique_ptr to make it moveable (private copy ctor)
    std::unique_ptr<Eigen::SparseLU<Eigen::SparseMatrix<T>>> solver =
        std::make_unique<typename decltype(solver)::element_type>();
    Eigen::SparseMatrix<T> systemMatrix;
    Eigen::Matrix<T, Eigen::Dynamic, 1> rhs;
    std::unordered_map<int, int> cellMapping; // cellSet -> solution index
    unsigned numCells = 0;
  };

  System sys;

public:
  ImplicitSolver(DenseCellSet<T, D> &cellSet, Coeff coeff, BC bc,
                 const std::string &dataId, T timeStabilityFactor = 1.)
      : cellSet(&cellSet), coeff(std::move(coeff)), bc(std::move(bc)),
        data(cellSet.getScalarData(dataId)),
        timeStabilityFactor(timeStabilityFactor) {
    initCellMapping();
    assembleRHS();
  }

  // Implicit diffusion time-step (backward Euler)
  void step(T dt) {
    setTimeStep(dt);

    // TODO: Add volumetric sources here: rhs += sources*dt
    sys.rhs = sys.solver->solve(sys.rhs);

    if (sys.solver->info() != Eigen::Success) {
      Logger::getInstance().addError("(Eigen) Matrix solving failed.").print();
    }

    // Write results to cellSet
    for (const auto &i : sys.cellMapping) {
      (*data)[i.first] = sys.rhs[i.second];
    }
  }

  void setTimeStep(T dt) {
    if (dt != dtCached) {
      assembleMatrix(dt);
      sys.solver->compute(sys.systemMatrix);
      dtCached = dt;
    }
  }

  T stableTimeStep(T dtRequested = 0) {
    // The implicit scheme is unconditionally stable
    if (dtRequested)
      return dtRequested;

    // Use the stable limit for the explicit scheme as default
    T dx = cellSet->getGridDelta();
    T dt = dx * dx / (coeff.max() * 2 * D) * timeStabilityFactor;
    return dt;
  }

private:
  void assembleMatrix(T dt) {
    auto &materials = *cellSet->getScalarData("Material");
    const T dx = cellSet->getGridDelta();
    const T dtdx2 = dt / (dx * dx);

    std::vector<Eigen::Triplet<T>> triplets; // (i,j,value) matrix entries
    triplets.reserve(2 * D * sys.numCells);

    for (const auto &ids : sys.cellMapping) {
      auto i = ids.first;
      auto cellMapIdx = ids.second;
      bool onBoundary = bc.isDirichletBoundary(i, *cellSet, materials);

      // Boundary conditions
      if (onBoundary) {
        triplets.push_back({cellMapIdx, cellMapIdx, 1.}); // diagonal
        continue;
      }

      int mat = static_cast<int>(materials[i]);

      const auto &neighbors = cellSet->getNeighbors(i);
      T D_sum = 0;
      for (auto n : neighbors) {
        if (n == -1)
          continue;
        int matN = static_cast<int>(materials[n]);
        if (coeff(matN) != 0) {
          T c = coeff.meanValue(mat, matN);
          D_sum += c;
          triplets.emplace_back(cellMapIdx, sys.cellMapping[n], -dtdx2 * c);
        }
      }

      triplets.emplace_back(cellMapIdx, cellMapIdx, 1. + dtdx2 * D_sum);
    }

    sys.systemMatrix.resize(sys.numCells, sys.numCells);
    sys.systemMatrix.setFromTriplets(triplets.begin(), triplets.end());
    sys.systemMatrix.makeCompressed();
  }

  void assembleRHS() {
    sys.rhs.resize(sys.numCells);

    // Add initial values (+ Dirichlet boundary conditions)
    for (const auto &ids : sys.cellMapping) {
      sys.rhs[ids.second] = (*data)[ids.first];
    }
  }

  void initCellMapping() {
    auto &materials = *cellSet->getScalarData("Material");
    sys.numCells = 0;
    for (int i = 0; i < cellSet->getNumberOfCells(); ++i) {
      if (coeff(materials[i]) != 0) {
        sys.cellMapping[i] = sys.numCells++;
      }
    }
  }
};

} // namespace cs
