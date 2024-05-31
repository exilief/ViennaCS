#pragma once

#include <csTracingParticle.hpp>

template <class T, int D>
class PlasmaDamageIon : public viennacs::Particle<PlasmaDamageIon<T, D>, T> {
public:
  PlasmaDamageIon(const T passedMeanEnergy = 100.,
                  const T passedMeanFreePath = 1.)
      : meanIonEnergy(passedMeanEnergy), meanFreePath(passedMeanFreePath) {}

  void initNew(viennacs::RNG &rngState) override final {
    std::uniform_real_distribution<T> uniDist;
    do {
      const auto rand1 = uniDist(rngState);
      const auto rand2 = uniDist(rngState);
      E = std::cos(M_PI * 2 * rand1) * std::sqrt(-2. * std::log(rand2)) *
              deltaIonEnergy +
          meanIonEnergy;
    } while (E < minEnergy);
  }

  std::pair<T, viennacs::Vec3D<T>>
  surfaceHit(const viennacs::Vec3D<T> &rayDir,
             const viennacs::Vec3D<T> &geomNormal, bool &reflect,
             viennacs::RNG &rngState) override final {
    auto cosTheta = -viennacs::DotProduct(rayDir, geomNormal);
    const T incAngle = std::acos(std::max(std::min(cosTheta, T(1)), T(0)));
    std::uniform_real_distribution<T> uniDist;

    T Eref_peak = 0;

    // Small incident angles are reflected with the energy fraction centered at
    // 0
    if (incAngle >= inflectAngle) {
      Eref_peak =
          Eref_max *
          (1 - (1 - A) * std::pow((M_PI_2 - incAngle) / (M_PI_2 - inflectAngle),
                                  n_r));
    } else {
      Eref_peak = Eref_max * A * std::pow(incAngle / inflectAngle, n_l);
    }
    // Gaussian distribution around the Eref_peak scaled by the particle energy
    T tempEnergy = Eref_peak * E;

    T NewEnergy;
    do {
      const auto rand1 = uniDist(rngState);
      const auto rand2 = uniDist(rngState);
      NewEnergy = tempEnergy +
                  (std::min((E - tempEnergy), tempEnergy) + E * 0.05) *
                      (1 - 2. * rand1) * std::sqrt(std::fabs(std::log(rand2)));
    } while (NewEnergy > E || NewEnergy <= 0.);

    auto impactEnergy = E - NewEnergy;

    if (NewEnergy > minEnergy) {
      reflect = true;
      auto direction = viennaray::ReflectionConedCosine<T, D>(
          rayDir, geomNormal, rngState, std::min(incAngle, minAngle));
      E = NewEnergy;
      return std::pair<T, viennacs::Vec3D<T>>{impactEnergy, direction};
    } else {
      reflect = false;
      return std::pair<T, viennacs::Vec3D<T>>{impactEnergy,
                                              viennacs::Vec3D<T>{0., 0., 0.}};
    }
  }

  T collision(
      viennacs::VolumeParticle<T> &particle, viennacs::RNG &rngState,
      std::vector<viennacs::VolumeParticle<T>> &particleStack) override final {
    using namespace viennacs;

    T fill = 0.;

    // inelastic losses through electrons (stopping power)
    particle.energy *=
        std::pow(nonLocalLosses, particle.distance * scaleFactor);

    if (particle.energy < displacementEnergyThreshold) {
      particle.energy = -1;
      return fill;
    }

    int numParticles = 0;
    if (particle.scattered < maxScatter) {
      std::uniform_int_distribution<> particleDist(1, maxScatter -
                                                          particle.scattered);
      numParticles = particleDist(rngState);
    }

    for (int i = 0; i < numParticles; i++) {
      T cosTheta, tmp, sinThetaSqr;
      Vec3D<T> direction;
      do {
        // random direction
        direction[0] = negUniDist(rngState);
        direction[1] = negUniDist(rngState);
        direction[2] = negUniDist(rngState);

        // normalize
        tmp = Norm(direction);
        direction = direction * (T(1) / tmp);

        // cos(angle)
        cosTheta = DotProduct(particle.direction, direction);
        // sin(angle)^2
        sinThetaSqr = 1 - cosTheta * cosTheta;
      } while (sinThetaSqr < mu * mu);

      // energy of scattered ion
      tmp = particle.energy *
            (cosTheta + std::sqrt(sinThetaSqr - mu * mu) * pre_fac) *
            (cosTheta + std::sqrt(sinThetaSqr - mu * mu) * pre_fac);

      // create new ion
      if (tmp > displacementEnergyThreshold) {
        particleStack.emplace_back(
            VolumeParticle<T>{particle.position, direction, tmp, 0.,
                              particle.cellId, particle.scattered - 1});
      }

      // energy transferred to atom (damage)
      tmp = particle.energy - tmp;
      if (tmp >= displacementEnergyThreshold) {
        // Kinchin-Pease damage model
        if (tmp < 2 * displacementEnergyThreshold / 0.8)
          fill += 1;
        else
          fill += 0.8 * tmp / (2 * displacementEnergyThreshold);
      }
    }

    particle.energy = -1; // kill current particle

    return fill;
  }

  T getSourceDistributionPower() const override final { return 1000.; }
  std::array<T, 2> getMeanFreePath() const override final {
    return {meanFreePath, meanFreePath / T(2)};
  }

private:
  static constexpr T A_p = 0.0337;
  static constexpr T A_Si = 7.;
  static constexpr T A_O = 2;
  static constexpr T A_SiO2 = 0.3;

  static constexpr T sqrt_Eth_p = 0.;
  static constexpr T sqrt_Eth_Si = 3.8729833462;
  static constexpr T sqrt_Eth_O = 3.8729833462;
  static constexpr T sqrt_Eth_SiO2 = 3.8729833462;
  static constexpr T Eref_max = 1.;

  const T meanFreePath = 1.;
  const T meanIonEnergy = 100.;
  const T deltaIonEnergy = meanIonEnergy / 10.;
  static constexpr T minEnergy = 1.; // Discard particles with energy < 1eV

  static constexpr T inflectAngle = 1.55334;
  static constexpr T minAngle = 1.3962634;
  static constexpr T n_l = 10.;
  static constexpr T n_r = 1.;

  static constexpr T A = 1. / (1. + (n_l / n_r) * (M_PI_2 / inflectAngle - 1.));

  T E;

  std::uniform_real_distribution<T> negUniDist =
      std::uniform_real_distribution<T>(-1., 1.);

  const T nonLocalLosses = 0.9;
  const int maxScatter = 10;
  static constexpr T scaleFactor = 1.;
  static constexpr T displacementEnergyThreshold = 15;
  static constexpr T mu = 28.0855 / 39.948;
  static constexpr T pre_fac = (1. / (1. + mu));
};