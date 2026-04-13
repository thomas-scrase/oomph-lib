#ifndef GENERALISED_FLUID_CONSTITUTIVE_LAW_HEADER
#define GENERALISED_FLUID_CONSTITUTIVE_LAW_HEADER

#include "generic.h"

namespace oomph
{
  class GeneralisedFluidConstitutiveLaw
  {
  public:
    GeneralisedFluidConstitutiveLaw() {}
    ~GeneralisedFluidConstitutiveLaw() {}

    void get_equation_of_state(const double& rho,
                               const double& p,
                               const double& e,
                               double& residual) const
    {
      residual = p - e * rho;
    }

    void get_deviatoric_stress(const DenseMatrix<double>& dudx,
                               DenseMatrix<double>& sigma) const
    {
      const unsigned dim = dudx.nrow();

      DenseMatrix<double> D(dim, dim, 0.0);
      double II = 0.0;
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          D(i, j) = 0.5 * (dudx(i, j) + dudx(j, i)); // D=(dudx+dudx^T) / 2
          II += 2.0 * D(i, j) * D(i, j); // II=2 D:D
        }
      }
      double gamma = sqrt(fmax(II, 1e-30)); // regularize small values

      const double eta_inf = 1.0;
      const double eta0 = 0.0;
      const double lambda = 1.0;
      const double n = 1.0;

      const double eta =
        eta0 +
        (eta_inf - eta0) * pow(1.0 + pow(lambda * gamma, 2.0), (n - 1.0) / 2.0);

      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          sigma(i, j) = eta * D(i, j);
        }
      }

      // We always need the deviatoric stress tensor
      double tr = 0.0;
      for (unsigned i = 0; i < dim; i++)
      {
        tr += sigma(i, i);
      }

      for (unsigned i = 0; i < dim; i++)
      {
        sigma(i, i) -= tr / 3.0;
      }
    }
  };
} // namespace oomph

#endif