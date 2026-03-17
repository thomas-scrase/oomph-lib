// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Header file for AnisotropicConstitutiveLaw objects

#ifndef OOMPH_ANISOTROPIC_CONSTITUTIVE_LAWS_HEADER
#define OOMPH_ANISOTROPIC_CONSTITUTIVE_LAWS_HEADER

// Config header
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "constitutive_laws.h"


namespace oomph
{
  class FibreReinforcedMooneyRivlin : public MooneyRivlin
  {
  public:
    FibreReinforcedMooneyRivlin(double* c1_pt, double* c2_pt, double* c3_pt)
      : MooneyRivlin(c1_pt, c2_pt), C3_pt(c3_pt)
    {
    }

    virtual ~FibreReinforcedMooneyRivlin() {}

    double W(const Vector<double>& I) override
    {
      return MooneyRivlin::W(I) + (*C3_pt) * I[3];
    }

    void derivatives(Vector<double>& I, Vector<double>& dWdI)
    {
      MooneyRivlin::derivatives(I, dWdI);
      dWdI[3] = (*C3_pt);
    }

    void get_I_anisotropic(const DenseMatrix<double>& g,
                           const DenseMatrix<double>& G,
                           const Vector<double>& fields,
                           Vector<double>& I,
                           Vector<DenseMatrix<double>>& dIdG) const
    {
      const unsigned dim = g.nrow();
#ifdef RANGE_CHECKING
      if (fields.size() < dim)
      {
        throw OomphLibError("Number of fields is not sufficient: " +
                              std::to_string(fields.size()) + "<" +
                              std::to_string(dim),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // We resize the invariants, the existing invariants should be unaffected
      I.resize(4);
      dIdG.resize(4, DenseMatrix<double>(dim));
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          I[3] += fields[i] * G(i, j) * fields[j];
          dIdG[3](i, j) = fields[i] * fields[j];
        }
      }
    }

    virtual void get_I_compressible(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& gup,
      const DenseMatrix<double>& Gup,
      const double& detg,
      const double& detG,
      const Vector<double>& fields,
      Vector<double>& I,
      Vector<DenseMatrix<double>>& dIdG) const override
    {
      StrainEnergyFunction::get_I_compressible(
        g, G, gup, Gup, detg, detG, fields, I, dIdG);
      get_I_anisotropic(g, G, fields, I, dIdG);
    }

    virtual void get_I_incompressible(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& gup,
      const DenseMatrix<double>& Gup,
      const double& detg,
      const double& detG,
      const Vector<double>& fields,
      Vector<double>& I,
      Vector<DenseMatrix<double>>& dIdG) const override
    {
      StrainEnergyFunction::get_I_incompressible(
        g, G, gup, Gup, detg, detG, fields, I, dIdG);
      get_I_anisotropic(g, G, fields, I, dIdG);
    }

    virtual void get_I_nearly_incompressible(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& gup,
      const DenseMatrix<double>& Gup,
      const double& detg,
      const double& detG,
      const Vector<double>& fields,
      Vector<double>& I,
      Vector<DenseMatrix<double>>& dIdG) const override
    {
      StrainEnergyFunction::get_I_nearly_incompressible(
        g, G, gup, Gup, detg, detG, fields, I, dIdG);
      get_I_anisotropic(g, G, fields, I, dIdG);
    }

  private:
    double* C3_pt;
  };

  class ThermalSofteningMooneyRivlin : public StrainEnergyFunction
  {
  public:
    ThermalSofteningMooneyRivlin(double* c1_pt, double* c2_pt)
      : StrainEnergyFunction(), C1_pt(c1_pt), C2_pt(c2_pt)
    {
    }

    virtual ~ThermalSofteningMooneyRivlin() {}

    double W(const Vector<double>& I) override
    {
      return (*C1_pt) * I[3] * (I[0] - 3.0) + (*C2_pt) * (I[1] - 3.0);
    }

    void derivatives(Vector<double>& I, Vector<double>& dWdI)
    {
      dWdI[0] = (*C1_pt) * I[3];
      dWdI[1] = (*C2_pt);
      dWdI[2] = 0.0;
      dWdI[3] = (*C1_pt) * (I[0] - 3.0);
    }

    void get_I_anisotropic(const DenseMatrix<double>& g,
                           const DenseMatrix<double>& G,
                           const Vector<double>& fields,
                           Vector<double>& I,
                           Vector<DenseMatrix<double>>& dIdG) const
    {
      const unsigned dim = g.nrow();
#ifdef RANGE_CHECKING
      if (fields.size() < 1)
      {
        throw OomphLibError("Number of fields is not sufficient",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // We resize the invariants, the existing invariants should be unaffected
      I.resize(4);
      dIdG.resize(4, DenseMatrix<double>(dim));
      I[3] = fields[0];
      dIdG[3] = DenseMatrix<double>(dim, dim, 0.0);
    }

    virtual void get_I_compressible(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& gup,
      const DenseMatrix<double>& Gup,
      const double& detg,
      const double& detG,
      const Vector<double>& fields,
      Vector<double>& I,
      Vector<DenseMatrix<double>>& dIdG) const override
    {
      StrainEnergyFunction::get_I_compressible(
        g, G, gup, Gup, detg, detG, fields, I, dIdG);
      get_I_anisotropic(g, G, fields, I, dIdG);
    }

    virtual void get_I_incompressible(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& gup,
      const DenseMatrix<double>& Gup,
      const double& detg,
      const double& detG,
      const Vector<double>& fields,
      Vector<double>& I,
      Vector<DenseMatrix<double>>& dIdG) const override
    {
      StrainEnergyFunction::get_I_incompressible(
        g, G, gup, Gup, detg, detG, fields, I, dIdG);
      get_I_anisotropic(g, G, fields, I, dIdG);
    }

    virtual void get_I_nearly_incompressible(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& gup,
      const DenseMatrix<double>& Gup,
      const double& detg,
      const double& detG,
      const Vector<double>& fields,
      Vector<double>& I,
      Vector<DenseMatrix<double>>& dIdG) const override
    {
      StrainEnergyFunction::get_I_nearly_incompressible(
        g, G, gup, Gup, detg, detG, fields, I, dIdG);
      get_I_anisotropic(g, G, fields, I, dIdG);
    }

  private:
    double* C1_pt;
    double* C2_pt;
    double* C3_pt;
  };

} // namespace oomph

#endif