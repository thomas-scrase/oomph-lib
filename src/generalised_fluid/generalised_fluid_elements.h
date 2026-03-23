#ifndef GENERALISED_FLUID_ELEMENTS_HEADER
#define GENERALISED_FLUID_ELEMENTS_HEADER

#include "../navier_stokes/navier_stokes_elements.h"
#include "generalised_fluid_constitutive_law.h"

namespace oomph
{
  template<unsigned DIM>
  class QGeneralisedFluidElement : public virtual QCrouzeixRaviartElement<DIM>
  {
  public:
    typedef void (*GeneralisedFluidThermalSourceFctPt)(const double& time,
                                                       const unsigned& ipt,
                                                       const Vector<double>& s,
                                                       const Vector<double>& x,
                                                       double& thermal_source);

    QGeneralisedFluidElement() : QCrouzeixRaviartElement<DIM>() {}

    virtual unsigned required_nvalue(const unsigned& n) const override
    {
      return QCrouzeixRaviartElement<DIM>::required_nvalue(n) + 2;
    }

    void set_incompressible()
    {
      Incompressible = true;
    }

    void set_compressible()
    {
      Incompressible = false;
    }

    GeneralisedFluidThermalSourceFctPt& thermal_source_fct_pt()
    {
      return Thermal_source_fct_pt;
    }

    GeneralisedFluidThermalSourceFctPt thermal_source_fct_pt() const
    {
      return Thermal_source_fct_pt;
    }

    GeneralisedFluidConstitutiveLaw*& constitutive_law_pt()
    {
      return Constitutive_law_pt;
    }

    GeneralisedFluidConstitutiveLaw* constitutive_law_pt() const
    {
      return Constitutive_law_pt;
    }

    void get_thermal_source(const double& time,
                            const unsigned& ipt,
                            const Vector<double>& s,
                            const Vector<double>& x,
                            double& thermal_source)
    {
      thermal_source = 0.0;
      if (Thermal_source_fct_pt)
      {
        (*Thermal_source_fct_pt)(time, ipt, s, x, thermal_source);
      }
    }

    virtual void get_equation_of_state(const double& time,
                                       const unsigned& ipt,
                                       const Vector<double>& s,
                                       const Vector<double>& x,
                                       const double& rho,
                                       const double& p,
                                       const double& e,
                                       double& residual) const
    {
      Constitutive_law_pt->get_equation_of_state(rho, p, e, residual);
    }

    // Abstracted to allow for general fields
    virtual void get_deviatoric_stress(const double& time,
                                       const unsigned& ipt,
                                       const Vector<double>& s,
                                       const Vector<double>& x,
                                       const DenseMatrix<double>& dudx,
                                       DenseMatrix<double>& sigma) const
    {
      Constitutive_law_pt->get_deviatoric_stress(dudx, sigma);
    }

    double drho_dt(const unsigned& l) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt = this->node_pt(l)->time_stepper_pt();

      // Initialise dudt
      double drhodt = 0.0;

      // Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Find the index at which the dof is stored
        const unsigned rho_nodal_index = this->rho_index_gfl();

        // Number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();
        // Loop over the timesteps
        for (unsigned t = 0; t < n_time; t++)
        {
          drhodt += time_stepper_pt->weight(1, t) *
                    this->nodal_value(t, l, rho_nodal_index);
        }
      }
      return drhodt;
    }

    double de_dt(const unsigned& l) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt = this->node_pt(l)->time_stepper_pt();

      // Initialise dudt
      double dedt = 0.0;

      // Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Find the index at which the dof is stored
        const unsigned e_nodal_index = this->e_index_gfl();

        // Number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();
        // Loop over the timesteps
        for (unsigned t = 0; t < n_time; t++)
        {
          dedt += time_stepper_pt->weight(1, t) *
                  this->nodal_value(t, l, e_nodal_index);
        }
      }
      return dedt;
    }

    double interpolated_rho(const Vector<double>& s)
    {
      // Find number of nodes
      unsigned n_node = this->nnode();
      // Local shape function
      Shape psi(n_node);
      // Find values of shape function
      this->shape(s, psi);
      // Index at which the nodal value is stored
      unsigned rho_nodal_index = rho_index_gfl();
      // Initialise value of rho
      double rho = 0.0;
      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        rho += this->nodal_value(l, rho_nodal_index) * psi[l];
      }
      return rho;
    }

    double interpolated_e(const Vector<double>& s)
    {
      // Find number of nodes
      unsigned n_node = this->nnode();
      // Local shape function
      Shape psi(n_node);
      // Find values of shape function
      this->shape(s, psi);
      // Index at which the nodal value is stored
      unsigned e_nodal_index = e_index_gfl();
      // Initialise value of rho
      double e = 0.0;
      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        e += this->nodal_value(l, e_nodal_index) * psi[l];
      }
      return e;
    }

    void output(std::ostream& outfile, const unsigned& nplot) override;

  protected:
    void fill_in_generic_residual_contribution_gfl(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag);

  public:
    void fill_in_contribution_to_residuals(Vector<double>& residuals) override
    {
      fill_in_generic_residual_contribution_gfl(
        residuals,
        GeneralisedElement::Dummy_matrix,
        GeneralisedElement::Dummy_matrix,
        0);
    }

    void fill_in_contribution_to_jacobian(
      Vector<double>& residuals, DenseMatrix<double>& jacobian) override
    {
      // Add the contribution to the residuals
      fill_in_contribution_to_residuals(residuals);
      // Allocate storage for the full residuals (residuals of entire element)
      unsigned n_dof = this->ndof();
      Vector<double> full_residuals(n_dof);
      // Get the residuals for the entire element
      this->get_residuals(full_residuals);
      this->fill_in_jacobian_from_internal_by_fd(residuals, jacobian, true);
      // Calculate the contributions from the external dofs
      //(finite-difference the lot by default)
      this->fill_in_jacobian_from_external_by_fd(residuals, jacobian, true);
      // Calculate the contributions from the nodal dofs
      this->fill_in_jacobian_from_nodal_by_fd(residuals, jacobian);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix) override
    {
      // fill_in_generic_residual_contribution_gfl(
      //   residuals, jacobian, mass_matrix, 2);
      throw OomphLibError("Not implemented mass matrix fill in",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

  protected:
    virtual inline unsigned rho_index_gfl() const
    {
      return DIM;
    }
    virtual inline unsigned e_index_gfl() const
    {
      return DIM + 1;
    }

  private:
    bool Incompressible = false;
    GeneralisedFluidThermalSourceFctPt Thermal_source_fct_pt;
    GeneralisedFluidConstitutiveLaw* Constitutive_law_pt =
      new GeneralisedFluidConstitutiveLaw;
  };

} // namespace oomph

#endif