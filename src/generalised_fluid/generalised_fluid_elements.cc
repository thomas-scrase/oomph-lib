#include "generalised_fluid_elements.h"

namespace oomph
{
  template<unsigned DIM>
  void QGeneralisedFluidElement<DIM>::fill_in_generic_residual_contribution_gfl(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix,
    unsigned flag)
  {
    // Return immediately if there are no dofs
    if (this->ndof() == 0) return;

    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Get continuous time from timestepper of first node
    double time = this->node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find out how many pressure dofs there are
    unsigned n_pres = this->npres_nst();

    // Find the indices at which the local velocities are stored
    unsigned u_nodal_index[DIM];
    unsigned rho_nodal_index = rho_index_gfl();
    unsigned e_nodal_index = e_index_gfl();
    for (unsigned i = 0; i < DIM; i++)
    {
      u_nodal_index[i] = this->u_index_nst(i);
    }

    // Set up memory for the shape and test functions
    Shape psif(n_node), testf(n_node);
    DShape dpsifdx(n_node, DIM), dtestfdx(n_node, DIM);

    // Set up memory for pressure shape and test functions
    Shape psip(n_pres), testp(n_pres);

    // Number of integration points
    unsigned n_intpt = this->integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Integers to store the local equations and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
        s[i] = this->integral_pt()->knot(ipt, i);
      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_nst(
        ipt, psif, dpsifdx, testf, dtestfdx);

      // Call the pressure shape and test functions
      this->pshape_nst(s, psip, testp);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate local values of the pressure and velocity components
      // Allocate
      double interpolated_p = 0.0;
      Vector<double> interpolated_u(DIM, 0.0);
      double interpolated_rho = 0.0;
      double interpolated_e = 0.0;
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);
      Vector<double> dudt(DIM, 0.0);
      double drhodt = 0.0;
      double dedt = 0.0;
      DenseMatrix<double> interpolated_dudx(DIM, DIM, 0.0);
      Vector<double> interpolated_drhodx(DIM, 0.0);
      Vector<double> interpolated_dedx(DIM, 0.0);

      Vector<double> G = this->g();

      // Calculate pressure
      for (unsigned l = 0; l < n_pres; l++)
        interpolated_p += this->p_nst(l) * psip[l];

      // Calculate velocities and derivatives:

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_rho += this->raw_nodal_value(l, rho_nodal_index) * psif[l];
        interpolated_e += this->raw_nodal_value(l, e_nodal_index) * psif[l];
        drhodt += drho_dt(l) * psif[l];
        dedt += de_dt(l) * psif[l];
        // Loop over directions
        for (unsigned i = 0; i < DIM; i++)
        {
          // Get the nodal value
          double u_value = this->raw_nodal_value(l, u_nodal_index[i]);
          interpolated_u[i] += u_value * psif[l];
          interpolated_x[i] += this->raw_nodal_position(l, i) * psif[l];
          dudt[i] += this->du_dt_nst(l, i) * psif[l];
          interpolated_drhodx[i] +=
            this->raw_nodal_value(l, rho_nodal_index) * dpsifdx(l, i);
          interpolated_dedx[i] +=
            this->raw_nodal_value(l, e_nodal_index) * dpsifdx(l, i);

          // Loop over derivative directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(i, j) += u_value * dpsifdx(l, j);
          }
        }
      }

      if (!this->ALE_is_disabled)
      {
        // Loop over nodes
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < DIM; i++)
          {
            mesh_velocity[i] += this->raw_dnodal_position_dt(l, i) * psif[l];
          }
        }
      }

      // Get the user-defined body force terms
      Vector<double> body_force(DIM);
      this->get_body_force_nst(time, ipt, s, interpolated_x, body_force);

      // Get the user-defined source function
      double source = this->get_source_nst(time, ipt, interpolated_x);

      double thermal_source;
      get_thermal_source(time, ipt, s, interpolated_x, thermal_source);

      Vector<double> thermal_flux(DIM, 0.0);
      for (unsigned i = 0; i < DIM; i++)
      {
        // Replace this with constitutive law for diffusion?
        thermal_flux[i] *= interpolated_dedx[i];
      }

      DenseMatrix<double> deviatoric_stress(DIM, DIM, 0.0);
      get_deviatoric_stress(
        time, ipt, s, interpolated_x, interpolated_dudx, deviatoric_stress);

      double equation_of_state_residual;
      get_equation_of_state(time,
                            ipt,
                            s,
                            interpolated_x,
                            interpolated_rho,
                            interpolated_p,
                            interpolated_e,
                            equation_of_state_residual);

      // Continuity equation (p residual) -> switch incompressible
      // Loop over the velocity components
      /*IF it's not a boundary condition*/
      for (unsigned l = 0; l < n_pres; l++)
      {
        local_eqn = this->p_local_eqn(l);
        if (local_eqn >= 0)
        {
          residuals[local_eqn] -= source * testp[l] * W;
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] -=
              interpolated_rho * interpolated_dudx(k, k) * testp[l] * W;
          }

          if (!Incompressible)
          {
            residuals[local_eqn] -= drhodt * testp[l] * W;
            for (unsigned k = 0; k < DIM; k++)
            {
              double tmp = interpolated_u[k];
              if (!this->ALE_is_disabled) tmp -= mesh_velocity[k];
              residuals[local_eqn] -=
                tmp * interpolated_drhodx[k] * testp[l] * W;
            }
          }
        }
      }

      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the velocity components
        for (unsigned i = 0; i < DIM; i++)
        {
          /*IF it's not a boundary condition*/
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[i]);
          if (local_eqn >= 0)
          {
            residuals[local_eqn] += body_force[i] * testf[l] * W;

            residuals[local_eqn] += testf[l] * G[i] * W;

            residuals[local_eqn] += interpolated_p * dtestfdx(l, i) * W;


            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn] -=
                deviatoric_stress(i, k) * dtestfdx(l, k) * W;
            }

            residuals[local_eqn] -= interpolated_rho * dudt[i] * testf[l] * W;

            for (unsigned k = 0; k < DIM; k++)
            {
              double tmp = interpolated_u[k];
              if (!this->ALE_is_disabled) tmp -= mesh_velocity[k];
              residuals[local_eqn] -=
                interpolated_rho * tmp * interpolated_dudx(i, k) * testf[l] * W;
            }
          }
        }

        // Fill in rho residual -> switch incompressible
        /*IF it's not a boundary condition*/
        local_eqn = this->nodal_local_eqn(l, rho_nodal_index);
        if (local_eqn >= 0)
        {
          if (Incompressible)
          {
            residuals[local_eqn] -= drhodt * testf[l] * W;

            for (unsigned k = 0; k < DIM; k++)
            {
              double tmp = interpolated_u[k];
              if (!this->ALE_is_disabled) tmp -= mesh_velocity[k];
              residuals[local_eqn] -=
                tmp * interpolated_drhodx[k] * testf[l] * W;
            }
          }
          else
          {
            residuals[local_eqn] -= equation_of_state_residual * testf[l] * W;
          }
        }

        // Fill in e residual
        /*IF it's not a boundary condition*/
        local_eqn = this->nodal_local_eqn(l, e_nodal_index);
        if (local_eqn >= 0)
        {
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] += thermal_flux[k] * dtestfdx(l, k) * W;

            for (unsigned kk = 0; kk < DIM; kk++)
            {
              residuals[local_eqn] -= deviatoric_stress(k, kk) *
                                      interpolated_dudx(k, kk) * testf[l] * W;
            }
          }

          residuals[local_eqn] -= interpolated_rho * dedt * testf[l] * W;

          for (unsigned k = 0; k < DIM; k++)
          {
            double tmp = interpolated_u[k];
            if (!this->ALE_is_disabled) tmp -= mesh_velocity[k];
            residuals[local_eqn] -=
              interpolated_rho * tmp * interpolated_dedx[k] * testf[l] * W;

            residuals[local_eqn] -=
              interpolated_rho * thermal_source * testf[l] * W;
          }
        }
      }
    }
  }

  template<unsigned DIM>
  void QGeneralisedFluidElement<DIM>::output(std::ostream& outfile,
                                             const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Coordinates
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << this->interpolated_x(s, i) << " ";
      }

      // Velocities
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << this->interpolated_u_nst(s, i) << " ";
      }

      // Pressure
      outfile << this->interpolated_p_nst(s) << " ";

      outfile << interpolated_rho(s) << " ";
      outfile << interpolated_e(s) << " ";

      outfile << std::endl;
    }
    outfile << std::endl;

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  // Template instantiations
  template class QGeneralisedFluidElement<2>;
  template class QGeneralisedFluidElement<3>;

} // namespace oomph