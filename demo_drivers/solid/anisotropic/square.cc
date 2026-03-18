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
// kruemelmonster
// Driver for a simple 2D anisotropic solid problem

// Generic oomph-lib routines
#include "generic.h"

// Include the mesh
#include "meshes/rectangular_quadmesh.h"

#include "constitutive.h"

#include "solid.h"

using namespace std;

using namespace oomph;

//==start_of_namespace================================================
/// Namespace for fish-shaped solution of 2D anisotropic pvd equation
//====================================================================
namespace GlobalParameters
{
  double c1 = 0.3;
  double c2 = 1.0;
  double c3 = 1.0;

  // StrainEnergyFunction* strain_energy_fct_pt =
  //   new FibreReinforcedMooneyRivlin(&c1, &c2, &c3);

  StrainEnergyFunction* strain_energy_fct_pt =
    new ThermalSofteningMooneyRivlin(&c1, &c2);
  // StrainEnergyFunction* strain_energy_fct_pt = new
  // GeneralisedMooneyRivlin(&c1,&c2);
  ConstitutiveLaw* constitutive_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(strain_energy_fct_pt);
  // ConstitutiveLaw* constitutive_law_pt = new GeneralisedHookean(&c1,&c2);

  double Gravity = 0.0;

  void gravity(const double& time, const Vector<double>& xi, Vector<double>& b)
  {
    b[0] = 0.0;
    b[1] = -Gravity;
  }

  void fields(const unsigned& ipt,
              const Vector<double>& s,
              const Vector<double>& xi,
              Vector<double>& fields)
  {
    // fields.resize(2);
    // fields[0] = cos(0.5 * MathematicalConstants::Pi * xi[0] / 4.0);
    // fields[1] = sin(0.5 * MathematicalConstants::Pi * xi[0] / 4.0);

    fields.resize(1);
    fields[0] = pow((4.0 - xi[0]) / 4.0, 3.0);
  }

} // namespace GlobalParameters


//==start_of_problem_class============================================
/// 2D anisotropic pvd problem
//====================================================================
template<class ELEMENT>
class AnisotropicSolidProblem : public Problem
{
public:
  /// Constructor: Pass number of elements and pointer to source function
  AnisotropicSolidProblem(const std::string& filename);

  /// Destructor (empty)
  ~AnisotropicSolidProblem()
  {
    delete mesh_pt();
  }

  /// Doc the solution, pass the number of the case considered,
  /// so that output files can be distinguished.
  void doc_solution(const unsigned& label);

  std::string Filename;

}; // end of problem class


//=====start_of_constructor===============================================
/// Constructor for 2D anisotropic pvd problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function.
//========================================================================
template<class ELEMENT>
AnisotropicSolidProblem<ELEMENT>::AnisotropicSolidProblem(
  const std::string& filename)
  : Filename(filename)
{
  // Build mesh and store pointer in Problem
  Problem::mesh_pt() =
    new ElasticRectangularQuadMesh<ELEMENT>(20, 10, 4.0, 2.0);

  // Set the boundary conditions for this problem: By default, all nodal
  // values are free -- we only need to pin the ones that have
  // Dirichlet conditions.

  for (unsigned l = 0; l < Problem::mesh_pt()->nboundary_node(2); l++)
  {
    dynamic_cast<SolidNode*>(mesh_pt()->boundary_node_pt(2, l))
      ->pin_position(0);
    dynamic_cast<SolidNode*>(mesh_pt()->boundary_node_pt(2, l))
      ->pin_position(1);
  }

  PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
    Problem::mesh_pt()->element_pt());

  // Loop over elements and set pointers to source function
  for (unsigned i = 0; i < Problem::mesh_pt()->nelement(); i++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* elem_pt =
      dynamic_cast<ELEMENT*>(Problem::mesh_pt()->element_pt(i));

    elem_pt->constitutive_law_pt() = GlobalParameters::constitutive_law_pt;

    elem_pt->set_incompressible();

    //  body_force_fct_pt
    elem_pt->body_force_fct_pt() = &GlobalParameters::gravity;

    //  fields_fct_pt
    elem_pt->fields_fct_pt() = &GlobalParameters::fields;
  }
  // dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_solid_pressure(0,0.0);

  // Setup equation numbering scheme
  assign_eqn_numbers();

} // end of constructor

//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void AnisotropicSolidProblem<ELEMENT>::doc_solution(const unsigned& label)
{
  using namespace StringConversion;

  // Number of plot points
  unsigned npts;
  npts = 5;

  // Output solution with specified number of plot points per element
  ofstream solution_file((Filename + to_string(label) + ".dat").c_str());
  mesh_pt()->output(solution_file, npts);
  solution_file.close();

} // end of doc


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 2D anisotropic pvd problem
//=====================================================================
int main()
{
  // Solve the thermal softening model
  {
    GlobalParameters::strain_energy_fct_pt = new ThermalSofteningMooneyRivlin(
      &GlobalParameters::c1, &GlobalParameters::c2);
    GlobalParameters::constitutive_law_pt =
      new IsotropicStrainEnergyFunctionConstitutiveLaw(
        GlobalParameters::strain_energy_fct_pt);

    AnisotropicSolidProblem<QPVDElementWithPressure<2>> problem(
      "ThermalSoftening");

    // Check whether the problem can be solved
    cout << "\n\n\nProblem self-test ";
    if (problem.self_test() == 0)
    {
      cout << "passed: Problem can be solved." << std::endl;
    }
    else
    {
      throw OomphLibError(
        "failed!", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    problem.doc_solution(0);
    for (unsigned i = 1; i < 10; i++)
    {
      // modify gravity
      problem.newton_solve();

      problem.doc_solution(i);

      GlobalParameters::Gravity += 1e-1;
    }
  }


  // Solve the fibre reinforced model
  {
    GlobalParameters::strain_energy_fct_pt = new FibreReinforcedMooneyRivlin(
      &GlobalParameters::c1, &GlobalParameters::c2, &GlobalParameters::c3);
    GlobalParameters::constitutive_law_pt =
      new IsotropicStrainEnergyFunctionConstitutiveLaw(
        GlobalParameters::strain_energy_fct_pt);

    AnisotropicSolidProblem<QPVDElementWithPressure<2>> problem(
      "FibreReinforced");

    // Check whether the problem can be solved
    cout << "\n\n\nProblem self-test ";
    if (problem.self_test() == 0)
    {
      cout << "passed: Problem can be solved." << std::endl;
    }
    else
    {
      throw OomphLibError(
        "failed!", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    problem.doc_solution(0);
    for (unsigned i = 1; i < 10; i++)
    {
      // modify gravity
      problem.newton_solve();

      problem.doc_solution(i);

      GlobalParameters::Gravity += 1e-1;
    }
  }

} // end of main
