#ifndef CHECK_SYMMETRIC_HPP
#define CHECK_SYMMETRIC_HPP

#include "MultigridUserVariables.hpp"
#include "PoissonParameters.H"
#include<iomanip>

// For now this checks symmetry on the low boundaries
int check_symmetric(LevelData<FArrayBox> &a_multigrid_vars, 
                     ProblemDomain &prob_domain, 
                     int comp, 
                     int iter,
                     int lev, 
                     const PoissonParameters &a_params
                    )
{
    if (a_params.periodic[0] == 1)
    {
        pout() << "Not checking symmetry because non-symmetric bdry conditions.\n";
        return 1;
    }
    // Loop over boxes
    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &a_multigrid_vars_box = a_multigrid_vars[dit()];
        Box b = a_multigrid_vars_box.box();
        b &= prob_domain.domainBox(); // rid the box of outer ghost cells
        IntVect low_end = b.smallEnd();
        // for (int dir = 0; dir < 3; dir++){pout() << low_end[dir] << endl;}
        // pout() << endl;
        // Stop if there are no outer ghost cells in this box
        // if (prob_domain.domainBox().contains(b))
        // {
        //     break;
        // }
        // Loop over cells in box
        BoxIterator bit(b);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            for (int dir = 0; dir < 3; dir++)
            {
                if (iv[dir] == 0) // if cell is boundary cell
                {
                    IntVect iv_ghost = iv;
                    iv_ghost[dir] -= 1; // set ghost cell
                    // evaluate psi at bdry and ghost cells
                    Real psi_bdry = a_multigrid_vars_box(iv, comp);
                    Real psi_ghost = a_multigrid_vars_box(iv_ghost, comp);
                    Real diff = sqrt((psi_bdry - psi_ghost) * (psi_bdry - psi_ghost));
                    //if (diff > 1.e-8) // compare values
                    if (psi_bdry != psi_ghost)
                    {
                        pout() <<   "Iter " << iter << ", level " << lev <<
                                    ", cell " << iv[0] << " " << iv[1] << " " << iv[2] << 
                                    ", " << "psi_bdry = " << psi_bdry << 
                                    " and psi_ghost = " << psi_ghost << endl;
                    }
                }
            }
        }
    }
    return 0;
}

int enforce_symmetric(LevelData<FArrayBox> &a_dpsi, 
                       ProblemDomain &prob_domain, 
                       int comp, 
                       int iter,
                       int lev, 
                       const PoissonParameters &a_params
                      )
{
    if (a_params.periodic[0] == 1)
    {
        pout() << "Not enforcing symmetry because non-symmetric bdry conditions.\n";
        return 1;
    }
    // Loop over boxes
    DataIterator dit = a_dpsi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &a_dpsi_box = a_dpsi[dit()];
        Box b = a_dpsi_box.box();
        b &= prob_domain.domainBox(); // rid the box of outer ghost cells
        BoxIterator bit(b);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            for (int dir = 0; dir < 3; dir++)
            {
                if (iv[dir] == 0)
                {
                    IntVect iv_ghost = iv;
                    iv_ghost[dir] -= 1; // set ghost cell
                    a_dpsi_box(iv_ghost, 0) = a_dpsi_box(iv, 0);
                }
                else if (iv[dir] == 20)
                {
                    IntVect iv_ghost = iv;
                    iv_ghost[dir] += 1; // set ghost cell
                    a_dpsi_box(iv_ghost, 0) = a_dpsi_box(iv, 0);
                }
            }
        }
    }
    return 0;
}

#endif