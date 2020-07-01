#ifndef CHECK_SYMMETRIC_HPP
#define CHECK_SYMMETRIC_HPP

#include "MultigridUserVariables.hpp"
#include "PoissonParameters.H"
#include<iomanip>

namespace 
{
    int high_end = 20;
}

int enforce_var_symmetric(LevelData<FArrayBox> &a_data, 
                          const int var,
                          const ProblemDomain &prob_domain,
                          const PoissonParameters &a_params, 
                          const int ilev)
{
    DataIterator dit = a_data.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &a_data_fabox = a_data[dit()];
        Box a_data_box = a_data_fabox.box();

        // We must check that there's any outer ghost cells in this box at all
        if (prob_domain.contains(a_data_box))
        {
            return 1;
        }

        bool we_have_ghosts = false; // Use this for consistency check
        for (int dir = 0; dir < SpaceDim; dir++)
        {
            // Must check that we have ghosts on the low side
            if (a_data_box.smallEnd(dir) < prob_domain.domainBox().smallEnd(dir))
            {
                we_have_ghosts = true;

                // Intersection with problemdomain to remove outer ghost cells
                a_data_box &= prob_domain.domainBox();

                // Construct box with exclusively outer ghost cells, ...
                Box lo_outer_ghosts = adjCellLo(a_data_box, dir, 3);
                // ... which is the box that we'll loop over
                BoxIterator bit(lo_outer_ghosts);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    IntVect iv_to = bit();
                    IntVect iv_from = iv_to;
                    iv_from[dir] += -2 * iv_from[dir] - 1;
                    a_data_fabox(iv_to, var) = a_data_fabox(iv_from, var);
                    // for (int dir = 0; dir < SpaceDim; dir++)
                    // {
                    //     pout() << iv[0] << " " << iv[1] << " " << iv[2] << endl;
                    // }
                    // if (iv[0] > -1 && iv[1] > -1 && iv[2] > -1)
                    // {
                    //     pout() << "Something is wrong!!\n";
                    // }
                } // end loop over IntVects
            }

            // Must check that we have ghosts on the high side
            if (a_data_box.bigEnd(dir) > prob_domain.domainBox().bigEnd(dir))
            {
                we_have_ghosts = true;
                a_data_box &= prob_domain.domainBox();

                Box hi_outer_ghosts = adjCellHi(a_data_box, dir, 3);
                BoxIterator bit(hi_outer_ghosts);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    IntVect iv_to = bit();
                    IntVect iv_from = iv_to;
                    iv_from[dir] -= 2 * (iv_from[dir] - prob_domain.domainBox().bigEnd(dir)) + 1;
                    a_data_fabox(iv_to, var) = a_data_fabox(iv_from, var);
                } // end loop over IntVects
            }
            if (we_have_ghosts == false)
            {
                pout() << "There's an issue: probdomain doesn't contain the box, but there are no ghosts\n";
            }
        }
    } // end loop over boxes
    return 0;
}

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
                else if (iv[dir] == high_end)
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