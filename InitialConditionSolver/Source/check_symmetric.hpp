#ifndef CHECK_SYMMETRIC_HPP
#define CHECK_SYMMETRIC_HPP

#include "MultigridUserVariables.hpp"

// For now this checks symmetry on the low boundaries
void check_symmetric(LevelData<FArrayBox> &a_multigrid_vars, 
                     int iter, 
                     ProblemDomain &prob_domain
                    )
{
    // Loop over boxes
    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &a_multigrid_vars_box = a_multigrid_vars[dit()];
        Box b = a_multigrid_vars_box.box();
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
                    Real psi_bdry = a_multigrid_vars_box(iv, c_psi_reg);
                    Real psi_ghost = a_multigrid_vars_box(iv_ghost, c_psi_reg);
                    Real diff = abs(psi_bdry - psi_ghost);
                    if (diff > 1.e-5) // compare values
                    {
                        pout() <<   "At iter " << iter <<
                                    " , at cell " << iv[0] << " " << iv[1] << " " << iv[2] << 
                                    ", " << "psi_bdry = " << psi_bdry << 
                                    " and psi_ghost = " << psi_ghost << endl;
                    }
                }
            }
        }
    }
}

#endif