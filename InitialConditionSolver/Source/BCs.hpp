#ifndef BCS_HPP
#define BCS_HPP

#include "BCFunc.cpp"

void SymmBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            Interval&       a_interval,
            int             a_order = 1)
{
  int isign = sign(a_side); // returns -1 for low, 1 for high
  RealVect facePos;

  // box, one cell thick, adjacant to valid. same as ghostboxlo.
  Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
  // intersection with box a_state (which has ghost cells, supposedly)
  toRegion &= a_state.box(); 

  Real* value = new Real[a_state.nComp()];

    // so this loop is only over outer ghost cells
  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
    {
      const IntVect& ivTo = bit();

      IntVect ivClose = ivTo -   isign*BASISV(a_dir);
      IntVect ivFar   = ivTo - 2*isign*BASISV(a_dir);

        // this is set to false in main
      if (!a_homogeneous) // if a_homogeneous if false
        {
            // finds face, since box is cell-centered
          getDomainFacePosition(facePos, ivClose, a_dx, a_dir, a_side);
          // fills value with right bdry conditions
          a_value(facePos.dataPtr(), &a_dir, &a_side, value);
        }

      for (int icomp = a_interval.begin(); icomp <= a_interval.end(); icomp++)
        {
            // comp value in bdry cell
          Real nearVal = a_state(ivClose, icomp);
          // comp value in next-to-bdry cell
          Real farVal  = a_state(ivFar,   icomp);

          Real inhomogVal = 0.;

          if (!a_homogeneous)
            {
              inhomogVal = value[icomp];
            }

          Real ghostVal=0;
           // inhomogVal = 3: we want this on the face
           // nearVal = 4
           // ghostVal = 2
          if (a_order == 1)
            {
              ghostVal = nearVal;
            }
          else if (a_order == 2)
            {
                MayDay::Error("quadratic interpolation not defined for symmetric boundary conditions");
            }
          else
            {
              MayDay::Error("bogus order argument");
            }

          a_state(ivTo, icomp) = ghostVal;
        }
    }

  delete[] value;
}

void SymmBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_dir,
            Side::LoHiSide  a_side)
{
    Interval stateInterval = a_state.interval();
    SymmBC(a_state, a_valid, a_dx, a_homogeneous, a_value, a_dir, a_side, 
        stateInterval);
}


#endif