subroutine error ( x, err )

!*****************************************************************************80
!
!! ERROR evaluates the error function.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    07 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) ERR, the function value.
!
  implicit none

  real ( kind = 8 ) c0
  real ( kind = 8 ) eps
  real ( kind = 8 ) er
  real ( kind = 8 ) err
  integer ( kind = 4 ) k
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  eps = 1.0D-15
  pi = 3.141592653589793D+00
  x2 = x * x

  if ( abs ( x ) < 3.5D+00 ) then

    er = 1.0D+00
    r = 1.0D+00

    do k = 1, 50
      r = r * x2 / ( k + 0.5D+00 )
      er = er + r
      if ( abs ( r ) <= abs ( er ) * eps ) then
        exit
      end if
    end do

    c0 = 2.0D+00 / sqrt ( pi ) * x * exp ( - x2 )
    err = c0 * er

  else

    er = 1.0D+00
    r = 1.0D+00
    do k = 1, 12
      r = - r * ( k - 0.5D+00 ) / x2
      er = er + r
    end do

    c0 = exp ( - x2 ) / ( abs ( x ) * sqrt ( pi ) )

    err = 1.0D+00 - c0 * er
    if ( x < 0.0D+00 ) then
      err = -err
    end if

  end if

  return
end
