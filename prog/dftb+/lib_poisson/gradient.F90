!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
module gradient


real(dp), allocatable :: gr(:,:)     !3,NNDIM
real(dp), allocatable :: hgrad(:,:)  !3,3*NNDIM
integer, allocatable :: conat(:)    !NNDIM+1
real(dp), allocatable :: convec(:,:) !3,NNDIM
logical :: constr


end module gradient
