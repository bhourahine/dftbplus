!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains MPI related einvironment settings
module mpienv
  use accuracy, only : lc
  use mpifx
  use message
  implicit none
  private

  public :: TMpiEnv, TMpiEnv_init


  !> Contains MPI related environment settings
  type :: TMpiEnv

    !> Global MPI communicator
    type(mpifx_comm) :: globalComm

    !> Communicator to access processes within current group
    type(mpifx_comm) :: groupComm

    !> Communicator to access equivalent processes in other groups
    type(mpifx_comm) :: interGroupComm

    !> Size of the process groups
    integer :: groupSize

    !> Number of processor groups
    integer :: nGroup

    !> Group index of the current process (starts with 0)
    integer :: myGroup

    !> Global rank of the processes in the given group
    integer, allocatable :: groupMembers(:)

    !> Whether current process is the global master
    logical :: tGlobalMaster

    !> Whether current process is the group master
    logical :: tGroupMaster

    !> Communicator to access processes within a replica
    type(mpifx_comm) :: intraReplicaComm

    !> Communicator to between replicas
    type(mpifx_comm) :: interReplicaComm

    !> Size of the group of replica
    integer :: replicaCommSize

    !> Replica group index the current process belongs to (starts with 0)
    integer :: myReplica

    !> Global rank of the processes in the given replica group
    integer, allocatable :: replicaMembers(:)

    !> Whether current process is the replica master
    logical :: tReplicaMaster

  end type TMpiEnv


contains

  !> Initializes MPI environment.
  subroutine TMpiEnv_init(this, nGroup, nReplicas)

    !> Initialised instance on exit
    type(TMpiEnv), intent(out) :: this

    !> Number of process groups to create
    integer, intent(in) :: nGroup

    !> Number of structure replicas
    integer, intent(in) :: nReplicas

    character(lc) :: tmpStr
    integer :: myRank, myGroup

    call this%globalComm%init()

    this%nGroup = nGroup * nReplicas
    this%groupSize = this%globalComm%size / this%nGroup

    if (this%nGroup * this%groupSize /= this%globalComm%size) then
      write(tmpStr, "(A,I0,A,I0,A)") "Number of groups (", this%nGroup,&
          & ") not compatible with number of processes (", this%globalComm%size, ")"
      call error(tmpStr)
    end if

    ! group this process belongs to
    this%myGroup = this%globalComm%rank / this%groupSize

    ! rank within the group
    myRank = mod(this%globalComm%rank, this%groupSize)

    ! communicator within this group
    call this%globalComm%split(this%myGroup, myRank, this%groupComm)

    ! array of global process ids within this group
    allocate(this%groupMembers(this%groupSize))
    call mpifx_allgather(this%groupComm, this%globalComm%rank, this%groupMembers)

    myGroup = myRank
    myRank = this%myGroup
    call this%globalComm%split(myGroup, myRank, this%interGroupComm)

    ! processors in a replica group
    this%replicaCommSize = this%globalComm%size / nReplicas
    ! replica group this process belongs to
    this%myReplica = this%globalComm%rank / max(nReplicas-1,1)
    ! rank within the replica group
    myRank = mod(this%globalComm%rank, this%replicaCommSize)
    ! communicator within this replica
    call this%globalComm%split(this%myReplica, myRank, this%intraReplicaComm)

    ! array of global process ids within this replica
    allocate(this%replicaMembers(this%replicaCommSize))
    call mpifx_allgather(this%intraReplicaComm, this%globalComm%rank, this%replicaMembers)

    ! communicator to equivalent processors in different replicas
    myGroup = myRank
    myRank = this%myReplica
    call this%globalComm%split(myGroup, myRank, this%interReplicaComm)

    this%tGlobalMaster = this%globalComm%master
    this%tGroupMaster = this%groupComm%master
    this%tReplicaMaster = this%intraReplicaComm%master

    if (this%tGlobalMaster .and. .not. this%tGroupMaster) then
      call error("Internal error: Global master process is not a group master process")
    end if

    if (this%tGlobalMaster .and. .not. this%tReplicaMaster) then
      call error("Internal error: Global master process is not a replica master process")
    end if

    if (this%tReplicaMaster .and. .not. this%tGroupMaster) then
      call error("Internal error: Replica master process is not a group master process")
    end if

  end subroutine TMpiEnv_init


end module mpienv
