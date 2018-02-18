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

    !> Rank within my group
    integer :: myGroupRank

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

    !> Rank within my replica
    integer :: myReplicaRank

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

    call this%globalComm%init()

    ! processors in a replica group
    this%replicaCommSize = this%globalComm%size / nReplicas

    if (nReplicas * this%replicaCommSize /= this%globalComm%size) then
      write(tmpStr, "(A,I0,A,I0,A)") "Number of replicas (", nReplicas,&
          & ") not compatible with number of processes (", this%globalComm%size, ")"
      call error(tmpStr)
    end if

    ! replica group this process belongs to
    this%myReplica = this%globalComm%rank / this%replicaCommSize

    ! rank within the replica group
    this%myReplicaRank = mod(this%globalComm%rank, this%replicaCommSize)

    ! communicator within this replica
    call this%globalComm%split(this%myReplica, this%myReplicaRank, this%intraReplicaComm)

    ! array of global process ids within this replica
    allocate(this%replicaMembers(this%replicaCommSize))
    call mpifx_allgather(this%intraReplicaComm, this%globalComm%rank, this%replicaMembers)

    ! communicator to equivalent processors in different replicas
    call this%globalComm%split(this%myReplicaRank, this%myReplica, this%interReplicaComm)

    ! groups within a replica
    this%nGroup = nGroup

    this%groupSize = this%intraReplicaComm%size / this%nGroup

    if (this%nGroup * this%groupSize /= this%intraReplicaComm%size) then
      write(tmpStr, "(A,I0,A,I0,A)") "Number of groups (", this%nGroup,&
          & ") not compatible with number of processes per replica (",&
          & this%globalComm%size/nReplicas, ")"
      call error(tmpStr)
    end if

    ! group this process belongs to
    this%myGroup = this%intraReplicaComm%rank / this%groupSize

    ! rank within the group
    this%myGroupRank = mod(this%intraReplicaComm%rank, this%groupSize)

    ! communicator within this group
    call this%intraReplicaComm%split(this%myGroup, this%myGroupRank, this%groupComm)

    ! array of global process ids within this group
    allocate(this%groupMembers(this%groupSize))
    call mpifx_allgather(this%groupComm, this%globalComm%rank, this%groupMembers)

    ! communicator to equivalent processors in different groups of the same replica
    call this%intraReplicaComm%split(this%myGroupRank, this%myGroup, this%interGroupComm)

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
