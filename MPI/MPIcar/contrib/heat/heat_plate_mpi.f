        include 'mpif.h'
        real, dimension (:,:), allocatable :: t,s,r
        real (kind=8):: timef,t2,t1
        real :: h,xl,yl,tol,h2,aa,sum,global_sum
        integer :: nx,ny,i,j,ij,nproc,iter,niter,n_left,npt,nxy
        integer :: ierror,nproc,my_rank,root
        integer, dimension(MPI_STATUS_SIZE) :: istatus
        integer, dimension(:), allocatable :: local_n
        real, dimension(:), allocatable :: local_a
c234567
       call mpi_init(ierror)
       call mpi_comm_size(MPI_COMM_WORLD, nproc, ierror)
       call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierror)
c
      xl = 10.
      yl = 10.
      tol = 1.0e-3
c
c initial set up done by my_rank = 0 (i.e. master) process only
c
      if(my_rank .eq. 0) then
      print *,'there are ',nproc,' processes'
      print *,'enter number of iterations and h:'
      read(5,*) niter,h
c
c  determine the workload for each processor
c
      nx = xl/h + 1
      ny = nx
c
      print *,'nx=ny= ',nx
c
      allocate(local_n(nproc))
      allocate(local_a(nproc))
c
      do i=1,nproc
        local_n(i) = ny/nproc 
      enddo
c
c if workload cannot be divided evenly
c to each process, assign to the
c first n_left processes with an
c additional piece of work
c
      if(mod(ny,nproc) .ne. 0) then
       n_left = ny - nproc*(ny/nproc)
        do i=1,n_left
          local_n(i) = local_n(i) + 1
        enddo
        local_a(1) = 0.0
        do i=2,nproc
          local_a(i) = local_a(i-1) + h*local_n(i-1)
        enddo
      endif

c
      print *,'process #,  local worklaod,  start location '
      do i=1,nproc
       print *,i-1, local_n(i), local_a(i)
      enddo

      endif

c
c  now all the processes begin to work together
c
      call MPI_Bcast(h,1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(niter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Scatter(local_n,1,MPI_INTEGER,
     &                 npt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Scatter(local_a,1,MPI_REAL,
     &                 aa,1,MPI_REAL,0,MPI_COMM_WORLD,ierror)
c
c  add one additional point to the end strips and 
c  add additional 2 points to the interior strips 
c  to facilitate overlap      
c
      xl = 10.0
      nx = xl/h + 1

      if( my_rank.eq.0 .or. my_rank.eq.nproc-1) then
          ny = npt + 1
      else
          ny = npt + 2
      endif
c
c
      allocate (t(nx,ny))
      allocate (s(nx,ny))
      allocate (r(nx,ny))
c
c      print *,'I am process ',my_rank,' nx,ny= ',nx,ny
c      print *,' aa= ',aa
c      print *,' h = ',h
c      print *,' niter= ',niter
c
      h2 = h*h
      nxy = nx*ny
c
c
c   source term
c
      do j=1,ny
       do i=1,nx
       t(i,j) = 20. 
       s(i,j) = 0.0
       r(i,j) = 0.0
       enddo
      enddo
c
c  fix the boundary conditions
c
c  domain is divided into horizontal stripes
c

c
c    along the west and east walls
c

      do j=1,ny
      t(1,j)  = 20.0
      t(nx,j) = 20.0
      enddo

c 
c   along the south and north walls
c

c234567
c
c  south boundary is owned by the first processor
c
      if( my_rank .eq. 0) then
      do i=1,nx
       t(i,1) = 20. 
      enddo
      endif
c
c  north boundary is owned by the last processor
c 
      if( my_rank .eq. nproc-1) then
      do i=1,nx
c
      xx = (i-1)*h
      if ( xx.ge.3.0 .and. xx.le.7.0 ) then
       t(i,ny) = 100.0
      else
       t(i,ny) = 20.0
      endif

      enddo
      endif


      call MPI_Barrier(MPI_COMM_WORLD,ierror)

      if( my_rank .eq. 0) then
       print *,'finish updating boundary conditions'
c       t1 = timef()
       t1 = MPI_Wtime()
      endif


      do iter=1,niter

c
c send and receive interfacial values from nearest neighbours
c
      if( my_rank.ne.0   .and.  my_rank.ne.nproc-1 ) then
c
c along the south side
c
        call MPI_Sendrecv(t(1,2),nx,MPI_REAL,my_rank-1,10,
     &                    t(1,1),nx,MPI_REAL,my_rank-1,20,
     &                    MPI_COMM_WORLD,istatus,ierror)
c
c along the north side
c
        call MPI_Sendrecv(t(1,ny-1),nx,MPI_REAL,my_rank+1,20,
     &                    t(1,ny),  nx,MPI_REAL,my_rank+1,10,
     &                    MPI_COMM_WORLD,istatus,ierror)

      endif

      if( my_rank .eq. 0 ) then
c
c along the north side
c
        root = 1
        call MPI_Sendrecv(t(1,ny-1),nx,MPI_REAL,root,20,
     &                    t(1,ny),  nx,MPI_REAL,root,10,
     &                    MPI_COMM_WORLD,istatus,ierror)

      endif


      if( my_rank .eq. nproc-1 ) then
c
c along the south side
c
        root = nproc-2
        call MPI_Sendrecv(t(1,2),nx,MPI_REAL,root,10,
     &                    t(1,1),nx,MPI_REAL,root,20,
     &                    MPI_COMM_WORLD,istatus,ierror)

      endif

     
      sum = 0.0

      do j=2,ny-1
       do i=2,nx-1
        r(i,j) = s(i,j)*h2 - (
     1                         t(i+1,j)+t(i,j+1)
     2                        -4.*t(i,j)
     3                        +t(i-1,j)+t(i,j-1)  )
        sum = sum + r(i,j)*r(i,j)
       enddo
      enddo 
c
c  update temperatures
c
      do ij=1,nxy
        t(ij,1) = t(ij,1) - 0.25*r(ij,1)
      enddo

c
c obtain global sum
c
      call MPI_Allreduce(sum,global_sum,1,MPI_REAL,
     &                   MPI_SUM,MPI_COMM_WORLD,ierror)
c
c
      global_sum = sqrt(global_sum)
c
c
c      if( my_rank .eq. 0) then
c      if (mod(iter,500) .eq. 0) then
c      print *,iter,global_sum
c      endif
c
c      write(2,*) iter, global_sum
c
c      endif
c
      if(global_sum .le. tol ) go to 1000
c
      enddo
c
1000     continue
c
c
        if( my_rank .eq. 0) then
c        t2 = timef()
        t2 = MPI_Wtime()
        print *,'no. of iterations took = ',iter
        print *,'final L2 norm of residual: ',global_sum
        print *,'wall clock time in seconds= ',(t2-t1)
        endif
c
        call MPI_FINALIZE(ierror)
c
      stop
      end

