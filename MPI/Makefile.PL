use ExtUtils::MakeMaker;
use FindBin;

$MPIDIR = "/usr/local/mpich/include";
foreach (qw(/usr/mpich/include /usr/local/mpich/include
	    /usr/include/mpi /usr/local/include/mpi
	    /usr/lib/mpich/include /usr/local/lib/mpich/include))
{
    if(-e "$_/mpi.h") { $MPIDIR = $_; last }
}

WriteMakefile(
    'NAME'	=> 'Parallel::MPIcar',
    'VERSION_FROM' => 'MPIcar.pm',
    'LD'        => 'mpicc',
    'CC'        => 'mpicc',      
    'LIBS'	=> [''],
    'DEFINE'	=> '-DFLOAT_HACK',
    'INC'	=> "-I$MPIDIR",
    'EXE_FILES' => ['mpi_maker'],
    'DESTDIR' => "$FindBin::Bin/../",
    'INSTALLDIRS' => "perl",
    'INSTALLMAN1DIR' => "lib/Parallel/MPIcar/docs",
    'INSTALLMAN3DIR' => "lib/Parallel/MPIcar/docs",
    'INSTALLARCHLIB' => "lib/",
    'INSTALLPRIVLIB' => "lib/",
    'INSTALLSCRIPT' => "bin/",
    'INSTALLBIN' => "bin/",
    PL_FILES    => { 't/01_mpi_init.PL'          => 't/01_mpi_init.pl',
                     't/02_mpi_finalize.PL'      => 't/02_mpi_finalize.pl',
                     't/03_mpi_comm_rank.PL'     => 't/03_mpi_comm_rank.pl',
                     't/04_mpi_comm_size.PL'     => 't/04_mpi_comm_size.pl',
                     't/05_mpi_initialized.PL'   => 't/05_mpi_initialized.pl',
                     't/06_mpi_sendrecv.PL'      => 't/06_mpi_sendrecv.pl',
                     't/06_sendrecv_int.PL'      => 't/06_sendrecv_int.pl',
                     't/06_sendrecv_intarray.PL' => 't/06_sendrecv_intarray.pl',
                     't/06_sendrecv_dbl.PL'      => 't/06_sendrecv_dbl.pl',
                     't/06_sendrecv_dblarray.PL' => 't/06_sendrecv_dblarray.pl',
                     't/06_sendrecv_flt.PL'      => 't/06_sendrecv_flt.pl',
                     't/06_sendrecv_fltarray.PL' => 't/06_sendrecv_fltarray.pl',
                     't/07_mpi_bcast.PL'         => 't/07_mpi_bcast.pl',
                     't/07_mpi_bcastarray.PL'    => 't/07_mpi_bcastarray.pl',
                     't/08_mpi_reduce.PL'        => 't/08_mpi_reduce.pl',
                     't/09_mpi_barrier.PL'       => 't/09_mpi_barrier.pl',
                     't/10_mpi_wtime.PL'         => 't/10_mpi_wtime.pl',
                     't/11_mpi_scattergather.PL' => 't/11_mpi_scattergather.pl',
                     't/11_mpi_scattergatherarray.PL' => 't/11_mpi_scattergatherarray.pl',
                     't/12_mpi_allreduce.PL'     => 't/12_mpi_allreduce.pl',
                     't/13_mpi_sendrecv.PL'     => 't/13_mpi_sendrecv.pl',
                   },
    dist        => { COMPRESS => 'gzip -9 -f', SUFFIX => 'gz' },
    clean       => { FILES => 't/*.pl PI[0-9][0-9][0-9]*'},
);
