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
#    'LIBS'	=> [''],
    'DEFINE'	=> '-DFLOAT_HACK',
    'INC'	=> "-I$MPIDIR",
#    'EXE_FILES' => [''],
    'DESTDIR' => "$FindBin::Bin/../../perl/",
    'INSTALLDIRS' => "perl",
    'INSTALLMAN1DIR' => "none",
    'INSTALLMAN3DIR' => "none",
    'INSTALLARCHLIB' => "lib/",
    'INSTALLPRIVLIB' => "lib/",
#    'INSTALLSCRIPT' => "bin/",
#    'INSTALLBIN' => "bin/",
    PL_FILES    => { 't/01_mpi_init'          => 't/01_mpi_init.pl',
                     't/02_mpi_finalize'      => 't/02_mpi_finalize.pl',
                     't/03_mpi_comm_rank'     => 't/03_mpi_comm_rank.pl',
                     't/04_mpi_comm_size'     => 't/04_mpi_comm_size.pl',
                     't/05_mpi_initialized'   => 't/05_mpi_initialized.pl',
                     't/06_mpi_sendrecv'      => 't/06_mpi_sendrecv.pl',
                     't/06_sendrecv_int'      => 't/06_sendrecv_int.pl',
                     't/06_sendrecv_intarray' => 't/06_sendrecv_intarray.pl',
                     't/06_sendrecv_dbl'      => 't/06_sendrecv_dbl.pl',
                     't/06_sendrecv_dblarray' => 't/06_sendrecv_dblarray.pl',
                     't/06_sendrecv_flt'      => 't/06_sendrecv_flt.pl',
                     't/06_sendrecv_fltarray' => 't/06_sendrecv_fltarray.pl',
                     't/07_mpi_bcast'         => 't/07_mpi_bcast.pl',
                     't/07_mpi_bcastarray'    => 't/07_mpi_bcastarray.pl',
                     't/08_mpi_reduce'        => 't/08_mpi_reduce.pl',
                     't/09_mpi_barrier'       => 't/09_mpi_barrier.pl',
                     't/10_mpi_wtime'         => 't/10_mpi_wtime.pl',
                     't/11_mpi_scattergather' => 't/11_mpi_scattergather.pl',
                     't/11_mpi_scattergatherarray' => 't/11_mpi_scattergatherarray.pl',
                     't/12_mpi_allreduce'     => 't/12_mpi_allreduce.pl',
                     't/13_mpi_sendrecv'     => 't/13_mpi_sendrecv.pl',
                   },
    dist        => { COMPRESS => 'gzip -9 -f', SUFFIX => 'gz' },
    clean       => { FILES => 't/*.pl PI[0-9][0-9][0-9]*'},
);
