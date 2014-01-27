package Parallel::Application;

use strict;
use Carp;
use vars qw(@ISA $VERSION %EXPORT_TAGS @EXPORT_OK);

require Exporter;

@ISA = qw(Exporter);
$VERSION = '0.01';

my @functions = qw(App_Init
                   App_Finalize
                   App_Comm_rank
                   App_Comm_size
                   App_Send
                   App_Recv);

my @constants = qw(MPI_ANY_SOURCE
                   MPI_ANY_TAG);

%EXPORT_TAGS = ( all => [@functions, @constants] );
@EXPORT_OK = ( @functions, @constants);

sub App_Init { croak "FATAL: App_Init must be overridden by the Parallel::Application child\n"; }
sub App_Finalize { croak "FATAL: App_Finalize must be overridden by the Parallel::Application child\n"; }
sub App_Comm_rank { croak "FATAL: App_Comm_rank must be overridden by the Parallel::Application child\n"; }
sub App_Comm_size { croak "FATAL: App_Comm_size must be overridden by the Parallel::Application child\n"; }
sub App_Send { croak "FATAL: App_Send must be overridden by the Parallel::Application child\n"; }
sub App_Recv { croak "FATAL: App_Recv must be overridden by the Parallel::Application child\n"; }
sub APP_ANY_TAG { croak "FATAL: APP_ANY_TAG must be overridden by the Parallel::Application child\n";  }
sub APP_ANY_SOURCE { croak "FATAL: APP_ANY_SOURCE must be overridden by the Parallel::Application child\n";  }

#do not override this method
sub run {
    my $self = shift;
    
    $self->pre_init();
    $self->setup();
    $self->decide();
    $self->_App_Init();
    $self->post_init();
    $self->_run();
    $self->pre_finalize();
    $self->_App_Finalize();
    $self->post_finalize();
}
