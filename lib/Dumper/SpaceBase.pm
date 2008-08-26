#<top> ------------------------------------------------------------------------
#     <Name>       SpaceBase.pm                    </Name> 
#     <class>      
#                  SpaceBase             </class>       
#     <Author>     M Hadi Islam            </Author>
#     <email>      hadi@genetics.utah.edu  </email>
#     <Does>       Conversion between space base
#                  (natural begin) to base base 
#                  number system  </Does>
#</top> ----------------------------------------------------------------------



package SpaceBase;

use Class::Struct;
use Data::Dumper;
use strict "vars";
use strict "refs";
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#---------------------- SUBS ---------------------------------------------------
#-------------------------------------------------------------------------------

sub tospace0base{
    my ($nbeg,$nend)=@_;
  
    if($nbeg < $nend){
	$nbeg--;
    }
    if($nbeg >  $nend){
	$nend--;
    }
    if($nbeg == $nend){
	$nbeg--;
    }

    return $nbeg,$nend;
}

#----------------------------------------------------------------------
sub tobase1base{
    my ($nbeg,$nend)=@_;
  
    if($nbeg < $nend){
	$nbeg++;
    }

    if($nbeg >  $nend){
	$nend++;
    }
    
    return $nbeg,$nend;
}



#<end>---------------------------------------------------------------------- 
1;
#</end>--------------------------------------------------------------------- 
