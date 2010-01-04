#-----Database Setup
DBI:Pg #interface type to database
dbname:makerweb #database name
host: #host on which database is found
port: #port on host to access database
username:cholt #username to connect to database
password:Annotate_this #password to connect to database

#-----Communication Options
admin_email:carson.holt@genetics.utah.edu #an e-mail address to send error and status information to
smtp_server:m2.genetics.utah.edu #outgoing e-mail. Required for users to receive notification of finished jobs.

#-----MAKER Server Specific Options
use_login:1 #whether to require login to access the web interface, 1 = yes, 0 = no
allow_guest:1 #enable guest accounts on the server, 1 = yes, 0 = no
allow_register:1 #allow users to register themselves vs. manually by admin, 1 = yes, 0 = no
MPI:0 #use mpi_maker instead of maker
max_cpus:2 #maximum number of cpus that can be dedicated to all MAKER jobs
job_cpus:1 #maximum number of cpus that can be used by a single MAKER job
max_submit_user:2000000 #maximum submission size for registered users (0 to disable limit)
max_submit_guest:200000 #maximum submission size for guest users (0 to disable limit)
persist_time_user:336 #time results persist for registered users, in hours (0 to disable limit)
persist_time_guest:72 #time results persist for guest users, in hours (0 to disable limit)
inactive_user:0 #time user account can be inactive before disabling, in days (0 to disable limit)
inactive_guest:14 #time guest account can be inactive before disabling, in days (0 to disable limit)
data_dir:/home/apache/MWS/data
cgi_dir:/data/var/www/cgi-bin/MWAS/ #web accesible directory where web interface CGI content is stored
html_dir:/data/var/www/html/MWAS/ #web accesible directory where web interface HTML conent is stored
APOLLO_ROOT:/usr/local/gmod/apollo #base directory for Apollo installation.  Used for building webstart of Apollo.
