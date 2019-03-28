# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <mpi.h>
# include <math.h>
# include "brownian_particle.h"
# include "Random_Numbers.h"
# include "stats_and_histogram.h"

using namespace std;

int main ( int argc, char *argv[] );
void timestamp ( );

int main ( int argc, char *argv[] )


{
  // MPI Variables
  MPI_Status status;
  MPI_Request request = MPI_REQUEST_NULL;
  int ierr;
  int p;
  int id;
  int master;

  // Brownian Particle Parameters
  double kT = 7.0
  , eta = 5.0
  , m = 11.0
  , dt = 0.001;

  int tau = (int) ceil(5*m/eta)
  , nsteps = 1000;

  double dT_out = tau*dt;

  brownian_particle b = brownian_particle(nsteps, dt, eta, kT, m, 0.0, 0.0, dT_out);

  // Brownian Statistics
  int nrepeats = 50000000
  , nrepeats_local = 1
  , nrepeats_0 = 1;

  stat_list s_l = e_x | e_v | e_x2 | e_v2 | e_xv;

  brownian_averages b_stat = brownian_averages(b, s_l);
  
  // For receiving information
  int num_samples_received = 0;
  double stats_received[nsteps];

  // File Writing
  char fname[] = "/home/btreece/Programs/BASIC_MPI/OUTPUT_50000000.txt";
  FILE *fp;

//-------------------//
//  Initialize MPI.  //
//-------------------//
  ierr = MPI_Init ( &argc, &argv );

  if ( ierr != 0 )
  {
    printf("\nSUM - Fatal error!\n  MPI_Init returned ierr = %i\n", ierr);
    exit ( 1 );
  }

//  Get the number of processes.
  MPI_Comm_size ( MPI_COMM_WORLD, &p );

//  Determine the rank of this process.
  MPI_Comm_rank ( MPI_COMM_WORLD, &id );

//  State the current time and prepare the work load
  if ( id == 0 )
  {
    timestamp ( );
    nrepeats_local = (int) nrepeats / p;
    nrepeats_0 = nrepeats - (p-1)*nrepeats_local;
  }

//  Set random seeds and check them
  srand(time(0)+id);
  printf("Node Id %i, Random Number %.4f\n", id, rand_uniform_0_1());

//  The master process broadcasts the work load 
//  to all the other processes.
  master = 0;
  MPI_Bcast ( &nrepeats_local, 1, MPI_INT, master, MPI_COMM_WORLD );

//  Give the master node its work load
  if ( id == 0 )
  { nrepeats_local = nrepeats_0; }

//  Hard syncronize
  MPI_Barrier(MPI_COMM_WORLD);

//  Each process runs it's trajectories
  printf("Node %i: workload %i\n", id, nrepeats_local);
  b_stat.add_samples(nrepeats_local, 1.0);

  MPI_Barrier(MPI_COMM_WORLD);
//  Each worker process sends its arrays back to the master process.
  if ( id != 0 )
  {
    MPI_Isend ( &b_stat.num_samples, 1, MPI_INT, master, 0, MPI_COMM_WORLD, &request );
    if (b_stat.xbar)
    {
      MPI_Isend ( b_stat.xbar, nsteps, MPI_DOUBLE, master, 1, MPI_COMM_WORLD, &request );
    }
    if (b_stat.x2bar)
    {
      MPI_Isend ( b_stat.x2bar, nsteps, MPI_DOUBLE, master, 2, MPI_COMM_WORLD, &request );
    }
    if (b_stat.vbar)
    {
      MPI_Isend ( b_stat.vbar, nsteps, MPI_DOUBLE, master, 3, MPI_COMM_WORLD, &request );
    }
    if (b_stat.v2bar)
    {
      MPI_Isend ( b_stat.v2bar, nsteps, MPI_DOUBLE, master, 4, MPI_COMM_WORLD, &request );
    }
    if (b_stat.xvbar)
    {
      MPI_Isend ( b_stat.xvbar, nsteps, MPI_DOUBLE, master, 5, MPI_COMM_WORLD, &request );
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else 
  {
    for ( int i = 1; i < p; i++ ) 
    {
      MPI_Recv ( &num_samples_received, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
      printf("Node %i, samples sent: %i\n", i, num_samples_received);
      if (b_stat.xbar)
      {
        MPI_Recv ( &stats_received, nsteps, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status );
        for (int j = 0; j < nsteps; j++)
        { b_stat.xbar[j] = ( b_stat.xbar[j]*b_stat.num_samples + stats_received[j]*num_samples_received ) / ( b_stat.num_samples + num_samples_received );}
      }
      if (b_stat.x2bar)
      {
        MPI_Recv ( &stats_received, nsteps, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status );
        for (int j = 0; j < nsteps; j++)
        { b_stat.x2bar[j] = ( b_stat.x2bar[j]*b_stat.num_samples + stats_received[j]*num_samples_received ) / ( b_stat.num_samples + num_samples_received );}
      }
      if (b_stat.vbar)
      {
        MPI_Recv ( &stats_received, nsteps, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status );
        for (int j = 0; j < nsteps; j++)
        { b_stat.vbar[j] = ( b_stat.vbar[j]*b_stat.num_samples + stats_received[j]*num_samples_received ) / ( b_stat.num_samples + num_samples_received );}
      }
      if (b_stat.v2bar)
      {
        MPI_Recv ( &stats_received, nsteps, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &status );
        for (int j = 0; j < nsteps; j++)
        { b_stat.v2bar[j] = ( b_stat.v2bar[j]*b_stat.num_samples + stats_received[j]*num_samples_received ) / ( b_stat.num_samples + num_samples_received );}
      }
      if (b_stat.xvbar)
      {
        MPI_Recv ( &stats_received, nsteps, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &status );
        for (int j = 0; j < nsteps; j++)
        { b_stat.xvbar[j] = ( b_stat.xvbar[j]*b_stat.num_samples + stats_received[j]*num_samples_received ) / ( b_stat.num_samples + num_samples_received );}
      }
      b_stat.num_samples += num_samples_received;
      printf("Master Node, %i samples\n", b_stat.num_samples);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Master Process Is Writing To File.\n\n");
    fp = fopen(fname, "a");
    b_stat.print_data_to_file(fp);
    fclose(fp);
  }

//  Terminate MPI.
  MPI_Finalize ( );

//  Terminate.
  if ( id == 0 )
  {
    cout << "\n";         
    cout << "Master process:\n";       
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );      
  }
  return 0;

}

void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n\n";

  return;
# undef TIME_SIZE
}
