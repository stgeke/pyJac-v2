/*

A cog-templated skeleton for pyJac kernel execution

OpenCL code adapted from:
    Based on https://www.olcf.ornl.gov/tutorials/opencl-vector-addition/
    and https://www.fixstars.com/en/opencl/book/OpenCLProgrammingBook/calling-the-kernel/

(C) Nicholas Curtis - 2018

Global declarations for Cog:
    - codegen: path to a serialized CallgenResult instance
    that may be loaded to generate this file
*/
#include <mpi.h>
#include "jacobian_main.hpp"
#include "read_initial_conditions.hpp"
#include "write_data.hpp"
#include "memcpy_2d.hpp"

Kernel::Kernel():
    initialized(false),
    compiled(false),
    problem_size(0),
    work_size(0)
{
    this->compiled=true;
}

Kernel::~Kernel()
{
}

void Kernel::resize(size_t problem_size, size_t work_size, bool do_not_compile)
{
    if (do_not_compile and !this->compiled)
    {
        // Assume that OpenCL kernel has previously been compiled (e.g., via a previous kernel call)
        this->compiled = true;
    }
    if (this->initialized && ((work_size != this->work_size) || (problem_size != this->problem_size)))
    {
        this->finalize_memory();
        this->init(problem_size, work_size);
    }
    else if (!this->initialized)
    {
        this->init(problem_size, work_size);
    }
}

size_t Kernel::per_run()
{
    return this->max_per_run < this->problem_size ? this->max_per_run : this->problem_size;
}

size_t Kernel::per_run(size_t problem_size)
{
    return this->max_per_run < problem_size ? this->max_per_run : problem_size;
}

size_t Kernel::this_run(size_t offset)
{
    size_t per_run = this->per_run();
    return this->problem_size - offset < per_run ? this->problem_size - offset : per_run;
}

/*
  Resets the program for a change in number of threads
*/
void Kernel::threadset(unsigned int num_threads)
{
    // get maximum allowed threads
    unsigned int max_threads = omp_get_max_threads();
    // check that # of threads < max allowed
    cassert(num_threads <= max_threads, "Can't use more than the maximum allowed threads by OpenMP.");
    // set number of threads
    omp_set_num_threads(num_threads);
    // and store
    this->work_size = num_threads;
}


const char* Kernel::_species_names[] = { "H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR" };
const char* Kernel::_rxn_strings[] = { "ThreeBodyReaction: 2 O + M <=> O2 + M", "ThreeBodyReaction: O + H + M <=> OH + M", "ElementaryReaction: O + H2 <=> H + OH", "ElementaryReaction: O + HO2 <=> OH + O2", "ElementaryReaction: O + H2O2 <=> OH + HO2", "ThreeBodyReaction: H + O2 + M <=> HO2 + M", "ElementaryReaction: H + 2 O2 <=> HO2 + O2", "ElementaryReaction: H + O2 + H2O <=> HO2 + H2O", "ElementaryReaction: H + O2 + AR <=> HO2 + AR", "ElementaryReaction: H + O2 <=> O + OH", "ThreeBodyReaction: 2 H + M <=> H2 + M", "ElementaryReaction: 2 H + H2 <=> 2 H2", "ElementaryReaction: 2 H + H2O <=> H2 + H2O", "ThreeBodyReaction: H + OH + M <=> H2O + M", "ElementaryReaction: H + HO2 <=> O + H2O", "ElementaryReaction: H + HO2 <=> O2 + H2", "ElementaryReaction: H + HO2 <=> 2 OH", "ElementaryReaction: H + H2O2 <=> HO2 + H2", "ElementaryReaction: H + H2O2 <=> OH + H2O", "ElementaryReaction: OH + H2 <=> H + H2O", "FalloffReaction: 2 OH (+M) <=> H2O2 (+M)", "ElementaryReaction: 2 OH <=> O + H2O", "ElementaryReaction: OH + HO2 <=> O2 + H2O", "ElementaryReaction: OH + H2O2 <=> HO2 + H2O", "ElementaryReaction: OH + H2O2 <=> HO2 + H2O", "ElementaryReaction: 2 HO2 <=> O2 + H2O2", "ElementaryReaction: 2 HO2 <=> O2 + H2O2", "ElementaryReaction: OH + HO2 <=> O2 + H2O" };
const char* Kernel::_order = "C";
const unsigned int Kernel::_nsp = 9;
const unsigned int Kernel::_nrxn = 28;

/*
Create C Kernel

Parameters
----------
problem_size : size_t
	The total number of conditions to execute this kernel over
work_size : size_t
	The number of OpenMP threads to use.
*/

void Kernel::init(size_t problem_size, size_t work_size)
{
    this->threadset(work_size);
    this->mem_init(problem_size, work_size);
    // mark initialized
    this->initialized = true;
    this->problem_size = problem_size;
    this->work_size = work_size;
}



/*
   Destroys the kernel object, will not be usable again
*/
void Kernel::finalize()
{
    if(this->initialized)
    {

        // mark deinit
        this->compiled = false;
        this->initialized = false;
    }
}


/*


Parameters
----------
problem_size : size_t
	The total number of conditions to execute this kernel over
work_size : size_t
	The number of OpenMP threads to use.
*/


void JacobianKernel::mem_init(size_t problem_size, size_t work_size)
{
    size_t per_run = this->d_per_run = this->per_run(problem_size);
    /* If we've run out of constant memory space, we will place converted
       global constant here */

    /* Alloc buffers */

    this->d_P_arr = (double*)malloc(per_run * sizeof(double));
    cassert(this->d_P_arr != NULL, "malloc failed");
    this->d_phi = (double*)malloc(10 * per_run * sizeof(double));
    cassert(this->d_phi != NULL, "malloc failed");
    this->d_jac = (double*)malloc(100 * per_run * sizeof(double));
    cassert(this->d_jac != NULL, "malloc failed");
    this->d_rwk = (double*)malloc(work_size*371 * sizeof(double));
    cassert(this->d_rwk != NULL, "malloc failed");
    /* and memset to zero */

    memset(this->d_P_arr, 0, per_run * sizeof(double));
    memset(this->d_phi, 0, 10 * per_run * sizeof(double));
    memset(this->d_jac, 0, 100 * per_run * sizeof(double));
    memset(this->d_rwk, 0, work_size*371 * sizeof(double));

    /* Transfer host constants here (if any), as we only need to do so once */

}

void JacobianKernel::finalize_memory()
{
    /* Free Memory */
    free(this->d_P_arr);
    free(this->d_phi);
    free(this->d_jac);
    free(this->d_rwk);
}

const std::size_t JacobianKernel::requiredMemorySize() const
{
    return work_size*371 * sizeof(double);
}

/*
Default constructor (requires resize before use)

Parameters
----------

*/

JacobianKernel::JacobianKernel()
{
    this->max_per_run = 21474836;
}


/*
Initializing contstructor

Parameters
----------
problem_size : size_t
	The total number of conditions to execute this kernel over
work_size : size_t
	The number of OpenMP threads to use.
do_not_compile : bool
	Unused -- incuded for consistent signatures.
*/

JacobianKernel::JacobianKernel(size_t problem_size, size_t work_size, bool do_not_compile)
{
    this->max_per_run = 21474836;
    this->compiled = do_not_compile;
this->init(problem_size, work_size);
}

/*
Execute the C kernel 'jacobian'

Parameters
----------
P_arr : double
	The array of pressures.
phi : double
	The state vector
jac : double
	The Jacobian of the time-rate of change of the state vector
*/

JacobianKernel::~JacobianKernel()
{
  this->finalize_memory();
  this->finalize();
}
void JacobianKernel::operator()(double* h_P_arr, double* h_phi, double* h_jac)
{
    cassert(this->initialized, "Must initialize kernel (e.g., via resize()) before use.");
    size_t per_run = this->d_per_run = this->per_run();

    for (size_t offset = 0; offset < this->problem_size; offset += per_run)
    {
        size_t this_run = this->this_run(offset);

        // Memory Transfers into the kernel, if any
        memcpy(this->d_P_arr, &h_P_arr[offset * 1], this_run * sizeof(double));
        memcpy(this->d_phi, &h_phi[offset * 10], 10 * this_run * sizeof(double));
        // run kernel
        jacobian_driver(this->d_per_run, this->d_P_arr, this->d_phi, this->d_jac, this->d_rwk);

        // Memory Transfers out
        memcpy(&h_jac[offset * 100], this->d_jac, 100 * this_run * sizeof(double));
    }
}


int main(int argc, char* argv[])
{

    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
    //check args
    cassert(argc >= 3, "Missing arguments...");

    //arglist is:
    //#0 - the program name
    //#1 - the problem size
    //#2 - the number of cores / threads [CPU/Accelerator] or number GPU blocks [GPU only]
    //#3 - whether to compile

    size_t problem_size = atoi(argv[1])/world_size;
    size_t work_size = atoi(argv[2]);
    int compile = 1;
    if (argc >= 4)
        compile = atoi(argv[3]);

    double* h_P_arr_local;
    double* h_phi_local;
    double* h_jac_local;
    h_P_arr_local = (double*)malloc(problem_size * sizeof(double));
    cassert(h_P_arr_local != NULL, "malloc failed");
    h_phi_local = (double*)malloc(10 * problem_size * sizeof(double));
    cassert(h_phi_local != NULL, "malloc failed");
    h_jac_local = (double*)malloc(100 * problem_size * sizeof(double));
    cassert(h_jac_local != NULL, "malloc failed");

    //read input data
    read_initial_conditions("data.bin", problem_size, h_P_arr_local, h_phi_local, 'C');
    JacobianKernel kernel;
    //first compile to binary
    double compilation_time = -1;
    if (compile)
    {
        StartTimer();
        kernel.compile();
        compilation_time = GetTimer();
    }

    StartTimer();
    kernel.resize(problem_size, work_size, !compile);
    double setup_time = GetTimer();
    MPI_Barrier(MPI_COMM_WORLD);
    StartTimer();
    kernel(h_P_arr_local, h_phi_local, h_jac_local);
    MPI_Barrier(MPI_COMM_WORLD);
    double runtime = GetTimer();

    problem_size = world_size*problem_size;
    if(world_rank == 0)
    printf("%zu,%.15le,%.15le,%.15le,%.15le\n", problem_size, compilation_time/1e3,
                setup_time/1e3, runtime/1e3, (double)problem_size/(runtime/1e3));

    // write output to file if supplied

    // local frees
    free(h_P_arr_local);
    free(h_phi_local);
    free(h_jac_local);
    
    MPI_Finalize();
    return 0;
}
