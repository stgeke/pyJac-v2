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

#ifndef KERNEL_H
#define KERNEL_H

#include "mechanism.hpp"
#include "error_check.hpp"
#include "timer.hpp"
#include "jacobian_driver.hpp"
// undefine work size in main to avoid name collisions
#undef work_size


#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>



//! \\brief The base kernel class
class Kernel
{

public:
    Kernel();
    virtual ~Kernel();

    /*
    Resize kernel's working data to fit the given sizes.

    Parameters
    ----------
    problem_size : size_t
    	The total number of conditions to execute this kernel over
    work_size : size_t
    	The number of OpenMP threads to use.
    do_not_compile : bool
    	Unused -- incuded for consistent signatures.
    */

    void resize(size_t problem_size, size_t work_size, bool do_not_compile=false);

    //! \\brief Returns the ordered list of species names in the chemical model
    static std::vector<const char*> speciesNames()
    {
        return std::vector<const char*>(_species_names, _species_names + 9);
    }

    //! \\brief Returns the list of reaction strings in the chemical model
    static std::vector<const char*> reactionStrings()
    {
        return std::vector<const char*>(_rxn_strings, _rxn_strings + 28);
    }

    //! \\brief Return the data-ordering used in this kernel, either 'C' (row-major)
    //!        or 'F' (column-major)
    static const char* order()
    {
        return _order;
    }

    //! \\brief Return the number of species in the mode;
    static unsigned int numSpecies()
    {
        return _nsp;
    }

    //! \\brief Return the number of species in the mode;
    static unsigned int numReactions()
    {
        return _nrxn;
    }


    /** \\brief Returns the total amount of working memory required per-thermochemical
      *        state for this kernel, in bytes.
      *
      * \note  This includes vectorization considerations.
      */
    virtual const std::size_t requiredMemorySize() const = 0;

    void compile(){}
    void threadset(unsigned int num_threads);

protected:
    size_t per_run();
    size_t per_run(size_t problem_size);
    size_t this_run(size_t offset);


    // flags indicating initialization status, etc.
    bool initialized;
    bool compiled;

    // past run sizes
    size_t d_per_run; // store for device per-run size
    size_t problem_size;
    size_t max_per_run;
    size_t work_size;

    // info variables

    // species names
    static const char* _species_names[];
    // reaction strings
    static const char* _rxn_strings[];
    // data order
    static const char* _order;
    // number of species
    static const unsigned int _nsp;
    // number of reactions
    static const unsigned int _nrxn;


    /*
    Create the C kernel.

    Parameters
    ----------
    problem_size : size_t
    	The total number of conditions to execute this kernel over
    work_size : size_t
    	The number of OpenMP threads to use.
    */

    void init(size_t problem_size, size_t work_size);

    // memory initialization / release accomplished in sub-classes
    virtual void mem_init(size_t problem_size, size_t work_size) = 0;
    virtual void finalize_memory() = 0;
    void finalize();
};

// and subclass(es)
class JacobianKernel : public Kernel
{
protected:
    // declare device buffers
    double* d_P_arr;
    double* d_phi;
    double* d_jac;
    double* d_rwk;
#ifdef PINNED
    // declare temporary pointers to hold mapped addresses
    double* h_temp_d;
#endif
    void mem_init(size_t problem_size, size_t work_size);
public:
    /*
    Base constructor -- no arguments necessary, for use with Cython.
    JacobianKernel::resize() must be called before use.
    */
    JacobianKernel();

    /*
    Initializing constructor.

    Parameters
    ----------
    problem_size : size_t
    	The total number of conditions to execute this kernel over
    work_size : size_t
    	The number of OpenMP threads to use.
    do_not_compile : bool
    	Unused -- incuded for consistent signatures.
    */
    JacobianKernel(size_t problem_size, size_t work_size, bool do_not_compile=false);
    virtual ~JacobianKernel();

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
    void operator()(double* h_P_arr, double* h_phi, double* h_jac);
    void finalize_memory();
    const std::size_t requiredMemorySize() const;
};


#endif
