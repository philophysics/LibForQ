LibForQ: A Fortran Library for Quantum Information Science

In this library we will continuously add and improve Fortran code to perform several numerical tasks one frequently 
needs when working in quantum information science. Among the functionalities already implemented, some examples are:
- Generators for random numbers, probability vectors, unitary matrices, state vectors, and density matrices;
- Trace, partial trace, and partial transpose;
- Entanglement, discord, and coherence quantifiers;
- Pauli group (PG), Generalized Gell Mann Matrices (GGMM), and Bloch vector and correlation matrix with GGMM;
- Some other matrix functions are also included, some examples are:
  > Norms, inner products, and distance and distinguishability measures;
  > Purity, entropies, information measures, and mutual information;
  > Some popular, and not so popular, quantum states;
  > Array display, Ginibre and identity matrices, adjoint, Kronecker product (KP), KP of n elements of the PG,
    Gram-Schmidt orthogonalization, projector, outer product, random permutation for a vector components, etc.
I think it is better to take a look at the code and see all functions included. The variables used and needed are 
explained there.

One can use this library simply by copying all the files to the main program folder and compiling them with:
  gfortran -lblas -llapack *.f90
To run the executable file a.out generated just use:
  ./a.out

Another, perhaps more convenient, manner of using the code is by creating a library. 
For that purpose you may follow the simple steps below:
1) Extract compressed file
2) Go to the associated folder
3) Create a library with the commands:
  gfortran -O3 -c rng_mt.f90
  gfortran -O3 -c *.f90
and
  ar cr libforq.a *.o *.mod 
To compile your main program using this library, copy libforq.a to your program's folder and use the command: 
  gfortran -lblas -llapack libforq.a main.f90
Even better, you can also copy the library to your computer's libraries folder, e.g. with:
  sudo cp libforq.a /usr/local/lib
and use it, anywhere in your computer, via
  gfortran -lblas -llapack -lforq main.f90
