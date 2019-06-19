# HArDCore2D
Hybrid Arbitrary Degree::Core 2D - Library to implement schemes with edge and cell polynomial unknowns on generic 2D polygonal meshes.

The complete documentation is available here:

https://jdroniou.github.io/HArDCore2D/

The purpose of HArD::Core2D is to provide easy-to-use tools to code hybrid schemes, such as the Hybrid High-Order method. The data structure is described using intuitive classes that expose natural functions we use in the mathematical description of the scheme. For example, each mesh element is a member of the class 'Cell', that gives access to its diameter, the list of its edges (themselves members of the class 'Edge' that describe the geometrical features of the edge), etc. Functions are also provided to compute the key elements usually required to implement hybrid schemes, such as mass matrices of local basis functions, stiffness matrices, etc. The approach adopted is that of a compromise between readibility/usability and efficiency. 

So, for example, when creating a mass matrix, the library requires the user to first compute the quadrature nodes and weights, then compute the basis functions at these nodes, and then assemble the mass matrix. This ensures a local control on the required degree of exactness of the quadrature rule, and also that basis functions are not evaluated several times at the same nodes (once computed and stored locally, the values at the quadrature nodes can be re-used several times). Each of these steps is however concentrated in one line of code, so the assembly of the mass matrix described above is actually done in three lines:

```
std::vector<HybridCore::qrule> quadT = hho.cell_qrule(iT, doe);<br>
std::vector<Eigen::ArrayXd> phi_quadT = hho.basis_quad('T', iT, quadT, nb_basis_functions);<br>
Eigen::MatrixXd MTT = hho.gram_matrix(phi_quadT, phi_quadT, nb_basis_functions, nb_basis_functions, quadT, true);
```

More details and examples are provided in the documentation.

This library has been developed with the direct help and indirect advice of several people. Many thanks to them: Daniel Anderson, Lachlan Grose, Tom Lemaitre, Daniele Di Pietro, Lorenzo Botti.


