# finite_difference

Computes finite difference matrices for the first and second derivative up to sixth order, including compact finite-difference schemes such as the fourth-order Padé scheme and sixth-order Lele scheme (S. K. Lele. Compact finite difference schemes with spectral-like resolution. Journal of Computational Physics, 103(1):16–42, November 1992.

## Example

The function shown in the example files is given by:

<img src="https://render.githubusercontent.com/render/math?math=f(x,y) = \sin x \ \sin y">

The partial derivatives are given by:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial f}{\partial x} = \cos x \ \sin y">

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial f}{\partial y} = \sin x \ \cos y">

The images below show the absolute error in log scale for both the partial derivatives in x and y direction.

### Second order

<img src="images/partialx_2.png" width="425"/> <img src="images/partialy_2.png" width="425"/>

### Fourth order

<img src="images/partialx_4.png" width="425"/> <img src="images/partialy_4.png" width="425"/>

### Sixth order

<img src="images/partialx_6.png" width="425"/> <img src="images/partialy_6.png" width="425"/>

### Padé fourth order

<img src="images/partialx_pade4.png" width="425"/> <img src="images/partialy_pade4.png" width="425"/>

### Lele sixth order

<img src="images/partialx_lele6.png" width="425"/> <img src="images/partialy_lele6.png" width="425"/>
