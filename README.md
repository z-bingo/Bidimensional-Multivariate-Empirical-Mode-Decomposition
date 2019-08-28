## Bidimensional-Multivariate-Empirical-Mode-Decomposition
Matlab codes of Bidimensional Multivariate Empirical Mode Decomposition (BMEMD).

### Introduction
BMEMD is a bidimensional and multivariate version of original EMD, which is capable of processing multi-images, such as image fusion, texture analysis and so on. More details about the BMEMD can be referred in our paper [Bidimensional Multivariate Empirical Mode Decomposition with Applications in Multi-Scale Image Fusion](https://ieeexplore.ieee.org/document/8805082?source=authoralert).

### Requirements
- Image Processing Toolbox (installed in Matlab)
- gridfitdir (attanched in this repo., add the path to your Matlab environment)

### How to use these codes?
#### Files and Directories
- bmemd.m
  - main code of proposed BMEMD
- bmemd_fusion.m
  - its application on multi-images fusion, several images are provided at path `./IMG`
- Texture_Generate.m
  - the code to generate synthetic texture images in paper

#### Usages (take only decomposition as the example)
```matlab
x: [n, h, w], a non-int array
q: a cell of length Q, the number of IMFs, and each array in the cell 
   share the same size with x representing the corresponding IMF of x
```
1. `q=bmemd(x)`
2. `q=bmemd(x, ndir)`, here, `ndir` is the number of projections
3. `q=bmemd(x,ndir,stop_crtit)`, `stop_crtit` means stopping conditions, can be choosen from `'stop'` and  `'fix_h'`. If `'stop'`, the default parameter is `[0.01, 0.1, 0.01]`, ohtherwise, `fix_h=2`
4. `q=bmemd(x,ndir,'stop',[x1,x2,x3])`, `[x1,x2,x3]` is parameter of stop criteria
5. `q=bmemd(x,nir,'fix_h',fix_h)`

### Reference
[1] N. Rehman and D. P. Mandic,, "Multivariate empirical mode decomposition," Proc. R. Soc. A, vol.466, no. 2117, pp. 1291-1302, 2010.  
[2] T. Tanaka and D. P. Mandic, “Complex empirical mode decomposition,” IEEE Signal Process. Lett., vol. 14, no. 2, pp. 101–104, Feb. 2007.  
[3] G. Rilling, P. Flandrin, P. Gonalves, and J. M. Lilly, “Bivariate empirical mode decomposition,” IEEE Signal Process. Lett., vol. 14, no.12, pp. 936–939, Dec. 2007.  
[4] N. Rehman and D. P. Mandic, “Empirical mode decomposition for trivariate signals,” IEEE Trans. Signal Process., vol. 58, no. 3, pp. 1059–1068, Mar. 2010.  
[5] J. C. Nunes, S. Guyot, and E. Delechelle, “Texture analysis based on local analysis of the bidimensional empirical mode decomposition,” Mach. Vis. Appl., vol. 16, no. 3, pp. 177–188, May 2005.

### Citation
```tex
@article{Xia2019BMEMD, 
	author		= {Y. Xia and B. Zhang and W. Pei and D. P. Mandic}, 
	journal		= {IEEE Access},
	title		= {Bidimensional Multivariate Empirical Mode Decomposition With Applications in Multi-Scale Image Fusion},
	year		= {2019},
	volume		= {7}, 
	pages		= {114261-114270}, 
	doi 		= {10.1109/ACCESS.2019.2936030}, 
	ISSN		= {2169-3536}, 
	month		= {Dec.}
}
```
