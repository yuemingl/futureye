<a href='http://code.google.com/p/scala-fem/'>Click here for ScalaFEM</a> -- A mathematical-look wrapper project of FuturEye written in Scala.

## Introduction ##

FuturEye is a Finite Element Methods Toolkit written in Java. It provides a concise, natural, easily understandable and mathematically appealing programming interface. It is also designed to solve inverse problems based on FEM. FuturEye, comes from the application in optical tomography achieved by solving inverse problems of the differential equation.

The essential components of FEM are abstracted out, such as nodes,
elements, meshes, degrees of freedom and shape function etc. The data and operations of
these classes are encapsulated together properly. The classes that
different from other existing object-oriented FEM softwares or
libraries are function classes. The behavior of the function classes
in Futureye is very similar to that in mathematical context. For
example algebra of functions, function derivatives and composition
of functions. Especially in FEM environment, shape functions,
Jacobin of coordinate transforms and numerical integration are all
based on the function classes. This feature leads to a more close
integration between theory and computer implementation.


FuturEye is designed to solve 1D,2D and 3D partial differential
equations(PDE) of scalar and/or vector valued unknown functions. It is motivated by solving inverse problems of partial differential equations in application of optical tomography. In order to solve inverse problems, usually some
forward problems must be solved first and many exterior data
manipulations should be performed during the solving processes.
There are many classes defined for those data operations. However,
the data processes are complicated in actual applications, we can
not write down all the tools for manipulating data. The design of
the basic classes allows operations to all aspect of data structure
directly or indirectly in an easily understanding way. This is
important for users who need to write their own operations or
algorithms on the bases of data structure in FuturEye. Some other
existing FEM softwares or libraries may over encapsulate the data
and operations for such inverse applications.

This toolkit can be used for various purposes:
  * Teaching: The feature of close relation to mathematical theory of FEM will help a student to understand basic FEM concepts, e.g. shape functions, the Jacobian and assembly process.
  * Research: FuturEye helps researchers quickly develop and test their models, experiment with data and algorithms. e.g. new equations, finite elements and solution methods without concerning too much about basic components in FEM programming.
  * Engineering: The performance and efficiency may be unsatisfactory for real applications,if a finite element class defined in a mathematical manner is without optimization.Thanks to the interface conception in Java, we can implement the same interface in many different ways, thus a carefully optimized finite element class can be used in applications with a huge number of elements.

## Futureye Users Distribution Around the World ##

<img width='540' height='288' src='https://lh6.googleusercontent.com/-iA032dYtHkU/UoMEbkf1gaI/AAAAAAAAAxY/kXLVPuOj9nE/w673-h415-no/users.png' />


## Object Recognition and Analysis from Images ##

<table width='970' height='290'>
<tr>
<td>
<img width='320' height='144' src='https://lh4.googleusercontent.com/-9keZTkOsMtM/UH9wUFbAILI/AAAAAAAAAvI/C7sMgqNQjNE/w484-h224-n-k/car1-1.png' />
</td>
<td>
<img width='320' height='144' src='https://lh3.googleusercontent.com/-aniDC0jfDCw/UH9wUFUXuHI/AAAAAAAAAvE/n8PnktxrNqE/w484-h224-n-k/car1-2.png' />
</td>
<td>Iphone demos :)</td>
</tr>
<tr>
<td>
<img width='320' height='144' src='https://lh3.googleusercontent.com/-GbwXVit_TlM/UH9wURc60TI/AAAAAAAAAvM/loMT-rXAxP8/w593-h273-n-k/hand1-1.png' />
</td>
<td>
<img width='320' height='144' src='https://lh4.googleusercontent.com/-VRmMoET0F_I/UH9wUrsmoDI/AAAAAAAAAvU/iGaeRwMfYB4/w594-h273-n-k/hand1-2.png' />
</td>
<td><img width='320' height='144' src='https://lh6.googleusercontent.com/-0n5nOmbANv0/UH9wVKoOtrI/AAAAAAAAAvc/clPAbAe2HBw/w654-h302-n-k/hand1-3.png' /></td>
</tr>
</table>

## Numerical Optimization for Inverse Problems & Brain Imaging ##
<table width='1024' height='360'>
<tr>
<td>
<img width='650' height='300' src='https://lh4.googleusercontent.com/-MsgxFIrber8/T-HaVnasgxI/AAAAAAAAAt4/oeGSx-y1L_U/s800/Language1.png' />
</td>
<td>
<img width='250' height='250' src='https://lh4.googleusercontent.com/-c1dcIyooXao/UGcf894N_KI/AAAAAAAAAu4/gu-GELAnBV0/s466/Brain+Imaging.png' />
</td>
</tr>
</table>


## Plane Elastic ##

<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh4.googleusercontent.com/-6BKnBWr-xUQ/TsaxTCHVeMI/AAAAAAAAAr4/LZeaDSZ_V_E/s288/elasticDamInit.png' />
</td>
<td>
<img src='https://lh3.googleusercontent.com/-CGqi7mLpOeo/TsaxTPkESVI/AAAAAAAAArk/K3mf8ixVBuU/s288/elasticDamRlt.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/-gvjWG6MvFXc/TsaxTIsLuQI/AAAAAAAAAro/8Cze5IQ2lf8/s288/elasticHoleInit.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/-IlcKQ47l-0A/TsaxThv5QsI/AAAAAAAAAr0/gUdNXBk9x_U/s288/elasticHoleRlt.png' />
</td>

</tr>
</table>

## Stokes Equation 1 ##

<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh5.googleusercontent.com/_Cil2MFH7iLM/TYV8oRWm0gI/AAAAAAAAAFw/ZlVl5dQGl04/s288/cylinder1.png' />
</td>
<td>
<img src='https://lh5.googleusercontent.com/_Cil2MFH7iLM/TYV8oxRYarI/AAAAAAAAAF0/P4SrH1x-Gqw/s288/cylinder2.png' />
</td>
<td>
<img src='https://lh4.googleusercontent.com/_Cil2MFH7iLM/TYV8o-ST-RI/AAAAAAAAAF4/djcZzGrJ2Fc/s288/cylinder4.png' />
</td>
</tr>
</table>

## Stokes Equation 2 ##

<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh3.googleusercontent.com/_Cil2MFH7iLM/TbpB8RvD_iI/AAAAAAAAAH0/n0YXrevRCoc/s288/box.png' />
</td>
<td>
<img src='https://lh4.googleusercontent.com/_Cil2MFH7iLM/TY_1viEOX0I/AAAAAAAAAGg/7TTar7IHuMk/s288/u_shape1.png' />
</td>
</tr>
</table>


## Adaptive FEM ##

<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TN19lZMjC9I/AAAAAAAAABw/T4uSGqOSHPE/s288/refine_r2.png.jpg' />
</td>
<td>
<img src='https://lh4.googleusercontent.com/_Cil2MFH7iLM/TN19mTKpbHI/AAAAAAAAAB8/T3DlYmqNyfQ/s288/refine_t2.png.jpg' />
</td>
<td>
<img src='https://lh5.googleusercontent.com/_Cil2MFH7iLM/TN19qJvDOCI/AAAAAAAAACk/elTca0hdmYE/s288/adaptive10.png.jpg' />
</td>
</tr>
</table>

## Convection diffusion ##

<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh3.googleusercontent.com/_Cil2MFH7iLM/TXqXepCO6fI/AAAAAAAAAEY/Ilapk8kOAgw/s288/ConDiff_t1.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TXqXeoiJgRI/AAAAAAAAAEU/R5vFs_t_5gY/s288/ConDiff_t2.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TXqXeiFB2CI/AAAAAAAAAEg/6A6f4AG9pNk/s288/ConDiff_t3.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TXqXfDWzrNI/AAAAAAAAAEc/lap5zRFPW74/s288/ConDiff_t4.png' />
</td>
<td>
<img src='https://lh3.googleusercontent.com/_Cil2MFH7iLM/TXqXfQ2XJUI/AAAAAAAAAEk/TdmfpSJ9GK0/s288/ConDiff_t5.png' />
</td>
</tr>
</table>


## Wave Equation ##

<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TXqZZh11ejI/AAAAAAAAAEs/sAscKtUytYE/s288/Wave0.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TXqZZwfs8xI/AAAAAAAAAEw/99rYRxDobIA/s288/Wave1.png' />
</td>
<td>
<img src='https://lh3.googleusercontent.com/_Cil2MFH7iLM/TXqZZtP2_aI/AAAAAAAAAE0/d7ozij8ZZwQ/s288/Wave2.png' />
</td>
<td>
<img src='https://lh4.googleusercontent.com/_Cil2MFH7iLM/TXqZaGmBL_I/AAAAAAAAAE4/ohy35upbgN0/s288/Wave3.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TXqZaEUrOqI/AAAAAAAAAFA/2Fuu7K70kfw/s288/Wave4.png' />
</td>
<td>
<img src='https://lh3.googleusercontent.com/_Cil2MFH7iLM/TXqZaBNIi0I/AAAAAAAAAE8/KaAUOme9wGE/s288/Wave5.png' />
</td>
<td>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TXqZahmRglI/AAAAAAAAAFE/LI_y3aV8hVQ/s288/Wave6.png' />
</td>
</tr>
</table>


## 3D Laplace ##

<table width='1024' height='288'>
<tr>
<td align='center'>
<img src='https://lh6.googleusercontent.com/_Cil2MFH7iLM/TZzHmMZN5pI/AAAAAAAAAG0/ZnHrd54OmiQ/s288/block2.png' />
</td>
<td align='center'>
<img src='https://lh3.googleusercontent.com/-id9ixuwimOE/TtPWfXZ2SUI/AAAAAAAAAsY/pUSBplUdu6U/s288/Laplace3D_Hex.png' />
</td>
<td align='center'>
<img src='https://lh3.googleusercontent.com/-FGrHJs2vCEQ/Ty3JeLAVzII/AAAAAAAAAss/HRj7bsKcfdU/s288/3D%2520model.png' />
</td>
</tr>
<tr>
<td align='center'>
Tetrahedron element<br>
</td>
<td align='center'>
Hexahedron element<br>
</td>
<td align='center'>
Hexahedron element<br>
</td>
</tr>
</table>


## Benchmark of Navier-Stokes Equations ##
<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh3.googleusercontent.com/-1dVIzXpwgac/TsafQjWW3WI/AAAAAAAAAqY/7qRiEKDEmcA/s400/NS1.png' />
</td>
<td>
<img src='https://lh3.googleusercontent.com/-flx640NVIf4/TsafQQUHOZI/AAAAAAAAAqQ/Y2zq_GU57gQ/s400/NS2.png' />
</td>
</tr>
<tr>
<td>
<img src='https://lh5.googleusercontent.com/--kWi6uDzZbY/TsafPwE_63I/AAAAAAAAAqA/W3XSIDxyFkk/s400/NS3.png' />
</td>
<td>
<img src='https://lh3.googleusercontent.com/-smYdom1bckU/TsafQFfkKVI/AAAAAAAAAqI/1xJK6MY5Kec/s400/NS4.png' />
</td>
</tr>
<tr>
<td>Stream line of Lid-driven cavity (Re=1000):</td>
<td><img src='https://lh5.googleusercontent.com/-MgfdT9KqmRI/Tsah4Xd15WI/AAAAAAAAAq0/D32uVqKnmTo/s288/NSBox.png' />
</td>
</tr>
</table>

## Inverse Problem: Near-Infrared Diffusion Optical Tomography (DOT) ##
<table width='1024' height='288'>
<tr>
<td>
<img src='https://lh3.googleusercontent.com/-R98pYWlTVRE/TsavryY5GnI/AAAAAAAAArQ/iIam-vbDVQc/s288/DOT1.png' />
</td>
<td>
<img src='https://lh5.googleusercontent.com/-QTCz_HJtGQY/Tsavroh0U7I/AAAAAAAAArA/GlEzcMR7-9U/s288/DOT2.png' />
</td>
<td>
<img src='https://lh4.googleusercontent.com/-XFbSayVIEi4/Tsavrp2P3FI/AAAAAAAAArE/74SF44b7nFE/s288/DOT3.png' />
</td>
<td>
<img src='https://lh4.googleusercontent.com/-4OgNXhubj7Q/TsavsB1Nb7I/AAAAAAAAArY/Pgc8d6oQd-s/s288/DOT4.png' />
</td>
</tr>
</table>