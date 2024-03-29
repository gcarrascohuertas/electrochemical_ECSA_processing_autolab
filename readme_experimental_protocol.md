# About

Here I provided a short-guide of experimental acquisition data. I summarized the main things you need to take into account but [gcarrascohuertas]( https://github.com/gcarrascohuertas) is not responsible for the bad use of this guide and its consequences over materials and people involved. Please, be careful and document all steps you perform. 

# Electrochemical setup

First you need following materials in order to perform experiments:

- Three-electrode cell

    - Platinum mesh as counter electrode (CE)). Here I used  a [Planium mesh from Goodfellow](http://www.goodfellow.com/A/Platinum-Mesh.html)
    - Ag/AgCl as reference electrode (RE). Here I used  a [Ag/AgCl/KCl (3 M)](https://www.basinc.com/products/MF-2056)
    - Your working electrode (WE). Here I used a a polished  [Glassy carbon electrode](https://www.basinc.com/products/MF-2012)   (GC, 3 mm diameter, BASI).

- Electrolyte
    - [K<sub>4</sub>[Fe(CN)<sub>6</sub>]](https://www.sigmaaldrich.com/ES/es/search/k4fe%5Bcn%5D6?focus=products&page=1&perPage=30&sort=relevance&term=K4Fe%5BCN%5D6&type=product_name)  0.004 M (0.000004 mol/cm<sup>3</sup>) dissolved in 0.1 M KCl (D = 7.6E−6 cm<sup>2</sup>/s)
   
# Electrochemical technique

Electrochemical techniques performed is cyclic voltammograms with a staircase profile (CV). Parameters: 

   - Start potential : 0.69 V.
   - Upper vertex potential: 0.7 V.
   - Lower vertex potential: -0.1 V.
   - Stop potential: 0.69 V. 
   - Number of scans: 4.
   - Step: 0.002 V.
   - Scan rates: 5 mV/seg , 10 mV/seg , 20 mV/seg , 40 mV/seg , 60 mV/seg ,80 mV/seg , 100 mV/seg , 125 mV/seg , 150 mV/seg , 200 mV/seg.
   
# Data treatment. Randles–Ševčík Analysis

Before you obtain results from data read this website for undestand Randles–Ševčík analysis:

https://en.wikipedia.org/wiki/Randles%E2%80%93Sevcik_equation

Once you understand this analysis. Following parameters were applied as Randles–Ševčík equation parameters:

- Faraday constant: 96485.33289 mol<sup>-1</sup>
- Diffusion coeficient for K<sub>4</sub>[Fe(CN)<sub>6</sub>] 0.0000073 cm<sup>2</sup>/s
- Concentration of K<sub>4</sub>[Fe(CN)<sub>6</sub>] 0.000004 mol/cm<sup>3</sup>
- Gas constant 8.3144598 J mol<sup>-1</sup> K<sup>-1</sup>
- Temperature 298 K



# Experimental NOVA software routine 

Once you finish setting-up all the electrochemical cell components. Create experimental sequence with the parameteres displayed in the  screenshots attached in this project.  

WARNING: Due to NOVA (Metrohm) software is licensed I can not upload sequence (.nox file). If you want NOVA (Metrohm) software sequence for ECSA measurements (.nox files) do no hesitate to contact me at :  gasparcarrascohuertas@gmail.com

# Bibliography

[A Practical Beginner’s Guide to Cyclic Voltammetry. Noémie Elgrishi, Kelley J. Rountree, Brian D. McCarthy, Eric S. Rountree, Thomas T. Eisenhart, and Jillian L. Dempsey
Journal of Chemical Education 2018 95 (2), 197-206. DOI: 10.1021/acs.jchemed.7b00361](https://pubs.acs.org/doi/10.1021/acs.jchemed.7b00361)
