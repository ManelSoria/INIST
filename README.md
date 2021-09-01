# INIST

This module aims to provide the access to the properties of certain species in low ranges of temperature that NASA 
polynomials can not reach (these polynomials could be acceded with [HGS](https://github.com/ManelSoria/HGS.git)).
The original script was developed and is maintained by [M. Soria](https://directori.upc.edu/directori/dadesPersona.jsp?id=1003031),
L. Frezza and [C. Fuster](https://www.linkedin.com/in/caleb-fuster-b6a566182/). 
The database is downloaded from the [NIST](https://webbook.nist.gov/chemistry/fluid/)

There is no more documentation in this module that the comments on the functions, but you can check the master thesis developed with the Matlab code for more info.
Request it to author of this package or download it from [UPCcommons](https://upcommons.upc.edu/discover?filtertype_1=author&filter_relational_operator_1=contains&filter_1=Caleb+Fuster&submit_apply_filter=)
when the institution publishes it.
___
## INIST Function

The function provides the access a multiples properties depending on the inputs. 
All the properties returned are per kilo (kg).


<details> 
 <summary> help(INIST) </summary>

    INIST - From Matlab INIST ((c) Manel Soria, Caleb Fuster, Lorenzo Frezza)
    Interpolation of Nonideal Idiosyncratic Splendiferous Tables
    (c) Caleb Fuster
    Data downloaded from NIST web page
    2021

    The return of this function is a dictionary with the required property
    ("t_ps" and "t_ph" can return Quality further than T if it is saturated ).

    Units: T(K), p(bar), h and u: kJ/kg, v: m^3/kg, rho: kg/m^3 s: kJ/kgK,
    a: m/s, cv and cp: kJ/kgK, JT: bar/K, mu: Pa.s, k: W/mK, MM: kg/mol
    SF: N.m
    1st argument: substance name
                  'Database' to return the list of database elements
    2nd and remaining arguments:
     critical temperature     'tcrit'
     critical pressure        'pcrit'
     critical volume          'vcrit'
     molecular mass           'MM'
     saturation temperature   'tsat_p', p
     saturation pressure      'psat_t', T
     saturated liquid properties as a function of pressure
        volume          'vl_p' , p
        energy          'ul_p' , p
        enthalpy        'hl_p' , p
        entropy         'sl_p' , p
        specific heat coeff at constant volume:
                        'cvl_p', p
        specific heat coeff at constant pressure:
                        'cpl_p', p
        sound speed     'al_p' , p
        viscosity       'mul_p', p
        density         'rl_p' , p
        conductivity    'kl_p' , p
     saturated vapour properties as a function of pressure
        volume          'vv_p' , p
        energy          'uv_p' , p
        enthalpy        'hv_p' , p
        entropy         'sv_p' , p
        specific heat coeff at constant volume:
                        'cvv_p', p
        specific heat coeff at constant pressure:
                        'cpv_p', p
        sound speed     'av_p' , p
        viscosity       'muv_p', p
        density         'rv_p' , p
        conductivity    'kv_p' , p
     saturated liquid properties as a function of temperature
        volume          'vl_t' , t
        energy          'ul_t' , t
        enthalpy        'hl_t' , t
        entropy         'sl_t' , t
        specific heat coeff at constant volume:
                        'cvl_t', t
        specific heat coeff at constant pressure:
                        'cpl_t', t
        sound speed     'al_t' , t
        viscosity       'mul_t', t
        density         'rl_t' , t
        conductivity    'kl_t' , t
     saturated vapour properties as a function of temperature
        volume          'vv_t' , t
        energy          'uv_t' , t
        enthalpy        'hv_t' , t
        entropy         'sv_t' , t
        specific heat coeff at constant volume:
                        'cvv_t', t
        specific heat coeff at constant pressure:
                        'cpv_t', t
        sound speed     'av_t' , t
        viscosity       'muv_t', t
        density         'rv_t' , t
        conductivity    'kv_t' , t
     non-saturated properties as a function of pressure and temperature
        volume          'v_pt' , p , t
        energy          'u_pt' , p , t
        enthaply        'h_pt' , p , t
        entrophy        's_pt' , p , t
        specific heat coeff at constant volume:
                        'cv_pt', p , t
        specific heat coeff at constant pressure:
                        'cp_pt', p , t
        sound speed     'a_pt' , p , t
        viscosity       'mu_pt', p , t
        density         'r_pt' , p , t
        conductivity    'k_pt',  p , t
     temperature as a function of ...
        pressure and entropy 't_ps', p ,s

     special functions:
          'minp'        returns the minimum isobar available
          'maxp'        idem max isobar
          'mint'        idem minimum temperature
          'maxt'        idem maximum temperature
          'isobars'     returns a vector with the available isobars
</details>

---
## INIST Database

The default database is conformed by the following highlighted species with their isobars between 0.001 up to 700 bars (depend on the species limits, use *<span style="color:GreenYellow">INIST("species", "isobars")</span>* to know the limits).
(Due to size of the database becomes slow to download, we restrict it)


<details>
<summary> Species </summary>
    <ul>
        <li> <span style="color:red"><i>O<sub>2</sub></i></span>
        <li> <span style="color:red"><i>H<sub>2</sub></i></span>
        <li> H<sub>2</sub>O
        <li> <span style="color:red"><i>R-134a</i></span> (R134a in the database)
        <li> CO
        <li> <span style="color:red">CO<sub>2</sub></span>
        <li> <span style="color:red"><i>N<sub>2</sub></i></span>
        <li> <span style="color:red">NH<sub>3</sub></span>
        <li> N<sub>2</sub>O
        <li> <span style="color:red">He</span>
        <li> <details>
            <summary> Hydrocarbons </summary>
                <ul>
                    <li> <span style="color:red">CH<sub>4</sub></span>
                    <li> CH<sub>4</sub>O
                    <li> C<sub>2</sub>H<sub>4</sub>
                    <li> C<sub>2</sub>H<sub>6</sub>
                    <li> C<sub>3</sub>H<sub>4</sub>
                    <li> C<sub>3</sub>H<sub>6</sub>
                    <li> <span style="color:red">C<sub>3</sub>H<sub>8</sub></span>
                    <li> <span style="color:red">C<sub>4</sub>H<sub>10</sub></span>
                    <li> C<sub>5</sub>H<sub>12</sub>
                    <li> C<sub>6</sub>H<sub>6</sub> (Benzene)
                    <li> C<sub>6</sub>H<sub>12</sub>
                    <li> C<sub>6</sub>H<sub>14</sub>
                </ul>
        </details>
    </ul>
</details>

This module contains a function (*INISTdatabase*) to download new data from all the previous species from the NIST page. It is not recommended downloading all the database due to it requires 
a couple of hours but if you require them only run the INISTdatabase without the restrictions:


If you need other species, you can adapt the previous function. Species idcas can be found in this [ZIP](https://webbook.nist.gov/chemistry/download/species.zip) provided by NIST. 
However, you can obtain it [here](https://webbook.nist.gov/chemistry/fluid/) following the next steps:


<ol>
    <li> Select the species</li>
    <li> Press Continue and put numbers in the range that the page request you</li>
    <li> Go to the URL and search ID=*& where asterisk is the letters and numbers that correspond to the idcas</li>
    <li> In this page, where Auxiliary Data, it can be found the reference to the enthalpy and entropy</li>
</ol>

