---
title: "Notes"
author: "Richard Zhao"
date: "2024-12-05"
knit: (function(inputFile, encoding) {
  rmarkdown::render(inputFile, encoding = encoding, output_dir = "docs") })
bibliography: phenHet.bib
---



# Phenomenological Incidence with Network

We would like to investigate how to generate better interpretations or insights to nonlinear/phenomenological incidence $f(S,I)$(instead of $SI$) in compartmental epidemiology model, using graph/network models for the social structure/heterogeneity.

## Literature/Scoping Review

(TO DO) A good start point is to conduct a literature/scoping review.

### Papers/authors we know, but might biased.

1.  MGM. Gomes
    a. [MGM. Gomes,EtAl.(2022)](https://doi.org/10.1016/j.jtbi.2022.111063) : Individual variation in susceptibility or exposure to SARS-CoV-2 lowers the herd immunity threshold. (PMID: 35189135)
    b. Argument about COVID with S. Bansal
    c. More?
2.  G. Dwyer & J. Dushoff
    a.  [G. Dwyer, J.S. Elkinton & J.P. Buonaccorsi (1997)](https://doi.org/10.1086/286089): Host Heterogeneity in Susceptibility and Disease Dynamics: Tests of a Mathematical Model.(PMID: 18811331)
    b.  [G. Dwyer, J. Dushoff, J.S. Elkinton & S.A. Levin (2000)](https://doi.org/10.1086/303379): Pathogen‐Driven Outbreaks in Forest Defoliators Revisited: Building Models from Experimental Data. (PMID: 10856195)
    c.  More?
    d.  (TO DO) Email Chain (Todd's claim)
    e.  [JD Notes](http://dushoff.github.io/notebook/outputs/powerPhenHet.wt.math)
3.  M.J. Keeling & B.T. Grenfell (pair approx)
    a.  [M.J. Keeling (1999)](https://doi.org/10.1098/rspb.1999.0716): The Effects of Local Spatial Structure on Epidemiological Invasions(<PMID:10343409>)
    b.  [D.A. Rand (1999)](https://doi.org/10.1002/9781444311501.ch4): Correlation Equations and Pair Approximations for Spatial Ecologies(Book Chapters, no hit in both query yet)
    c.  [B. Finkenstädt, B.T. Grenfell (2002)](https://doi.org/10.1111/1467-9876.00187): Time series modelling of childhood diseases: A dynamical systems approach. (No hit on PubMed or MathSciNet)
4.  R.M. Granich
    a.  [R.M. Granich EtAl(2009)](https://doi.org/10.1016/S0140-6736(08)61697-9): Universal voluntary HIV testing with immediate antiretroviral therapy as a strategy for elimination of HIV transmission: a mathematical model(PMID: 19038438)
    b.  More?
5.  Power-law forms
    a.  [E.B. Wilson & J. Worcester(1944)](https://doi.org/10.1073/pnas.31.1.24): The Law of Mass Action in Epidemiology(PMID: 16588678)
    b.  [W. Liu, H.W. Hethcote & S.A. Levin (1987)](https://doi.org/10.1007/BF00277162): Dynamical behavior of epidemiological models with nonlinear incidence rates(PMID: 3668394, MR0908379)
    c.  [W. Liu, S.A. Levin & Y. Iwasa(1986)](https://doi.org/10.1007/BF00276956): Influence of nonlinear incidence rates upon the behavior of SIRS epidemiological models(PMID: 3958634, MR0829132)
6.  C. Rose
    a.  [C. Rose EtAl. sculpting paper](https://doi.org.libaccess.lib.mcmaster.ca/10.1016/j.jtbi.2021.110839): Heterogeneity in susceptibility dictates the order of epidemic models (PMID: 34314731, MR4297041)
    b.  Others?
7.  H. Berestycki
    a.  [H. Berestycki EtAl.](https://link.springer.com/article/10.1007/s00285-022-01861-w): Epidemic modeling with heterogeneity and social diffusion (PMID: 36964799, MR4568209)
    b.  Others?
8.  A. Korobeinikov
    a.  [A. Korobeinikov(2006)](https://link.springer.com/article/10.1007/s11538-005-9037-9): Lyapunov Functions and Global Stability for SIR and SIRS Epidemiological Models with Non-Linear Transmission (PMID: 16794947)
    b.  [A. Korobeinikov(2007)](https://link.springer.com/article/10.1007/s11538-007-9196-y): Global Properties of Infectious Disease Models with Nonlinear Incidence (PMID: 17443392, MR2224783)
9.  A.S. Novozhilov
    a.  [Novozhilov(2008)](https://www.sciencedirect.com/science/article/pii/S0025556408001211?via%3Dihub): On the spread of epidemics in a closed heterogeneous population(PMID: 18722386, MR2462419)

### A more unbiased query to literature

Potential Search Key:

-   Epidem\* & Infect\* & (Heterogeneity \| Structure \| Network) & (Nonlinear Incidence \| Dynamics)
    -   Might require finer tuning for different platform
    -   Need to filter the results to focus more on Incidence not R_0 or Final infection size.

#### Searching Platforms:

-   [Web of Science](https://www-webofscience-com.libaccess.lib.mcmaster.ca/wos/alldb/basic-search)
-   [Scopus](https://www.scopus.com/search/form.uri?display=basic#basic)
-   [Google Scholar](https://scholar.google.com/?hl=en&as_sdt=0,5)
-   [Open Alex](https://openalex.org/)
-   [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
-   [MathSciNet](https://mathscinet.ams.org/mathscinet/publications-search)

#### Searching Jornal

-   [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
	1.  [epidem\* AND (heterogeneity OR structure OR network)](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+%28heterogeneity+OR+structure+OR+network%29): 283,962
	    -  Include: not finished 10k limitation
	2.  [epidem\* AND infect\* AND (heterogeneity OR structure OR network)](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29): 70,839
	    - Include: not finished 10k limitation
	3.  [2+ AND dynamics](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynamics): 6,683
	    - Include: 2b
	4.  [2+ AND dynamics AND nonlinear](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynamics+AND+nonlinear): 276
	    -   [2+ AND dynam\* AND nonlinear](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+nonlinear)
	5.  [2+ AND dynam\* AND nonlinear AND incidence](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+nonlinear+AND+incidence):195
	    -   [2+ AND dynam\* AND nonlinear incidence](https://pubmed.ncbi.nlm.nih.gov/?term=epidem%2A+AND+infect%2A+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam%2A+AND+nonlinear+incidence&sort=relevance)
	6.  [2+ AND dynam\* AND "nonlinear incidence"](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+%22nonlinear+incidence%22): 17

-   **New Query: Trying to hit every selected paper first**

	7.  [(epidem\* \| disease \| infect\*) AND (dynam\* \| equilibri\*) AND (inciden\* \| susceptib\* \| spread \| transmission) AND (hetero\* \| network \| variation \| nonlinear \| non-linear)](https://pubmed.ncbi.nlm.nih.gov/?term=%28epidem*+%7C+disease+%7C+infect*%29+AND+%28dynam*+%7C+equilibri*%29+AND+%28inciden*+%7C+susceptib*+%7C+spread+%7C+transmission%29+AND+%28hetero*+%7C+network+%7C+variation+%7C+nonlinear+%7C+non-linear%29%22&size=200): 14015
	8.  [7 without network](https://pubmed.ncbi.nlm.nih.gov/?term=%28epidem*+OR+disease+OR+infect*%29+AND+%28dynam*+OR+equilibri*%29+AND+%28inciden*+OR+susceptib*+OR+spread+OR+transmission%29+AND+%28hetero*+OR+variation+OR+nonlinear+OR+non-linear%29&size=200):10218
	9.  [(epidem\* \| disease \| infect\*) AND (dynam\* \| equilibri\*) AND (inciden\* \| susceptib\* \| spread \| transmission) AND (hetero\* \| **(individual variation)** \| nonlinear \| non-linear)](https://pubmed.ncbi.nlm.nih.gov/?term=%28epidem*+OR+disease+OR+infect*%29+AND+%28dynam*+OR+equilibri*%29+AND+%28inciden*+OR+susceptib*+OR+spread+OR+transmission%29+AND+%28hetero*+OR+%28individual+variation%29+OR+nonlinear+OR+non-linear%29&size=200): 7656 with all hitting
	10. 9 & (power): 226, hit 6a, 9a
	11. 9 & (exponential): 132, hit 6a

-   [MathSciNet](https://mathscinet.ams.org/mathscinet/publications-search)

	1.  [epidem\* AND (heterogeneity OR structure OR network)](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29&page=1&size=20&sort=newest&facets=): 4,492
	2.  [epidem\* AND infect\* AND (heterogeneity OR structure OR network)](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29&page=1&size=20&sort=newest&facets=): 1,674
	3.  [2+ AND dynamics](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynamics&page=1&size=20&sort=newest&facets=): 978
	4.  [2+ AND dynam\* AND nonlinear](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20nonlinear&page=1&size=20&sort=newest&facets=): 238
	5.  [2+ AND dynam\* AND nonlinear incidence](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20nonlinear%20incidence&page=1&size=20&sort=newest&facets=):43
	6.  [2+ AND dynam\* AND "nonlinear incidence"](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20%22nonlinear%20incidence%22&page=2&size=20&sort=newest&facets=):25

Both easyPM and Rentrenz cannot download record or fetch ids more than 10,000 from pubmed/NCBI.

- [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) Unix based api is required
- Amy's [pubmed_notes](pubmed_notes.md)

#### By-path: Author Terms in query

```         
AND ((Gomes, Mgm[Author]) OR (Dwyer, G[Author]) OR (Dushoff, J[Author]) OR (Elkinton, Js[Author]) OR (Keeling, Mj[Author]) OR (Grenfell, Bt[Author]) OR (Granich, Rm[Author]) OR (Wilson, EB[Author]) OR (Worcester, J[Author]) OR (Hethcote, HW[Author]) OR (Levin, SA[Author]) OR (Liu, W[Author]) OR (Berestycki, H[Author]) OR (Rose, C[Author]) OR (Korobeinikov, A[Author]) OR (Novozhilov, AS[Author]))
```

-   [PubMed](https://pubmed.ncbi.nlm.nih.gov/)
	1.  [1+Authors](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+%28Gomes%2C+Mgm%5BAuthor%5D+OR+Dwyer%2C+G%5BAuthor%5D+OR+Dushoff%2C+J%5BAuthor%5D+OR+Keeling+MJ%5BAuthor%5D+OR+Grenfell%2C+Bt%5BAuthor%5D+OR+Granich%2C+Rm%5BAuthor%5D+OR+Wilson+EB%5BAuthor%5D+OR+Worcester+J%5BAuthor%5D+OR+Hethcote+HW%5BAuthor%5D+OR+Levin+SA%5BAuthor%5D+OR+Liu%2C+W%5BAuthor%5D%29&sort=): 619 results, Include 2b, 3b, 6a, 7a
	2.  [2+Authors](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+%28Gomes%2C+Mgm%5BAuthor%5D+OR+Dwyer%2C+G%5BAuthor%5D+OR+Dushoff%2C+J%5BAuthor%5D+OR+Keeling+MJ%5BAuthor%5D+OR+Grenfell%2C+Bt%5BAuthor%5D+OR+Granich%2C+Rm%5BAuthor%5D+OR+Wilson+EB%5BAuthor%5D+OR+Worcester+J%5BAuthor%5D+OR+Hethcote+HW%5BAuthor%5D+OR+Levin+SA%5BAuthor%5D+OR+Liu%2C+W%5BAuthor%5D%29&sort=): 309 results, Include 2b, 6a, 7a
	3.  [7+Authors](https://pubmed.ncbi.nlm.nih.gov/?term=%28epidem*+%7C+disease+%7C+infect*%29+AND+%28dynam*+%7C+equilibri*%29+AND+%28inciden*+%7C+susceptib*+%7C+spread+%7C+transmission%29+AND+%28hetero*+%7C+network+%7C+variation+%7C+nonlinear+%7C+non-linear%29+AND+%28%28Gomes%2C+Mgm%5BAuthor%5D%29+OR+%28Dwyer%2C+G%5BAuthor%5D%29+OR+%28Dushoff%2C+J%5BAuthor%5D%29+OR+%28Elkinton%2C+Js%5BAuthor%5D%29+OR+%28Keeling%2C+Mj%5BAuthor%5D%29+OR+%28Grenfell%2C+Bt%5BAuthor%5D%29+OR+%28Granich%2C+Rm%5BAuthor%5D%29+OR+%28Wilson%2C+EB%5BAuthor%5D%29+OR+%28Worcester%2C+J%5BAuthor%5D%29+OR+%28Hethcote%2C+HW%5BAuthor%5D%29+OR+%28Levin%2C+SA%5BAuthor%5D%29+OR+%28Liu%2C+W%5BAuthor%5D%29+OR+%28Berestycki%2C+H%5BAuthor%5D%29+OR+%28Rose%2C+C%5BAuthor%5D%29+OR+%28Korobeinikov%2C+A%5BAuthor%5D%29%29&sort=): 198 results (14015 without authors), Include 1a, 2a, 2b, 3a, 4a, 5b, 5c, 6a, 7a, 8a, 8b, 9a
        -   Rentrez result: almost hit every selected paper
        -   3b is a book chapter with no record, 5a is too old with no abstract for hitting
    2.  [8+Authors]():164 results (10218 without authors), Include 1a, 2a, 2b, 3a, 4a, 5b, 5c, 6a, 7a, 8a, 8b, 9a

-   [MathSciNet](https://mathscinet.ams.org/mathscinet/publications-search)

Idea: Specify fields to better refine the result (TBD)

### Checking the references:

Remove duplicates in two journals and check if the papers interested is in each result of combinations of keywords.

### Literature next steps: 
- (TO DO) MeSH Field
- (TO DO) Take a random sample from (e.g. query #9) and check the quality of the query
- (TO DO) Summarize the query journal more properly
- ?? Build a citation network from core references listed in this note
- ?? Screen a random sample 
	- 300-500 (with AI)
	- 100 by ourselves
- (TO DO)Refine the query more

## Next Steps

### ==Idea: A table of heterogeneity assumptions==

Considering how heterogeneity could affect incidence: 
- Direct: Susceptibility $\times$ Contact Rate $\times$ Infectivity.
- Indirect: Recovery (affect the duration of infectiousness), Recalcitrance(Immunity or VD)

| Susceptibility    | Contact Rate                                                | Infectivity                                                     | Recovery    | Recalcitrance (Immunity or VD) | Incidence Term      | Literature |
| ----------------- | ----------------------------------------------------------- | --------------------------------------------------------------- | ----------- | ------------------------------ | ------------------- | ---------- |
| Gamma Distrib.    | Fully Mix                                                   | Hom                                                     | Hom |                                | $(S/S_0)^p I$       | JD & Todd  |
| Parametric Het | Fully Mix                                                   | Hom                                                     | Hom |                                | $f(S/S_0) I$        | Novozhilov |
| Parametric Het | Fully Mix                                                   | Parametric Het, Perfect inheriting from infector to infectee | Hom |                                | $f(S/S_0) g(I/I_0)$ | Novozhilov |
|                   |                                                             |                                                                 |             |                                |                     |            |
| Hom       | Degree distribution, Interchangeable nodes from same compartment  | Hom                                                     |             |                                |                     | MFSH Network    |
|                   |                                                             |                                                                 |             |                                |                     |            |
|                   |                                                             |                                                                 |             |                                |                     |            |

[I.Z. Kiss, J.C. Miller & P.L. Simon(2017)](https://link.springer.com/book/10.1007/978-3-319-50806-1)
For Configuration Network model, assuming susceptibility and infectivity is independent: 
- Increase heterogeneity in infectivity decrease the Epidemic probability but will not affect the attack rate/final infection size. 
- Increase heterogeneity in susceptibility decrease the the attack rate but will not affect the epidemic probability.

### Empirical Study on Network

Start with empirical networks from simulation/data and try to investigate the expression of incidence

### Constructional Study

Start with certain incidence expression or certain network structure and try to construct the corresponding network. (?? such that the result would be numerically close enough)

## Miscellaneous Ideas

- Just came across [tuschhoffHeterogeneityCorrelationHost2024](https://www.medrxiv.org/content/10.1101/2024.12.10.24318805v1) in Google's "recommended papers" alert.
- [H. Berestycki EtAl.](https://link.springer.com/article/10.1007/s00285-022-01861-w) paper (Todd & Jonathan)
- [MGM. Gomes,EtAl.(2022)](https://doi.org/10.1016/j.jtbi.2022.111063) paper (and the argument with Bansal)
- [JD's Talk about heterogeneity]([https://www.youtube.com/watch?v=9OhB3WfSpS8](https://www.youtube.com/watch?v=9OhB3WfSpS8 "https://www.youtube.com/watch?v=9OhB3WfSpS8")) and Thesis
- [I.Z. Kiss, J.C. Miller & P.L. Simon(2017)](https://link.springer.com/book/10.1007/978-3-319-50806-1)'s Claim about heterogeneity of transmission and recovery on Network.
-  Construction of possible $f(S,I/I_0)$ incidence term.
- [M.C. Bootsma, K.M.D. Chan, O. Diekmann, H. Inaba(2024)](https://doi.org/10.5206/mase/16718)


## TO DO List

- [H. Berestycki EtAl.](https://link.springer.com/article/10.1007/s00285-022-01861-w) paper (Todd & Jonathan)
- [MGM. Gomes,EtAl.(2022)](https://doi.org/10.1016/j.jtbi.2022.111063) paper
- Finish [M.C. Bootsma, K.M.D. Chan, O. Diekmann, H. Inaba(2024)](https://doi.org/10.5206/mase/16718) 
- Exponential incidence term for $e^S$ in Novozhilov framework
- Effective contact/incidence rate $(0,1)$: [final-size repo](https://github.com/toddlparsons/final-size)
- Tails truncation

- Long term: empirical network and incidence term
## References
