---
title: "Notes"
author: "Richard Zhao"
date: "2024-12-05"
knit: (function(inputFile, encoding) {
  rmarkdown::render(inputFile, encoding = encoding, output_dir = "docs") })
output:
  html_document:
    keep_md: true
bibliography: phenHet.bib
---



# Phenomenological Incidence with Network

We would like to investigate how to generate better interpretations or insights to nonlinear/phenomenological incidence $f(S,I)$(instead of $SI$) in compartmental epidemiology model, using graph/network models for the social structure/heterogeneity.

## Literature/Scoping Review

(TO DO) A good start point is to conduct a literature/scoping review.

### Papers/authors we know, but might biased.

1.  MGM. Gomes
    a.  [MGM. Gomes,EtAl.(2022)](https://doi.org/10.1016/j.jtbi.2022.111063) : Individual variation in susceptibility or exposure to SARS-CoV-2 lowers the herd immunity threshold. (PMID: 32511451)
    b.  (??) Argument about COVID with S. Bansal
    c.  More?
2.  G. Dwyer & J. Dushoff
    a.  [G. Dwyer, J.S. Elkinton & J.P. Buonaccorsi (1997)](https://doi.org/10.1086/286089): Host Heterogeneity in Susceptibility and Disease Dynamics: Tests of a Mathematical Model.(PMID: 18811331)
    b.  [G. Dwyer, J. Dushoff, J.S. Elkinton & S.A. Levin (2000)](https://doi.org/10.1086/303379): Pathogen‐Driven Outbreaks in Forest Defoliators Revisited: Building Models from Experimental Data. (PMID: 10856195)
    c.  More?
    d.  (TO DO) Email Chain (Todd's claim)
    e.  [JD Notes](http://dushoff.github.io/notebook/outputs/powerPhenHet.wt.math)
3.  M.J. Keeling & B.T. Grenfell (pair approx)
    a.  [M.J. Keeling (1999)](https://doi.org/10.1098/rspb.1999.0716): The Effects of Local Spatial Structure on Epidemiological Invasions(<PMID:10343409>)
    b.  [D.A. Rand (1999)](https://doi.org/10.1002/9781444311501.ch4): Correlation Equations and Pair Approximations for Spatial Ecologies(Book Chapters, no hit in both query yet)
    c.  More?
4.  R.M. Granich
    a.  [R.M. Granich EtAl(2009)](https://doi.org/10.1016/S0140-6736(08)61697-9): Universal voluntary HIV testing with immediate antiretroviral therapy as a strategy for elimination of HIV transmission: a mathematical model(PMID: 19038438)
    b.  More?
5.  Power-law forms
    a.  [E.B. Wilson & J. Worcester(1944)](https://doi.org/10.1073/pnas.31.1.24): The Law of Mass Action in Epidemiology(PMID: 16588678)
    b.  [W. Liu, H.W. Hethcote & S.A. Levin (1987)](https://doi.org/10.1007/BF00277162): Dynamical behavior of epidemiological models with nonlinear incidence rates(PMID: 3668394, MR0908379)
    c.  [W. Liu, S.A. Levin & Y. Iwasa(1986)](https://doi.org/10.1007/BF00276956): Influence of nonlinear incidence rates upon the behavior of SIRS epidemiological models(PMID: 3958634, MR0829132)
6.  C. Rose
    a.  [C. Rose EtAl. sculpting paper](https://doi.org.libaccess.lib.mcmaster.ca/10.1016/j.jtbi.2021.110839)
7.  H. Berestycki
    a.  [H. Berestycki EtAl.](https://link.springer.com/article/10.1007/s00285-022-01861-w)

## \### Others to consider adding above

### A more unbiased query to literature

Potential Search Key:

-   Epidem\* & Infect\* & (Heterogeneity \| Structure \| Network) & (Nonlinear Incidence \| Dynamics )
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
    -   Include: not finished 10k limitation
2.  [epidem\* AND infect\* AND (heterogeneity OR structure OR network)](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29): 70,839
    -   Include: not finished 10k limitation
3.  [2+ AND dynamics](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynamics): 6,683
    -   Include: 2b
4.  [2+ AND dynamics AND nonlinear](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynamics+AND+nonlinear): 276
    -   [2+ AND dynam\* AND nonlinear](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+nonlinear)
5.  [2+ AND dynam\* AND nonlinear AND incidence](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+nonlinear+AND+incidence):195
    -   [2+ AND dynam\* AND nonlinear incidence](https://pubmed.ncbi.nlm.nih.gov/?term=epidem%2A+AND+infect%2A+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam%2A+AND+nonlinear+incidence&sort=relevance)
6.  [2+ AND dynam\* AND "nonlinear incidence"](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+%22nonlinear+incidence%22): 17

-   [MathSciNet](https://mathscinet.ams.org/mathscinet/publications-search)

1.  [epidem\* AND (heterogeneity OR structure OR network)](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29&page=1&size=20&sort=newest&facets=): 4,492
2.  [epidem\* AND infect\* AND (heterogeneity OR structure OR network)](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29&page=1&size=20&sort=newest&facets=): 1,674
3.  [2+ AND dynamics](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynamics&page=1&size=20&sort=newest&facets=): 978
4.  [2+ AND dynam\* AND nonlinear](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20nonlinear&page=1&size=20&sort=newest&facets=): 238
5.  [2+ AND dynam\* AND nonlinear incidence](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20nonlinear%20incidence&page=1&size=20&sort=newest&facets=):43
6.  [2+ AND dynam\* AND "nonlinear incidence"](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20%22nonlinear%20incidence%22&page=2&size=20&sort=newest&facets=):25

Both easyPM and Rentrenz cannot download record or fetch ids more than 10,000 from pubmed/NCBI.

-   (TODO?) [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) Unix based api is requried

### By-path: Author Terms in query

```         
AND ((Gomes, Mgm[Author]) OR (Dwyer, G[Author]) OR (Dushoff, J[Author]) OR (Elkinton, Js[Author]) OR (Keeling, Mj[Author]) OR (Grenfell, Bt[Author]) OR (Granich, Rm[Author]) OR (Wilson, EB[Author]) OR (Worcester, J[Author]) OR (Hethcote, HW[Author]) OR (Levin, SA[Author]) OR (Liu, W[Author]))`
```

-   [PubMed](https://pubmed.ncbi.nlm.nih.gov/)

    1.  [1+Authors](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+%28Gomes%2C+Mgm%5BAuthor%5D+OR+Dwyer%2C+G%5BAuthor%5D+OR+Dushoff%2C+J%5BAuthor%5D+OR+Keeling+MJ%5BAuthor%5D+OR+Grenfell%2C+Bt%5BAuthor%5D+OR+Granich%2C+Rm%5BAuthor%5D+OR+Wilson+EB%5BAuthor%5D+OR+Worcester+J%5BAuthor%5D+OR+Hethcote+HW%5BAuthor%5D+OR+Levin+SA%5BAuthor%5D+OR+Liu%2C+W%5BAuthor%5D%29&sort=): 619 results, Include 2b, 3b
    2.  [2+Authors](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+%28Gomes%2C+Mgm%5BAuthor%5D+OR+Dwyer%2C+G%5BAuthor%5D+OR+Dushoff%2C+J%5BAuthor%5D+OR+Keeling+MJ%5BAuthor%5D+OR+Grenfell%2C+Bt%5BAuthor%5D+OR+Granich%2C+Rm%5BAuthor%5D+OR+Wilson+EB%5BAuthor%5D+OR+Worcester+J%5BAuthor%5D+OR+Hethcote+HW%5BAuthor%5D+OR+Levin+SA%5BAuthor%5D+OR+Liu%2C+W%5BAuthor%5D%29&sort=): 309 results, Include 2b

-   [MathSciNet](https://mathscinet.ams.org/mathscinet/publications-search)

### Checking the references:

Remove duplicates in two journals and check if the papers interested is in each result of combinations of keywords.

## Next Steps

### Empirical Study

Start with empirical/simulated networks from data and try to investigate the expression of incidence

### Constructional Study

Start with certain incidence expression or certain network structure and try to construct the corresponding network. (?? such that the result would be numerically close enough)

## Miscellaneous

Just came across @tuschhoffHeterogeneityCorrelationHost2024 in Google's "recommended papers" alert

## References
