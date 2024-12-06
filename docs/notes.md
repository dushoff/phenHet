---
title: "Note"
author: "Richard Zhao"
date: "2024-12-05"
knit: (function(inputFile, encoding) {
  rmarkdown::render(inputFile, encoding = encoding, output_dir = "docs") })
output:
  html_document:
    keep_md: true
---



# Phenomenological Incidence with Network

We would like to investigate how to generate better interpretations or insights to nonlinear/phenomenological incidence $f(S,I)$(instead of $SI$) in compartmental epidemiology model, using graph/network models for the social structure/heterogeneity.

## Literature/Scoping Review

(TO DO) A good start point is to conduct a literature/scoping review.

### Papers/authors we know, but might biased.

-   MGM. Gomes
    -   [MGM. Gomes,EtAl.(2022)](https://doi.org/10.1016/j.jtbi.2022.111063) : Individual variation in susceptibility or exposure to SARS-CoV-2 lowers the herd immunity threshold.
    -   (??) Argument about COVID with S. Bansal
    -   More?
-   G. Dwyer & J. Dushoff
    -   [G. Dwyer, J.S. Elkinton & J.P. Buonaccorsi (1997)](https://doi.org/10.1086/286089): Host Heterogeneity in Susceptibility and Disease Dynamics: Tests of a Mathematical Model.
    -   [G. Dwyer, J. Dushoff, J.S. Elkinton & S.A. Levin (2000)](https://doi.org/10.1086/303379): Pathogen‐Driven Outbreaks in Forest Defoliators Revisited: Building Models from Experimental Data.”
    -   More?
    -   (TO DO) Email Chain (Todd's claim)
    -   [JD Notes](http://dushoff.github.io/notebook/outputs/powerPhenHet.wt.math)
-   M.J. Keeling & B.T. Grenfell (pair approx)
    -   (??) Any specific paper, especially for Grenfell
    -   [M.J. Keeling (1999)](https://doi.org/10.1098/rspb.1999.0716): The Effects of Local Spatial Structure on Epidemiological Invasions
    -   [D.A. Rand (1999)](https://doi.org/10.1002/9781444311501.ch4): Correlation Equations and Pair Approximations for Spatial Ecologies
    -   More?
-   R.M. Granich
    -   [R.M. Granich EtAl(2009)](https://doi.org/10.1016/S0140-6736(08)61697-9): Universal voluntary HIV testing with immediate antiretroviral therapy as a strategy for elimination of HIV transmission: a mathematical model'
    -   More?
-   Power-law forms
    -   [E.B. Wilson & J. Worcester(1944)](https://doi.org/10.1073/pnas.31.1.24): The Law of Mass Action in Epidemiology
    -   [W. Liu, H.W. Hethcote & S.A. Levin (1987)](https://doi.org/10.1007/BF00277162): Dynamical behavior of epidemiological models with nonlinear incidence rates
    -   [W. Liu, S.A. Levin & Y. Iwasa(1986)](https://doi.org/10.1007/BF00276956): Influence of nonlinear incidence rates upon the behavior of SIRS epidemiological models

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
-   [mathscinet](https://mathscinet.ams.org/mathscinet/publications-search)

#### Searching Jornal

-   [PubMed](https://pubmed.ncbi.nlm.nih.gov/)

1.  [epidem\* AND (heterogeneity OR structure OR network)](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+%28heterogeneity+OR+structure+OR+network%29): 283,962
2.  [epidem\* AND infect\* AND (heterogeneity OR structure OR network)](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29): 70,839
3.  [2+ AND dynamics](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynamics): 6,683
4.  [2+ AND dynamics AND nonlinear](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynamics+AND+nonlinear): 276
    -   [2+ AND dynam\* AND nonlinear](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+nonlinear)
5.  [2+ AND dynam\* AND nonlinear AND incidence](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+nonlinear+AND+incidence):195
    -   [2+ AND dynam\* AND nonlinear incidence](https://pubmed.ncbi.nlm.nih.gov/?term=epidem%2A+AND+infect%2A+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam%2A+AND+nonlinear+incidence&sort=relevance)
6.  [2+ AND dynam\* AND "nonlinear incidence"](https://pubmed.ncbi.nlm.nih.gov/?term=epidem*+AND+infect*+AND+%28heterogeneity+OR+structure+OR+network%29+AND+dynam*+AND+%22nonlinear+incidence%22): 17

-   [mathscinet](https://mathscinet.ams.org/mathscinet/publications-search)

1.  [epidem\* AND (heterogeneity OR structure OR network)](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29&page=1&size=20&sort=newest&facets=): 4,492
2.  [epidem\* AND infect\* AND (heterogeneity OR structure OR network)](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29&page=1&size=20&sort=newest&facets=): 1,674
3.  [2+ AND dynamics](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynamics&page=1&size=20&sort=newest&facets=): 978
4.  [2+ AND dynam\* AND nonlinear](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20nonlinear&page=1&size=20&sort=newest&facets=): 238
5.  [2+ AND dynam\* AND nonlinear incidence](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20nonlinear%20incidence&page=1&size=20&sort=newest&facets=):43
6.  [2+ AND dynam\* AND "nonlinear incidence"](https://mathscinet.ams.org/mathscinet/publications-search?query=epidem%2a%20AND%20infect%2a%20AND%20%28heterogeneity%20OR%20structure%20OR%20network%29%20AND%20dynam%2a%20AND%20%22nonlinear%20incidence%22&page=2&size=20&sort=newest&facets=):25

### Cheking the references:
Remove duplicates in two journals and check if the papers interested is in each result of combinations of keywords.


## Next Steps

### Empirical Study

Start with empirical/simulated networks from data and try to investigate the expression of incidence

### Constructional Study

Start with certain incidence expression or certain network structure and try to construct the corresponding network. (?? such that the result would be numerically close enough)
