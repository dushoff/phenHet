[resources.md](./resources.md): Summary notes for related literatures and previous investigations:

[approaches.md](approaches.md): A quick summary document 

[JR_NegBinom_Result.md](JR_NegBinom_Result.md): The main background and proceeding document for current MSV random network. RZ is working on summarize recent progress to update this document.   
- Compiled as html here: https://dushoff.github.io/phenHet/JR_NegBinom_Result.html (Hyperlink in .html version not working due to folder issue)



## simulator info

* `NetSimulator.cpp` is a vertex-based simulator 
* `edgelist.cpp` is an edge-based simulator
* `scaleExamples.R` is set up to use *either* machinery, based on what's in the Makefile: compare `scaleEdges.Rout` with `scaleExamples.Rout`. `compareSims.R` (not yet plumbed) is the beginning of machinery to compare the two approaches.



