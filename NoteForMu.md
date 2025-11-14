# 1. Class and Probability of Edges
### 1.1 Whole network
Let us ignore the vertices and their status for a moment, and focus on edges and their binary status, which means whether the edge has t
transmitted the infection or not. 

For simplicity later, I use the terminology for occupation status from percolation theory, at any time $t$:
- If the transmission has already passed through the edge (i.e. the edge has already transmitted the infection), the edge is **occupied**.
- If the transmission has not yet passed through the edge (i.e. the edge has not yet transmitted the infection), the edge is **unoccupied**. 
For our SIR model with no reinfection, an occupied edge can not switch back to unoccupied.

At any moment $t$, MSV framework classify all edges in the whole network (with proportion/probability $=1$) into two large classes using its core variable $\phi(t)$:
- $\phi(t)$ proportion of edges has not yet transmitted the infection (unoccupied)
- $1-\phi(t)$ proportion of edges has already transmitted the infection (occupied)
As we are modeling **random** networks, the proportions for edges can also be treated as probabilities:
- $\phi(t)=\mathbb{P}(\text{A randomly selected edge is occupied at } t)$ 
- $1-\phi(t)=\mathbb{P}(\text{A randomly selected edge is unoccupied at } t)$

### 1.2 A Random Focal Vertex $a$
Lets uniformly randomly pick a focal vertex $a$, **regardless of the its status**, with a random variable $K$ as its degree. For this random vertex $a$ in random network, choose a random edge $E$ connected to vertex $a$ is the same as choose a random edge from the whole network.

So at any time $t$, this edge $E$ should also have probability $\phi(t)$ being unoccupied and probability $1-\phi(t)$ being occupied.

### 1.3 A Neighbour Vertex $b$ of $a$, with unoccupied edge $E$
For the **unoccupied** probability $\phi(t)$ of $E$, we can further parse it finer based on the status (S, I or R) of neighbour vertex $b$ connect to $E$ on the other side than $a$.

For this neighbour vertex $b$, as we know it is connected to an unoccupied edge $E$, so it is not the same as a random chosen vertex from whole network.
Therefore, the probability that neighbour vertex $b$ is susceptible, infected and recovered respectfully should not be $S(t), I(t), R(t)$, i.e. the proportion of each compartment.

### 1.4 Two perspective
==I think here is the key issue and confusion point: ==
The vertex $b$ has two equivalent identities in random network, but looks different from the perspective of vertex or edge: 
1. $b$ is a random vertex connected to an unoccupied edge $E$. ($a$ is not involved).
2. From randomly chosen focal vertex $a$, $E$ is a random unoccupied edge connected to $a$, then $b$ is the neighbour connect to $a$ through $E$.

NOTE: In probabilities, AND is denoted by $\wedge$ and condition is denoted by $|$.

MSV framework use $\phi_S(t), \phi_I(t), \phi_R(t)$ to denote the probability that a vertex connected to unoccupied edge $E$ **AND** being susceptible, infected or recovered.
So based on the two perspective, we have two event with the same probability.
$$\begin{align}
\phi_S(t) = & \ \mathbb{P}(\text{A random vertex } b \text{ is susceptible  at time } t \wedge b \text{ attached to a unoccupied edge } E )
\\
= & \ \mathbb{P}(\text{For random vertex } a \text{, one of its random edge } E \text{ is unoccupied } 
\\ 
& \wedge E \text{ attached to a suseptible neighbour } b)
\end{align}$$Similar for $\phi_I$ and $\phi_R$.

Use the same language for $\phi(t)$, we get the partition:
$$\begin{align}
\phi(t) &= \mathbb{P}(\text{A random edge } E \text{ is unoccupied at } t )
\\
& = \mathbb{P}(\text{For random vertex } a \text{, one of its random edge } E \text{ is unoccupied})
\\
& = \phi_S(t) + \phi_I(t) + \phi_R(t)
\end{align}$$
And further we have the full partition of all edges in the network:
$$\\
1 =  (1-\phi(t)) + \phi(t) = (1-\phi(t)) + (\phi_S(t) + \phi_I(t) + \phi_R(t))
$$
### 1.5 Definition of $\phi_S$
The expression for $\phi_S$ is derived from the definition without $a$ involved:
$$\begin{align}
\phi_S(t) & =  \mathbb{P}(\text{A random vertex } b \text{ is susceptible  at time } t \wedge b \text{ attached to a unoccupied edge } E )
\\
& = \frac{G'_p(\phi(t))}{\delta}=\frac{G'_p(\phi(t))}{G'_p(1)}
\\
& =  \sum_k \frac{p_k k}{\delta} \phi(t)^{k-1}=\frac{1}{\delta}\sum_k k (p_k \phi^{k-1})
\\
& =\frac{1}{\delta}\sum_k k (\frac{ \phi \times p_k \times\phi^{k-1}}{\phi})
\\
& = \frac{1}{\delta}\sum_k k (\frac{\mathbb{P}(E \text{ is unoccupied} \wedge b \text{ has degree k }  \wedge b \text{ has no other occupied edge}) }{\mathbb{P}(E \text{ is unoccupied})})
\end{align}$$
So indeed as Todd point out, this is conditioned on the focal edge $E$ is unoccupied, but only from $b$ and $E$'s perspective.


