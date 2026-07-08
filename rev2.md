This document contains my (the author) second revision of the paper. This review is mostly concentrated on the writing clarity and compactness of the arguments presented. I wish for the following edits to be made to improve the overall quality of the paper to be of quality publishable in a top journal such as Automatica or IEEE TAC. The following points are my main concerns:
:
1. **Reducing repetition and redundancy**: The paper should be more concise and clearly written without unnecessary repetition of ideas or phrases. Each argument should be presented only once, and any repeated information should be removed (or kept if the context really requires it).

2. **Reducing the side ideas**: The paper should focus on the main arguments and avoid digressing into side ideas that do not directly support the thesis. Any tangential discussions should be minimized or removed to maintain a clear and focused narrative.

3. **Some ideas to reduce their size in the paper**: The following ideas should be condensed or summarized to reduce their size in the paper:
   - The linearization of the cost and the f.o.t. discussion in the MPC-eKF case can be avoided and the cost immediately assumed quadratic. In a sentence or more, we can mention that this quadratic cost can come from local linearization.
   - The \ell_1 discussion together with its epigraph approximation (for smoothness and convergence of solvers) can be reduced into a single paragraph (or more) as a small subsection.
   - The paper keep comparing each example gradients to finite-difference gradients. However, the gradients we find are exact and do not need to be verified. You can mention this verification only at one point in the paper and remove the repeated comparisons.

4. **The presented numerical examples should be self-sufficient**: Pull out the necessary information from the Appendix and put them in the main text, presenting each example and how we approached it concisely. Then delete the numerical example details in the Appendix.

5. **The code for reproducibility**: The paper should mention that the code for reproducing the results is available in a public repository (e.g., GitHub) and provide a link to it (you can use a placeholder for the link for now).
