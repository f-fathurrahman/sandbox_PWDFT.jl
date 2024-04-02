# ElkDFTWrapper

A simple wrapper to Elk DFT package.

Mainly used for writing variables or arrays to Julia's serialized data.

Some convention (hopefully I consistently used this):

- Elk's original subroutines are all prefixed with function `elk_` followed with
  the subroutine name.

- "Getter" functions are not exported.
  They must be called with `ElkDFTWrapper.getter_function`.

- No specific convention for other function names.

