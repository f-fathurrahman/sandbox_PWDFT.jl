# # Rational numbers
#
# In julia rational numbers can be constructed with the `//` operator.
# Lets define two rational numbers, `x` and `y`:

## Define variable x and y
x = 1//3;
y = 2//5;

# When adding `x` and `y` together we obtain a new rational number:

z = x + y;

println("This is z = $z")

#=
This is some markdown with equation
$$
  \alpha + \beta = \Gamma
$$
=#

# And this is another calculation
xx = sin(z)

# ## Example plot
using Plots, PlotThemes
theme(:dark)
x = range(0, stop = 6π, length = 1000)
y1 = sin.(2*x)
y2 = cos.(x)
plot(x, [y1, y2], fmt=:svg, size=(400,300))

# ## Another plot
#=
This is another plot.
=#
plot(fmt=:svg, size=(400,300))
histogram!(randn(10_000))
